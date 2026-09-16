from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import h5py
import numpy as np
import pandas as pd
from scipy.io import loadmat, whosmat

from .config import ensure_config
from .provenance import file_hash, digest

SCAN_ENTITIES = ("task", "acq", "run", "echo", "rec", "dir")


def subject_code(value):
    value = str(value).strip().removeprefix("sub-").upper()
    if re.fullmatch(r"\d+_S_\d+", value):
        value = "ADNI" + value.replace("_", "")
    if not re.fullmatch(r"ADNI\d+S\d+", value):
        raise ValueError(f"Invalid ADNI subject code: {value}")
    return value


def session_code(value):
    value = str(value).strip().removeprefix("ses-")
    if value.lower() == "bl":
        return "M000"
    if re.fullmatch(r"[mM]\d+", value):
        return f"M{int(value[1:]):03d}"
    raise ValueError(f"Unsupported visit code: {value}")


def entities(path):
    # Inspect basename first; directory subject/session must agree if present.
    path = Path(path)
    values = dict(re.findall(r"(?:^|_)([A-Za-z]+)-([^_./]+)", path.name))
    for key in ("sub", "ses"):
        directory = [p[len(key) + 1:] for p in path.parts[:-1] if p.startswith(key + "-")]
        candidates = ([values[key]] if key in values else []) + directory
        normalise = subject_code if key == "sub" else session_code
        codes = {normalise(v) for v in candidates}
        if len(codes) != 1:
            raise ValueError(f"Missing or conflicting {key} identifiers: {path}")
        values[key] = codes.pop()
    return {key: values.get(key, "") for key in ("sub", "ses", *SCAN_ENTITIES)}


def mat_array(path, variable):
    """Decode MATLAB storage order, then leave scientific axes untouched."""
    if h5py.is_hdf5(path):
        with h5py.File(path, "r") as f:
            if variable not in f or not isinstance(f[variable], h5py.Dataset):
                raise ValueError(f"Missing numeric MATLAB dataset {variable}")
            ds = f[variable]
            if "MATLAB_class" not in ds.attrs:
                raise ValueError("HDF5 dataset lacks MATLAB_class; storage order is unverified")
            if ds.dtype.kind not in "fiu":
                raise ValueError("Expected a real numeric MATLAB dataset")
            a = np.asarray(ds)
            return a.transpose(tuple(reversed(range(a.ndim))))
    data = loadmat(path, variable_names=[variable])
    if variable not in data or not isinstance(data[variable], np.ndarray):
        raise ValueError(f"Missing numeric MATLAB variable {variable}")
    a = data[variable]
    if a.dtype.kind not in "fiu":
        raise ValueError("Expected a real numeric MATLAB array")
    return a


def read_tensor(path, config):
    a = mat_array(path, config.inputs.tensor_variable)
    if a.ndim != 3:
        raise ValueError(f"Expected a three-dimensional tensor, got {a.shape}")
    a = np.transpose(a, config.inputs.tensor_axes)
    p = config.cohort.parcels
    if a.shape != (p, p, 16):
        raise ValueError(f"Expected {(p, p, 16)}, got {a.shape}; specify axes explicitly")
    mask = ~np.eye(p, dtype=bool) if config.transform.exclude_diagonal else np.ones((p, p), bool)
    if not np.isfinite(a[mask]).all():
        raise ValueError("Non-finite tensor entries in included connections")
    return np.asarray(a, dtype=np.float64)


def series_metadata(path, parcels, variable=None):
    if h5py.is_hdf5(path):
        with h5py.File(path, "r") as f:
            entries = []
            for name, ds in f.items():
                if isinstance(ds, h5py.Dataset) and ds.dtype.kind in "fiu" and "MATLAB_class" in ds.attrs:
                    entries.append((name, tuple(reversed(ds.shape)), "double"))
    else:
        entries = whosmat(path)
    numeric = {"double", "single", "int8", "uint8", "int16", "uint16", "int32", "uint32", "int64", "uint64"}
    choices = []
    for name, shape, kind in entries:
        if variable is not None and name != variable:
            continue
        if kind in numeric and len(shape) == 2 and shape.count(parcels) == 1:
            axis = 1 if shape[0] == parcels else 0
            choices.append((name, int(shape[axis]), tuple(shape)))
    if len(choices) != 1:
        raise ValueError(f"Expected one unambiguous parcel-by-time matrix; candidates={entries}")
    return choices[0]


def read_atlas(path, parcels):
    atlas = pd.read_csv(path, dtype={"parcel_id": str})
    required = {"parcel_index", "parcel_id", "name", "network_id"}
    if not required <= set(atlas):
        raise ValueError(f"Atlas requires {sorted(required)} in tensor order")
    if len(atlas) != parcels or not np.array_equal(atlas.parcel_index, np.arange(parcels)):
        raise ValueError("Atlas must contain exactly parcels rows, in explicit 0-based tensor order")
    if atlas.parcel_id.isna().any() or atlas.parcel_id.duplicated().any() or atlas.name.isna().any():
        raise ValueError("Atlas parcel IDs must be unique and nonmissing; names are required")
    if not atlas.network_id.isin(range(1, 9)).all():
        raise ValueError("Atlas network_id must be an integer from 1 to 8")
    return atlas


def clinical_table(path):
    df = pd.read_csv(path, low_memory=False)
    if not {"PTID", "VISCODE"} <= set(df):
        raise ValueError("Clinical table requires PTID and VISCODE")
    if "Years.bl" in df and "Years_bl" not in df:
        df = df.rename(columns={"Years.bl": "Years_bl"})
    df["clinical_row"] = np.arange(len(df))
    issues = []
    subject, session = [], []
    for idx, row in df.iterrows():
        try:
            subject.append(subject_code(row.PTID))
            session.append(session_code(row.VISCODE))
        except ValueError as e:
            if len(subject) > len(session):
                subject.pop()
            subject.append(None)
            session.append(None)
            issues.append({"clinical_row": int(idx), "reason": str(e)})
    df["subject"], df["session"] = subject, session
    for col in ("DX", "PTGENDER", "PTEDUCAT", "AGE", "Years_bl", "ICV", "MMSE", "MOCA", "EXAMDATE"):
        if col not in df:
            df[col] = np.nan
    for col in ("PTEDUCAT", "AGE", "Years_bl", "ICV", "MMSE", "MOCA"):
        df[col] = pd.to_numeric(df[col], errors="coerce").replace([np.inf, -np.inf], np.nan)
    df.loc[df.ICV <= 0, "ICV"] = np.nan
    df["group"] = df.DX.replace({"Dementia": "AD"})
    df["sex"] = df.PTGENDER.astype("string").str.strip().str.lower().replace({"male": "Male", "female": "Female"})
    df.loc[~df.sex.isin(["Male", "Female"]), "sex"] = pd.NA
    df["education"] = df.PTEDUCAT
    df["visit_age"] = df.AGE + df.Years_bl
    df["EXAMDATE"] = pd.to_datetime(df.EXAMDATE, errors="coerce").dt.strftime("%Y-%m-%d")
    for score, bins in (("MMSE", [-1, 19, 27, 30]), ("MOCA", [-1, 19, 28, 30])):
        labels = ["0-19", "20-27", "28-30"] if score == "MMSE" else ["0-19", "20-28", "29-30"]
        valid = df[score].where(df[score].between(0, 30) & (df[score] % 1 == 0))
        df[f"{score}_category"] = pd.cut(valid, bins, labels=labels).astype("string")
    return df, issues


def etiv_lookup(root):
    values, files = {}, []
    if not root:
        return values, files
    pattern = re.compile(r"#\s*Measure\s+.*EstimatedTotalIntraCranialVol.*?,\s*([\d.eE+-]+)\s*,\s*mm\^3")
    for path in sorted(Path(root).glob("sub-*/stats/aseg.stats")):
        files.append(path)
        try:
            subject = subject_code(path.parts[-3])
            found = [float(m.group(1)) for line in path.read_text().splitlines() if (m := pattern.match(line))]
            if len(found) == 1 and np.isfinite(found[0]) and found[0] > 0:
                values.setdefault(subject, []).append(found[0])
        except (ValueError, OSError):
            continue
    # Conflicting values are deliberately not usable for recovery.
    return {s: v[0] for s, v in values.items() if len(set(v)) == 1}, files


@dataclass
class Inventory:
    scans: pd.DataFrame
    atlas: pd.DataFrame
    sources: dict
    issues: list

    def summary(self):
        df = self.scans
        return {"discovered": len(df), "eligible_visits": int(df.eligible.sum()),
                "selected_scans": int(df.selected.sum()),
                "exclusions": df.loc[~df.eligible, "reason"].value_counts().to_dict(),
                "selected_groups": df.loc[df.selected, "group"].fillna("missing").value_counts().to_dict(),
                "length_rule": "inclusive time-series-length filter; acquisition protocol is not verified",
                "issues": self.issues}


def source_files(config):
    inp, c = config.inputs, config.cohort
    for key in ("tensors", "series"):
        if not Path(getattr(inp, key)).is_dir() or not getattr(inp, key):
            raise FileNotFoundError(f"Missing {key} directory: {getattr(inp, key)}")
    for key in ("clinical", "atlas"):
        if not getattr(inp, key) or not Path(getattr(inp, key)).is_file():
            raise FileNotFoundError(f"Missing {key} file: {getattr(inp, key)}")
    for key in ("freesurfer", "resolutions"):
        value = getattr(inp, key)
        if value and not Path(value).exists():
            raise FileNotFoundError(f"Configured {key} does not exist: {value}")
    desc = "-CCS" if c.estimator == "ccs" else ""
    tensors = sorted(Path(inp.tensors).rglob(inp.tensor_glob or f"*{c.parcels}_desc{desc}-phiidtensor.mat"))
    series = sorted(Path(inp.series).rglob(inp.series_glob or f"*ExtSchaefer{c.parcels}_ROI_activity.mat"))
    if not tensors:
        raise ValueError("No tensor files match the configured atlas and estimator")
    stats = sorted(Path(inp.freesurfer).glob("sub-*/stats/aseg.stats")) if inp.freesurfer else []
    extra = [Path(inp.clinical), Path(inp.atlas)] + ([Path(inp.resolutions)] if inp.resolutions else [])
    return tensors, series, sorted(set(tensors + series + stats + extra))


def fingerprints(config):
    _, _, files = source_files(config)
    return {str(p.resolve()): {"sha256": file_hash(p), "bytes": p.stat().st_size} for p in files}


def inventory(config, sources=None):
    config = ensure_config(config)
    inp, c = config.inputs, config.cohort
    tensor_files, series_files, _ = source_files(config)
    atlas = read_atlas(inp.atlas, c.parcels)
    clinical, issues = clinical_table(inp.clinical)
    clinical_groups = {key: sub for key, sub in clinical.dropna(subset=["subject", "session"]).groupby(["subject", "session"])}
    etiv, _ = etiv_lookup(inp.freesurfer)
    lookup = {}
    for path in series_files:
        try:
            ent = entities(path)
            lookup.setdefault((ent["sub"], ent["ses"]), []).append((path, ent))
        except ValueError as e:
            issues.append({"file": str(path), "reason": str(e)})
    resolutions = {}
    if inp.resolutions:
        resolution_file = Path(inp.resolutions)
        df = pd.read_csv(resolution_file, keep_default_na=False)
        if "tensor_file" not in df:
            raise ValueError("Resolution table requires tensor_file")
        for _, row in df.iterrows():
            key = str((resolution_file.parent / str(row.tensor_file)).resolve())
            if key in resolutions:
                raise ValueError(f"Duplicate resolution for {key}")
            record = row.to_dict()
            if record.get("series_file"):
                record["series_file"] = str((resolution_file.parent / record["series_file"]).resolve())
            resolutions[key] = record
        unmatched = set(resolutions) - {str(p.resolve()) for p in tensor_files}
        if unmatched:
            raise ValueError(f"Resolution references undiscovered tensors: {sorted(unmatched)}")
    rows = []
    for path in tensor_files:
        row = {"scan_id": digest(str(path.resolve()))[:20], "tensor_file": str(path.resolve()),
               "series_file": "", "n_time": np.nan, "eligible": False, "selected": False,
               "reason": "", "subject": "", "session": "", "month": np.nan, "group": np.nan}
        resolution = resolutions.get(str(path.resolve()), {})
        try:
            ent = entities(path)
            row.update({"subject": ent["sub"], "session": ent["ses"], **{k: ent[k] for k in SCAN_ENTITIES}})
            row["month"] = int(ent["ses"][1:])
            if str(resolution.get("use_scan", "true")).lower() in {"false", "0", "no"}:
                raise ValueError("excluded_by_resolution")
            candidates = [(p, e) for p, e in lookup.get((ent["sub"], ent["ses"]), [])
                          if all(not ent[k] or not e[k] or ent[k] == e[k] for k in SCAN_ENTITIES)]
            if resolution.get("series_file"):
                candidates = [(p, e) for p, e in candidates if str(p.resolve()) == resolution["series_file"]]
            if len(candidates) != 1:
                raise ValueError("missing_series" if not candidates else "ambiguous_series")
            series_path, se = candidates[0]
            row["series_file"] = str(series_path.resolve())
            row["series_identity_partial"] = any(bool(ent[k]) != bool(se[k]) for k in SCAN_ENTITIES)
            try:
                var, count, shape = series_metadata(series_path, c.parcels, inp.series_variable)
            except (ValueError, OSError, TypeError) as e:
                raise ValueError(f"unreadable_or_ambiguous_series: {e}") from e
            row.update(n_time=count, series_variable=var, series_shape="x".join(map(str, shape)))
            candidates_clin = clinical_groups.get((ent["sub"], ent["ses"]), pd.DataFrame())
            if candidates_clin.empty:
                raise ValueError("missing_clinical_record")
            if resolution.get("clinical_row", "") != "":
                candidates_clin = candidates_clin.loc[candidates_clin.clinical_row == int(resolution["clinical_row"])]
            if candidates_clin.empty:
                raise ValueError("missing_clinical_record")
            unique = candidates_clin.drop(columns="clinical_row").drop_duplicates()
            if len(unique) != 1:
                raise ValueError("conflicting_clinical_records")
            cr = candidates_clin.iloc[0]
            row.update({k: cr[k] for k in ("DX", "group", "sex", "education", "AGE", "Years_bl", "visit_age", "ICV", "MMSE", "MOCA", "MMSE_category", "MOCA_category", "EXAMDATE", "clinical_row")})
            row["clinical_duplicate_count"] = len(candidates_clin)
            row["ICV_source"] = "ADNIMERGE" if pd.notna(cr.ICV) else "missing"
            if pd.isna(cr.ICV) and ent["sub"] in etiv:
                row["ICV"], row["ICV_source"] = etiv[ent["sub"]], "FreeSurfer_eTIV"
            if not c.time_min <= count <= c.time_max:
                raise ValueError("below_length_window" if count < c.time_min else "above_length_window")
            try:
                read_tensor(path, config)
            except (ValueError, OSError, KeyError, TypeError) as e:
                raise ValueError(f"invalid_tensor: {e}") from e
            row["eligible"] = True
            row["reason"] = "eligible"
        except ValueError as e:
            row["reason"] = str(e)
        rows.append(row)
    scans = pd.DataFrame(rows)
    # A source series cannot identify multiple distinct tensors under one estimator.
    eligible = scans.loc[scans.eligible]
    duplicates = eligible.duplicated(["subject", "session"], keep=False)
    idx = eligible.index[duplicates]
    scans.loc[idx, ["eligible", "reason"]] = [False, "ambiguous_same_visit_scans"]
    eligible = scans.loc[scans.eligible]
    if c.visit == "latest":
        selected = eligible.groupby("subject")["month"].idxmax()
    else:
        visit = "M000" if c.visit == "baseline" else session_code(c.visit)
        selected = eligible.index[eligible.session == visit]
    scans.loc[selected, "selected"] = True
    return Inventory(scans, atlas, sources if sources is not None else fingerprints(config), issues)
