from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
import numpy as np
import pandas as pd

from .config import ensure_config
from .inputs import fingerprints, inventory, read_tensor
from .numerics import (Design, complete_cases, normalise_tdmi, fit_coefficients,
                       apply_adjustment, regional_sums, network_means, global_tdmi)
from .provenance import atomic_json, code_hash, digest, environment, seal, transaction, utc_now, verify


@dataclass
class Prepared:
    """Shared acquisition store with a separate adjusted model for each grouping."""
    path: Path
    config: object
    grouping: str = "diagnosis"

    def __post_init__(self):
        if self.grouping not in {"diagnosis", "moca"}:
            raise ValueError("Grouping must be diagnosis or moca")

    @property
    def groups(self):
        return ("CN", "MCI", "AD") if self.grouping == "diagnosis" else ("29-30", "20-28", "0-19")

    @property
    def pairs(self):
        from itertools import combinations
        return tuple(combinations(self.groups, 2))

    @property
    def group_column(self):
        return "group" if self.grouping == "diagnosis" else "MOCA_category"

    @property
    def model_path(self):
        return self.path if self.grouping == "diagnosis" else self.path / "moca"

    def for_group(self, grouping):
        return Prepared(self.path, self.config, grouping)

    def target(self, name):
        return target_id(name, self.grouping)

    def finish(self, result):
        result.metadata.update(grouping=self.grouping, grouping_column=self.group_column,
                               groups=self.groups, reference_group=self.groups[0])
        # Observation tables retain actual diagnosis alongside the analysis label.
        identities = self.table("inventory").set_index("scan_id")
        for table in result.tables.values():
            if "scan_id" in table:
                for column, source in (("diagnosis", "group"), ("MOCA_category", "MOCA_category")):
                    table[column] = table.scan_id.map(identities[source])
                table["grouping"] = self.grouping
        for recipe in result.figures:
            prefix, rest = recipe["id"].split(".", 1)
            recipe["id"] = self.target(prefix) + "." + rest
        return result

    @property
    def manifest(self):
        return json.loads((self.path / "manifest.json").read_text())

    def table(self, name):
        root = self.model_path if name in {"subjects", "selected"} else self.path
        table = pd.read_csv(root / f"{name}.csv", keep_default_na=True)
        if self.group_column in table:
            table["analysis_group"] = table[self.group_column]
        return table

    def array(self, name):
        root = self.path if name.startswith("raw") else self.model_path
        return np.load(root / f"{name}.npy", mmap_mode="r", allow_pickle=False)


def target_id(name, grouping):
    if grouping == "diagnosis":
        return name
    if name == "cn-reference":
        return "moca-reference"
    if name.startswith("taxa-"):
        return name.removesuffix("-dx") + "-moca"
    return name + "-moca"


def prepare(config):
    config = ensure_config(config)
    sources = fingerprints(config)
    settings = config.to_dict()
    runtime = environment()
    key = digest({"inputs": sources, "configuration": {k: settings[k] for k in ("inputs", "cohort", "transform")},
                  "implementation": code_hash(), "python": runtime["python"], "packages": runtime["packages"]})
    destination = Path(config.cache) / key
    if destination.exists():
        if not verify(destination):
            raise ValueError(f"Prepared cache is incomplete or corrupt: {destination}; move it aside and retry")
        return Prepared(destination, config)
    inv = inventory(config, sources=sources)
    visits = inv.scans.loc[inv.scans.eligible].copy().reset_index(drop=True)
    if visits.empty:
        raise ValueError(f"No eligible scans. Audit: {inv.summary()}")
    visits["tensor_row"] = np.arange(len(visits))
    selected = visits.loc[visits.selected].copy().reset_index(drop=True)
    if selected.empty:
        raise ValueError("No scans satisfy the requested visit selection")
    p, t = config.cohort.parcels, config.transform
    shape = (len(visits), p, p, 16)
    with transaction(destination) as out:
        inv.scans.to_csv(out / "inventory.csv", index=False)
        inv.atlas.to_csv(out / "atlas.csv", index=False)
        visits.to_csv(out / "visits.csv", index=False)
        raw = np.lib.format.open_memmap(out / "raw.npy", mode="w+", dtype="float64", shape=shape)
        raw_region = np.lib.format.open_memmap(out / "raw_regional.npy", mode="w+", dtype="float64", shape=(len(visits), p, 16))
        nets = sorted(inv.atlas.network_id.unique())
        raw_network = np.lib.format.open_memmap(out / "raw_network.npy", mode="w+", dtype="float64", shape=(len(visits), len(nets), 16))
        tdmi_values = []
        for i, row in visits.iterrows():
            a = read_tensor(row.tensor_file, config)
            raw[i] = a
            raw_region[i] = regional_sums(a, t.exclude_diagonal, t.regional_axis)
            raw_network[i] = network_means(a, inv.atlas.network_id, t.exclude_diagonal, t.network_edges)
            tdmi_values.append(global_tdmi(a, t.exclude_diagonal))
        raw.flush(); raw_region.flush(); raw_network.flush()
        np.save(out / "raw_tdmi.npy", np.asarray(tdmi_values))
        audits = {}
        for grouping in ("diagnosis", "moca"):
            model_out = out if grouping == "diagnosis" else out / "moca"
            audits[grouping] = _prepare_model(model_out, selected, raw, inv.atlas, config, grouping)
        atomic_json(out / "audit.json", {**inv.summary(), **audits["diagnosis"], "groupings": audits})
        seal(out, {"kind": "preparation", "cache_id": key, "created_at": utc_now(),
                   "configuration": settings, "sources": sources, "environment": runtime,
                   "dimensions": {"parcels": p, "atoms": 16, "visits": len(visits), "subjects": audits["diagnosis"]["adjusted_subjects"], "moca_subjects": audits["moca"]["adjusted_subjects"]},
                   "network_ids": list(map(int, nets))})
    return Prepared(destination, config)


def _prepare_model(out, selected, raw, atlas, config, grouping):
    out.mkdir(parents=True, exist_ok=True)
    view = Prepared(out, config, grouping)
    selected = selected.copy()
    valid, reasons = complete_cases(selected, config.transform)
    valid &= selected[view.group_column].isin(view.groups).to_numpy()
    selected["adjustment_exclusion"] = [r if pd.notna(g) and g in view.groups else (r + f";missing_or_unsupported_{view.grouping}").strip(";")
                                        for r, g in zip(reasons, selected[view.group_column])]
    selected["adjusted_eligible"] = valid
    subjects = selected.loc[valid].copy().reset_index(drop=True)
    subjects["adjusted_row"] = np.arange(len(subjects))
    selected.to_csv(out / "selected.csv", index=False)
    subjects.to_csv(out / "subjects.csv", index=False)
    cohort_counts = {"adjusted_subjects": len(subjects),
              "adjustment_exclusions": selected.loc[~valid, "adjustment_exclusion"].value_counts().to_dict(),
              "selected_groups": selected[view.group_column].fillna("missing").value_counts().to_dict()}
    atomic_json(out / "cohort_audit.json", cohort_counts)
    p, t = config.cohort.parcels, config.transform
    nets = sorted(atlas.network_id.unique())
    normalisation_audit = []
    # Empty adjusted cohorts are valid preparation: raw/longitudinal targets can still run.
    adjusted = np.lib.format.open_memmap(out / "adjusted.npy", mode="w+", dtype="float64", shape=(len(subjects), p, p, 16))
    if len(subjects):
        design = Design.fit(subjects, t)
        X = design.matrix(subjects)
        for i, row in subjects.iterrows():
            a = np.asarray(raw[int(row.tensor_row)])
            count = 0
            if t.tdmi == "tdmi":
                a, count = normalise_tdmi(a, t.epsilon, t.exclude_diagonal)
            adjusted[i] = a
            normalisation_audit.append({"scan_id": row.scan_id, "invalid_tdmi_edges": count})
        Y = adjusted.reshape(len(subjects), -1)
        beta = np.lib.format.open_memmap(out / "coefficients.npy", mode="w+", dtype="float64", shape=(X.shape[1], p*p*16))
        nobs = np.empty(p*p*16, dtype=int)
        for start in range(0, Y.shape[1], t.block_features):
            end = min(Y.shape[1], start + t.block_features)
            b, counts = fit_coefficients(X, Y[:, start:end])
            beta[:, start:end], nobs[start:end] = b, counts
            Y[:, start:end] = apply_adjustment(Y[:, start:end], X, b, design.reference, t.adjustment)
        beta.flush()
        np.save(out / "coefficient_nobs.npy", nobs.reshape(p, p, 16))
        np.save(out / "design_matrix.npy", X)
        imputed = []
        if t.missing == "median":
            for i, row in subjects.iterrows():
                for name in t.confounds:
                    source = "AGE" if name == "visit_age" and t.age == "baseline" else name
                    if pd.isna(row[source]):
                        imputed.append({"scan_id": row.scan_id, "confound": name,
                                        "replacement": "missing_category" if name == "sex" else design.medians[name]})
        atomic_json(out / "design.json", {**design.to_dict(), "rank": int(np.linalg.matrix_rank(X)),
                    "parameter_count": X.shape[1], "rank_deficient": np.linalg.matrix_rank(X) < X.shape[1],
                    "imputations": imputed, "adjustment": t.adjustment})
    else:
        np.save(out / "coefficients.npy", np.empty((0, p*p*16)))
        np.save(out / "coefficient_nobs.npy", np.zeros((p, p, 16), int))
        np.save(out / "design_matrix.npy", np.empty((0, 0)))
        atomic_json(out / "design.json", {"status": "insufficient_data", "columns": []})
    adjusted.flush()
    region = np.lib.format.open_memmap(out / "regional.npy", mode="w+", dtype="float64", shape=(len(subjects), p, 16))
    network = np.lib.format.open_memmap(out / "network.npy", mode="w+", dtype="float64", shape=(len(subjects), len(nets), 16))
    for i in range(len(subjects)):
        region[i] = regional_sums(adjusted[i], t.exclude_diagonal, t.regional_axis)
        network[i] = network_means(adjusted[i], atlas.network_id, t.exclude_diagonal, t.network_edges)
    region.flush(); network.flush()
    atomic_json(out / "normalisation_audit.json", normalisation_audit)
    return cohort_counts
