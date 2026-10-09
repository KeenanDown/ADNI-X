from __future__ import annotations

from itertools import combinations
import json
import warnings
import numpy as np
import pandas as pd
from scipy import stats

from .constants import (ATOMS, GROUPS, PAIRS, NETWORKS, TAXA_FOUR, TAXA_SIX,
                        TAXA_PRINCIPAL, MOCA_GROUPS, MMSE_GROUPS)
from .inference import validate_symmetric_atoms, hc3_coefficients
from .numerics import (compare, fdr, omnibus, hierarchy, correlation, cn_standardise,
                       interval, permutation_p, connection_mask)
from .provenance import rng_for
from .results import Result, InsufficientData


def context(prepared, minimum=1):
    df, atlas = prepared.table("subjects"), prepared.table("atlas")
    if len(df) < minimum:
        raise InsufficientData(f"Need at least {minimum} adjusted participants; found {len(df)}")
    return df, atlas, prepared.config.statistics


def analysis_view(prepared, grouping):
    """Select the independently fitted model before reading adjusted arrays."""
    target = {"group": "diagnosis", "MOCA_category": "moca"}.get(grouping)
    return prepared.for_group(target) if target and prepared.grouping != target else prepared


def grouped_context(prepared, grouping):
    df, atlas, settings = context(prepared, 2)
    grouping = grouping or prepared.group_column
    if grouping == "MMSE_category":
        # Categorise MMSE only when explicitly requested; default preparation
        # retains the source score without calculating an MMSE category.
        score = df.MMSE.where(df.MMSE.between(0, 30) & (df.MMSE % 1 == 0))
        df = df.copy()
        df[grouping] = pd.cut(score, [-1, 19, 27, 30],
                             labels=["0-19", "20-27", "28-30"]).astype(object)
    groups = {"group": GROUPS, "MOCA_category": MOCA_GROUPS,
              "MMSE_category": MMSE_GROUPS}[grouping]
    df = df.loc[df[grouping].isin(groups)].copy()
    if len(df) < 2:
        raise InsufficientData(f"Need at least two participants with {grouping}")
    return df, atlas, settings, groups


def correct(table, family_columns, settings):
    table = table.copy()
    table["q"], table["significant"] = np.nan, False
    grouping = family_columns[0] if len(family_columns) == 1 else family_columns
    groups = table.groupby(grouping, dropna=False).groups if family_columns else {"all": table.index}
    for key, idx in groups.items():
        q, reject = fdr(table.loc[idx, "p"], settings.alpha, settings.correction)
        table.loc[idx, "q"], table.loc[idx, "significant"] = q, reject
    return table


def _compare_columns(a, b, welch):
    # Vectorised equivalent of compare(), used for millions of edge tests.
    na, nb = np.isfinite(a).sum(axis=0), np.isfinite(b).sum(axis=0)
    with warnings.catch_warnings(), np.errstate(all="ignore"):
        warnings.simplefilter("ignore", RuntimeWarning)
        ma, mb = np.nanmean(a, axis=0), np.nanmean(b, axis=0)
        va, vb = np.nanvar(a, axis=0, ddof=1), np.nanvar(b, axis=0, ddof=1)
        statistic, p = stats.ttest_ind(b, a, axis=0, equal_var=not welch, nan_policy="omit")
        pooled = np.sqrt(((na - 1)*va + (nb - 1)*vb) / (na + nb - 2))
        d = (mb - ma) / pooled
        g = d * (1 - 3 / (4*(na + nb - 2) - 1))
    valid = (na >= 2) & (nb >= 2) & ((va + vb) > 0)
    return {"n_a": na, "n_b": nb, "mean_a": ma, "mean_b": mb, "difference": mb - ma,
            "t": np.where(valid, statistic, np.nan), "p": np.where(valid, p, np.nan),
            "cohens_d": np.where(valid, d, np.nan), "hedges_g": np.where(valid, g, np.nan),
            "status": np.where((na < 2) | (nb < 2), "insufficient_samples", np.where(valid, "ok", "zero_variance"))}


def group_matrices(prepared, grouping=None):
    prepared = analysis_view(prepared, grouping)
    grouping = grouping or prepared.group_column
    df, atlas, settings, groups = grouped_context(prepared, grouping)
    data = prepared.array("adjusted")
    p = len(atlas)
    symmetry = validate_symmetric_atoms(data, df.index.to_numpy(),
        prepared.config.transform.exclude_diagonal, prepared.config.transform.block_features)
    means = np.full((len(groups), p, p, 16), np.nan)
    for gi, group in enumerate(groups):
        indices = df.index[df[grouping] == group].to_numpy()
        if len(indices):
            total, count = np.zeros((p, p, 16)), np.zeros((p, p, 16), int)
            for i in indices:
                valid = np.isfinite(data[i])
                total += np.where(valid, data[i], 0)
                count += valid
            np.divide(total, count, out=means[gi], where=count > 0)
    mask = connection_mask(p, prepared.config.transform.exclude_diagonal, upper=True)
    ii, jj = np.where(mask)
    flat = data.reshape(len(data), -1)
    rows = []
    for atom in (0, 15):
        columns = (ii*p + jj)*16 + atom
        for a, b in combinations(groups, 2):
            ia = df.index[df[grouping] == a].to_numpy()
            ib = df.index[df[grouping] == b].to_numpy()
            for start in range(0, len(columns), prepared.config.transform.block_features):
                end = start + prepared.config.transform.block_features
                cols = columns[start:end]
                comparisons = _compare_columns(flat[np.ix_(ia, cols)], flat[np.ix_(ib, cols)], settings.welch)
                rows.append(pd.DataFrame({"atom": ATOMS[atom], "group_a": a, "group_b": b,
                    "parcel_i": ii[start:end], "parcel_j": jj[start:end], **comparisons}))
    table = correct(pd.concat(rows, ignore_index=True), ["atom", "group_a", "group_b"], settings)
    result = Result(tables={"edges": table, "atlas": atlas,
                           "group_counts": df.groupby(grouping).size().rename("n").reset_index().rename(columns={grouping: "group"})},
                    arrays={"group_means": means}, metadata={"groups": groups, "grouping": grouping,
                        "cohort": "adjusted", "family": "unique ROI pairs per atom/contrast",
                        "edge_convention": "unique_upper", "symmetry": symmetry})
    name = "group-matrices"
    for atom in (0, 15):
        result.figure(f"{name}.{ATOMS[atom]}", "group_matrices", atom=atom,
                      group_a=groups[0], group_b=groups[-1])
    return prepared.finish(result)


def coefficients(prepared):
    df, atlas, settings = context(prepared, 2)
    model = json.loads((prepared.model_path / "design.json").read_text())
    beta = np.asarray(prepared.array("coefficients"))
    names = ["Intercept", *model["columns"]]
    X = np.asarray(prepared.array("design_matrix"))
    p, transform = len(atlas), prepared.config.transform
    symmetry = validate_symmetric_atoms(prepared.array("adjusted"), df.index.to_numpy(),
        transform.exclude_diagonal, transform.block_features)
    ii, jj = np.where(connection_mask(p, transform.exclude_diagonal, upper=True))
    raw = prepared.array("raw").reshape(-1, p*p*16)
    tensor_rows = df.tensor_row.to_numpy(int)
    tables = []
    for atom in (0, 15):
        block_size = max(1, transform.block_features//16) if transform.tdmi == "tdmi" else transform.block_features
        for start in range(0, len(ii), block_size):
            i, j = ii[start:start+block_size], jj[start:start+block_size]
            edges = i*p+j
            columns = edges*16+atom
            Y = np.array(raw[np.ix_(tensor_rows, columns)])
            if transform.tdmi == "tdmi":
                # The same raw 16-atom denominator used during preparation.
                atom_columns = (edges[:, None]*16+np.arange(16)).ravel()
                original = raw[np.ix_(tensor_rows, atom_columns)].reshape(len(df), len(edges), 16)
                denominator = original.sum(axis=-1)
                valid = np.isfinite(denominator) & (abs(denominator) >= transform.epsilon)
                normalised = np.full(Y.shape, np.nan)
                np.divide(Y, denominator, out=normalised, where=valid)
                Y = normalised
            fitted = beta[:, columns]
            inference = hc3_coefficients(X, Y, fitted)
            for k, name in enumerate(names[1:], 1):
                tables.append(pd.DataFrame({"atom": ATOMS[atom], "coefficient": name,
                    "parcel_i": i, "parcel_j": j, "estimate": fitted[k],
                    "nobs": inference["nobs"], "df_resid": inference["df_resid"],
                    **{key: inference[key][k] for key in
                       ("standard_error", "t", "p", "ci_lower", "ci_upper", "status")}}))
    columns = ["atom", "coefficient", "parcel_i", "parcel_j", "estimate", "nobs", "df_resid",
               "standard_error", "t", "p", "ci_lower", "ci_upper", "status"]
    inference_table = correct(pd.concat(tables, ignore_index=True) if tables else
                              pd.DataFrame(columns=columns), ["atom"], settings)
    result = Result(tables={"atlas": atlas, "design": pd.DataFrame(prepared.array("design_matrix"), columns=names)},
                    arrays={"coefficients": beta.reshape(len(beta), len(atlas), len(atlas), 16),
                            "nobs": np.asarray(prepared.array("coefficient_nobs"))},
                    metadata={"design": model, "names": names, "cohort": "adjusted",
                        "symmetry": symmetry, "coefficient_inference": "HC3; two-sided approximate Student-t; 95% CI",
                        "family": "all non-intercept coefficients x unique ROI pairs, per atom",
                        "tested_atoms": ["rtr", "sts"], "intercept_tested": False})
    result.tables["inference"] = inference_table
    for atom in (0, 15):
        result.figure(f"coefficients.{ATOMS[atom]}", "coefficients", atom=atom)
        result.figure(f"coefficients.{ATOMS[atom]}.significant", "coefficients", atom=atom, significant=True)
    return prepared.finish(result)


def network_effects(prepared, grouping=None):
    prepared = analysis_view(prepared, grouping)
    grouping = grouping or prepared.group_column
    df, atlas, settings, groups = grouped_context(prepared, grouping)
    values = prepared.array("network")
    nets = sorted(atlas.network_id.unique())
    rows, observations = [], []
    for ni, net in enumerate(nets):
        for atom, atom_name in enumerate(ATOMS):
            for i, row in df.iterrows():
                observations.append({"scan_id": row.scan_id, "subject": row.subject, "group": row.analysis_group,
                                     "network_id": net, "network": NETWORKS[net], "atom": atom_name, "value": values[i, ni, atom]})
                observations[-1]["group"] = row[grouping]
            for a, b in combinations(groups, 2):
                ia, ib = df.index[df[grouping] == a], df.index[df[grouping] == b]
                rows.append({"network_id": net, "network": NETWORKS[net], "atom": atom_name, "group_a": a, "group_b": b,
                             **compare(values[ia, ni, atom], values[ib, ni, atom], settings.welch)})
    table = correct(pd.DataFrame(rows), ["group_a", "group_b"], settings)
    result = Result(tables={"effects": table, "observations": pd.DataFrame(observations)},
                    metadata={"cohort": "adjusted", "groups": groups, "grouping": grouping,
                              "family": "all 16 atoms x all networks per contrast"})
    name = "network-effects"
    for a, b in combinations(groups, 2):
        for effect in ("cohens_d", "hedges_g"):
            for subset in ("all-atoms", "rtr-sts"):
                result.figure(f"{name}.{b}-vs-{a}.{effect}.{subset}", "effects_grid",
                              table="effects", group_a=a, group_b=b, effect=effect,
                              atoms=list(ATOMS) if subset == "all-atoms" else ["rtr", "sts"])
    return prepared.finish(result)


def regional_effects(prepared):
    df, atlas, settings = context(prepared, 2)
    values = prepared.array("regional")
    rows = []
    for atom in (0, 15):
        for a, b in prepared.pairs:
            for region in range(len(atlas)):
                rows.append({"parcel_index": region, "parcel_id": atlas.parcel_id.iloc[region],
                             "network_id": atlas.network_id.iloc[region], "atom": ATOMS[atom], "group_a": a, "group_b": b,
                             **compare(values[df.analysis_group == a, region, atom], values[df.analysis_group == b, region, atom], settings.welch)})
    table = correct(pd.DataFrame(rows), ["atom", "group_a", "group_b"], settings)
    counts = []
    for (atom, a, b, net), sub in table.groupby(["atom", "group_a", "group_b", "network_id"]):
        for direction, sign in (("increase", 1), ("decrease", -1)):
            n = int((sub.significant & (sub.difference*sign > 0)).sum())
            counts.append({"atom": atom, "group_a": a, "group_b": b, "network_id": net,
                           "direction": direction, "n_significant": n, "n_regions": len(sub), "proportion": n / len(sub)})
    net_rows = []
    for net in sorted(atlas.network_id.unique()):
        idx = np.flatnonzero(atlas.network_id == net)
        for atom in (0, 15):
            net_values = values[:, idx, atom].mean(axis=1)
            for a, b in prepared.pairs:
                net_rows.append({"network_id": net, "atom": ATOMS[atom], "group_a": a, "group_b": b,
                                 **compare(net_values[df.analysis_group == a], net_values[df.analysis_group == b], settings.welch)})
    net_table = correct(pd.DataFrame(net_rows), ["atom", "group_a", "group_b"], settings)
    result = Result(tables={"effects": table, "network_counts": pd.DataFrame(counts), "network_regional_effects": net_table, "atlas": atlas},
                    metadata={"cohort": "adjusted", "family": "regions per atom/contrast; network regional averages per atom/contrast"})
    for atom in ("rtr", "sts"):
        for a, b in prepared.pairs:
            result.figure(f"regional-effects.{atom}.{b}-vs-{a}", "significance_strip", atom=atom, group_a=a, group_b=b)
        result.figure(f"regional-effects.{atom}.network-counts", "network_counts", atom=atom)
    return prepared.finish(result)


def hierarchies(prepared, grouping=None):
    prepared = analysis_view(prepared, grouping)
    grouping = grouping or prepared.group_column
    df, atlas, settings, names = grouped_context(prepared, grouping)
    values = np.asarray(prepared.array("regional"))[df.index.to_numpy()]
    labels = df[grouping]
    means, ranks, rows, correlations = {}, {}, [], []
    for group in names:
        sub = values[labels == group]
        means[group] = sub.mean(axis=0) if len(sub) else np.full(values.shape[1:], np.nan)
        ranks[group] = hierarchy(means[group])
        for r in range(len(atlas)):
            rows.append({"group": group, "parcel_index": r, "parcel_id": atlas.parcel_id.iloc[r],
                         "network_id": atlas.network_id.iloc[r], "hierarchy": ranks[group][r],
                         "rtr": means[group][r, 0], "sts": means[group][r, 15], "n": len(sub)})
    name = "hierarchies"
    rng = rng_for(settings.seed, prepared.target(name))
    pairs = list(combinations(names, 2))
    samples = np.full((settings.resamples, len(pairs)), np.nan)
    grouped = {g: values[labels == g] for g in names}
    for rep in range(settings.resamples):
        sampled_ranks = {}
        for g, sub in grouped.items():
            sampled_ranks[g] = hierarchy(sub[rng.integers(len(sub), size=len(sub))].mean(axis=0)) if len(sub) >= 2 else np.full(len(atlas), np.nan)
        for pi, (a, b) in enumerate(pairs):
            samples[rep, pi] = correlation(sampled_ranks[a], sampled_ranks[b])
    for pi, (a, b) in enumerate(pairs):
        lower, upper = interval(samples[:, pi])
        value = correlation(ranks[a], ranks[b])
        correlations.append({"group_a": a, "group_b": b, "pearson_r": value, "ci_lower": lower, "ci_upper": upper,
                             "valid_resamples": np.isfinite(samples[:, pi]).sum(),
                             "status": "ok" if np.isfinite(value) and np.isfinite(lower) else "insufficient_or_constant_data"})
    profile = pd.DataFrame(rows)
    for atom in ("rtr", "sts"):
        baseline = profile.loc[profile.group == names[0]].set_index("parcel_index")[atom]
        profile[f"{atom}_deviation"] = profile[atom] - profile.parcel_index.map(baseline)
        endpoints = profile.loc[profile.group.isin([names[0], names[-1]])].groupby("parcel_index")[atom].mean()
        profile[f"{atom}_midpoint"] = profile[atom] - profile.parcel_index.map(endpoints)
    result = Result(tables={"profiles": profile, "correlations": pd.DataFrame(correlations), "atlas": atlas},
                    arrays={"bootstrap_correlations": samples}, metadata={"groups": names, "pairs": pairs,
                        "cohort": "adjusted", "grouping": grouping, "profile_baseline": names[0],
                        "deviation_definition": "group mean minus per-parcel baseline group mean",
                        "profile_interpretation": "descriptive; ordering adds no inferential evidence"})
    for measure in ("hierarchy", "rtr", "sts", "rtr_deviation", "sts_deviation", "rtr_midpoint", "sts_midpoint"):
        result.figure(f"{name}.{measure}", "profiles", measure=measure, anchor=names[0],
                      order="difference" if "_" in measure else "reference")
    result.figure(f"{name}.correlations", "correlation_intervals")
    return prepared.finish(result)


def taxa(prepared, six=False, grouping=None, principal=False):
    prepared = analysis_view(prepared, grouping)
    grouping = grouping or prepared.group_column
    df, atlas, settings, groups = grouped_context(prepared, grouping)
    values = np.asarray(prepared.array("network"))
    taxonomy = TAXA_PRINCIPAL if principal else TAXA_SIX if six else TAXA_FOUR
    nets = sorted(atlas.network_id.unique())
    observations, tests, posthoc = [], [], []
    for taxon, atoms in taxonomy.items():
        vals = values[:, :, list(atoms)].sum(axis=2)
        for ni, net in enumerate(nets):
            samples = [vals[df.index[df[grouping] == g], ni] for g in groups]
            for i, subject in df.iterrows():
                observations.append({"scan_id": subject.scan_id, "subject": subject.subject,
                                     "group": subject[grouping], "taxon": taxon, "network_id": net,
                                     "network": NETWORKS[net], "value": vals[i, ni]})
            tests.append({"taxon": taxon, "network_id": net, **omnibus(samples, settings.welch)})
            for ai, bi in combinations(range(len(groups)), 2):
                posthoc.append({"taxon": taxon, "network_id": net, "group_a": groups[ai], "group_b": groups[bi],
                                **compare(samples[ai], samples[bi], settings.welch)})
    tests = correct(pd.DataFrame(tests), [], settings)
    posthoc = correct(pd.DataFrame(posthoc), [], settings)
    posthoc = posthoc.merge(tests[["taxon", "network_id", "q"]].rename(columns={"q": "omnibus_q"}), on=["taxon", "network_id"], validate="many_to_one")
    suffix = "mmse" if grouping == "MMSE_category" else "dx"
    name = f"taxa-{'principal' if principal else 'six' if six else 'four'}-{suffix}"
    result = Result(tables={"observations": pd.DataFrame(observations), "omnibus": tests, "posthoc": posthoc},
                    metadata={"groups": groups, "grouping": grouping, "cohort": "adjusted", "taxonomy": taxonomy,
                              "display_order": list(taxonomy),
                              "family": "all taxon/network panels; all panel/contrast posthoc tests separately"})
    result.figure(f"{name}.panels", "taxa_panels")
    result = prepared.finish(result)
    if grouping == "MMSE_category":
        result.metadata.update(grouping="mmse", grouping_column=grouping,
                               groups=groups, reference_group=groups[0])
        result.tables["observations"]["grouping"] = "mmse"
    return result


def export_label(group):
    """MATLAB identifiers must begin with a letter and contain no hyphens."""
    return group if group in {"CN", "MCI", "AD"} else "MOCA_" + group.replace("-", "_")


def cn_reference(prepared):
    df, atlas, _ = context(prepared, 2)
    values = np.asarray(prepared.array("regional"))
    cn = values[df.analysis_group == prepared.groups[0]]
    if len(cn) < 2:
        raise InsufficientData(f"Reference standardisation needs at least two {prepared.groups[0]} participants")
    z, mean, sd, valid = cn_standardise(values, cn)
    rows, subject_rows, pcs, meta = [], [], [], []
    groups = {}
    for group in prepared.groups:
        group_values = z[df.analysis_group == group]
        groups[group] = group_values.mean(axis=0) if len(group_values) else np.full(values.shape[1:], np.nan)
        for r in range(len(atlas)):
            rows.append({"group": group, "parcel_index": r, "parcel_id": atlas.parcel_id.iloc[r],
                         "network_id": atlas.network_id.iloc[r], "rtr": groups[group][r, 0], "sts": groups[group][r, 15]})
    for i, row in df.iterrows():
        subject_rows.append({"scan_id": row.scan_id, "subject": row.subject, "group": row.analysis_group,
                             "rtr": np.mean(z[i, :, 0]), "sts": np.mean(z[i, :, 15])})
    for group in prepared.groups[1:]:
        xy = groups[group][:, [0, 15]]
        keep = np.isfinite(xy).all(axis=1)
        if keep.sum() < 3:
            meta.append({"group": group, "status": "insufficient_reference_variance"})
            continue
        centre = xy[keep].mean(axis=0)
        eigval, eigvec = np.linalg.eigh(np.cov(xy[keep].T))
        if eigval.sum() <= 0 or eigval[-1] - eigval[0] <= np.finfo(float).eps * max(1, abs(eigval[-1])):
            meta.append({"group": group, "status": "degenerate_PCA"})
            continue
        pc1 = eigvec[:, -1]
        if pc1[1] < 0 or (pc1[1] == 0 and pc1[0] < 0):
            pc1 = -pc1
        scores = np.full(len(atlas), np.nan)
        scores[keep] = (xy[keep] - centre) @ pc1
        meta.append({"group": group, "status": "ok", "centre": centre.tolist(), "pc1": pc1.tolist(),
                     "explained_variance": float(eigval[-1] / eigval.sum())})
        for r in range(len(atlas)):
            pcs.append({"group": group, "parcel_index": r, "parcel_id": atlas.parcel_id.iloc[r],
                        "network_id": atlas.network_id.iloc[r], "PC1": scores[r]})
    result = Result(tables={"regions": pd.DataFrame(rows), "subjects": pd.DataFrame(subject_rows),
                           "pca_scores": pd.DataFrame(pcs, columns=["group", "parcel_index", "parcel_id", "network_id", "PC1"]), "atlas": atlas},
                    arrays={"standardised": z, "reference_mean": mean, "reference_sd": sd, "valid_reference": valid},
                    metadata={"cohort": "adjusted", "pca": meta, "reference": f"{prepared.groups[0]}, per-region/atom population SD"})
    vectors = {f"{export_label(g)}_{ATOMS[a]}": groups[g][:, a, None] for g in prepared.groups[1:] for a in (0, 15)}
    vectors.update(parcel_index=np.arange(len(atlas))[:, None], parcel_id=atlas.parcel_id.astype(str).to_numpy(dtype=object)[:, None],
                   network_id=atlas.network_id.to_numpy()[:, None])
    for group in prepared.groups[1:]:
        sub = result.tables["pca_scores"].query("group == @group")
        vectors[f"{export_label(group)}_PC1"] = sub.PC1.to_numpy()[:, None] if len(sub) else np.full((len(atlas), 1), np.nan)
    result.matlab["regional_brain_exports"] = vectors
    result.figure("cn-reference.regions", "reference_scatter", table="regions")
    result.figure("cn-reference.subjects", "reference_scatter", table="subjects")
    for group in prepared.groups[1:]:
        result.figure(f"cn-reference.{group}.networks", "reference_scatter", table="regions", group=group)
        result.figure(f"cn-reference.{group}.pca", "pca", group=group)
    return prepared.finish(result)


def alterations(prepared):
    df, atlas, settings = context(prepared, 2)
    values = np.asarray(prepared.array("regional"))
    nets = sorted(atlas.network_id.unique())
    tables, corr_rows, arrays = [], [], {}
    for a, b in prepared.pairs:
        ia, ib = np.flatnonzero(df.analysis_group == a), np.flatnonzero(df.analysis_group == b)
        for atom in (0, 15):
            key = f"{ATOMS[atom]}_{b}_vs_{a}"
            rng = rng_for(settings.seed, f"{prepared.target('alterations')}.{key}")
            if min(len(ia), len(ib)) < 2:
                for net in nets:
                    tables.append({"atom": ATOMS[atom], "group_a": a, "group_b": b, "network_id": net,
                                   "observed": np.nan, "p": np.nan, "ci_lower": np.nan, "ci_upper": np.nan,
                                   "status": "insufficient_samples"})
                continue
            va, vb = values[ia, :, atom], values[ib, :, atom]
            baseline = va.mean(axis=0)
            delta = vb.mean(axis=0) - baseline
            masks = [atlas.network_id.to_numpy() == net for net in nets]
            obs = np.array([delta[m].mean() for m in masks])
            null = np.full((settings.permutations, len(nets)), np.nan)
            bootstrap = np.full((settings.resamples, len(nets)), np.nan)
            corr_null = np.full(settings.permutations, np.nan)
            combined = np.concatenate([va, vb])
            for rep in range(settings.permutations):
                shuffled = rng.permutation(delta)
                null[rep] = [shuffled[m].mean() for m in masks]
                perm = rng.permutation(len(combined))
                base = combined[perm[:len(va)]].mean(axis=0)
                change = combined[perm[len(va):]].mean(axis=0) - base
                corr_null[rep] = correlation(change, base)
            for rep in range(settings.resamples):
                change = vb[rng.integers(len(vb), size=len(vb))].mean(axis=0) - va[rng.integers(len(va), size=len(va))].mean(axis=0)
                bootstrap[rep] = [change[m].mean() for m in masks]
            for ni, net in enumerate(nets):
                lower, upper = interval(bootstrap[:, ni])
                tables.append({"atom": ATOMS[atom], "group_a": a, "group_b": b, "network_id": net,
                               "observed": obs[ni], "p": permutation_p(obs[ni], null[:, ni], centre=float(delta.mean())),
                               "ci_lower": lower, "ci_upper": upper, "n_a": len(va), "n_b": len(vb),
                               "status": "ok" if np.isfinite(obs[ni]) else "nonfinite_features"})
            corr_obs = correlation(delta, baseline)
            corr_rows.append({"atom": ATOMS[atom], "group_a": a, "group_b": b, "pearson_r": corr_obs,
                              "p": permutation_p(corr_obs, corr_null, centre=float(np.nanmean(corr_null))),
                              "status": "ok" if np.isfinite(corr_obs) else "constant_or_missing_features"})
            arrays[f"regional_shuffle_{key}"] = null
            arrays[f"subject_bootstrap_{key}"] = bootstrap
            permutation_name = "diagnosis_permutation" if prepared.grouping == "diagnosis" else "moca_permutation"
            arrays[f"{permutation_name}_{key}"] = corr_null
            arrays[f"delta_{key}"] = delta
    network = correct(pd.DataFrame(tables), ["atom", "group_a", "group_b"], settings)
    corrs = correct(pd.DataFrame(corr_rows), [], settings) if corr_rows else pd.DataFrame(columns=["atom", "group_a", "group_b", "pearson_r", "p", "q", "status"])
    result = Result(tables={"networks": network, "correlations": corrs, "atlas": atlas}, arrays=arrays,
                    metadata={"cohort": "adjusted", "exploratory": True,
                              "regional_null": "exchangeable regional changes; does not preserve spatial structure",
                              "correlation_null": "analysis-group-label permutation, centred at its empirical null mean; shared baseline retained",
                              "bootstrap": "stratified participant resampling; conditional on fitted confound adjustment"})
    for atom in ("rtr", "sts"):
        for a, b in prepared.pairs:
            result.figure(f"alterations.{atom}.{b}-vs-{a}", "network_intervals", atom=atom, group_a=a, group_b=b)
    return prepared.finish(result)


def tdmi(prepared):
    df, atlas = prepared.table("selected"), prepared.table("atlas")
    settings = prepared.config.statistics
    df = df.loc[df.analysis_group.isin(prepared.groups)].copy()
    if df.empty:
        raise InsufficientData(f"No selected scans with a target {prepared.grouping} label")
    df["value"] = prepared.array("raw_tdmi")[df.tensor_row.to_numpy(int)]
    observations = df[["scan_id", "subject", "session", "analysis_group", "n_time", "value"]].rename(columns={"analysis_group": "group"})
    summary = observations.groupby("group").value.agg(n="count", mean="mean", SD="std", median="median").reset_index()
    cn = summary.loc[summary.group == prepared.groups[0], "mean"]
    if len(cn):
        summary[f"difference_from_{prepared.groups[0]}"] = summary["mean"] - cn.iloc[0]
        summary[f"percent_difference_from_{prepared.groups[0]}"] = 100*(summary["mean"] / cn.iloc[0] - 1) if cn.iloc[0] != 0 else np.nan
    samples = [observations.loc[observations.group == g, "value"].to_numpy() for g in prepared.groups]
    tests = pd.DataFrame([{"method": "Welch" if w else "ordinary", "primary": w == settings.welch,
                           **omnibus(samples, w)} for w in (True, False)])
    pairs = correct(pd.DataFrame([{"group_a": a, "group_b": b, **compare(observations.loc[observations.group == a, "value"],
                          observations.loc[observations.group == b, "value"], settings.welch)} for a, b in prepared.pairs]), [], settings)
    result = Result(tables={"observations": observations, "summary": summary, "omnibus": tests, "pairwise": pairs,
                           "length_audit": prepared.table("inventory")},
                    metadata={"cohort": "selected eligible scans with analysis-group labels; no confound-completeness restriction",
                              "measure": "sum of 16 RAW atoms, mean across included edges", "tdmi_normalisation": "none"})
    result.figure("tdmi.groups", "group_distribution", table="observations", value="value")
    result.figure("tdmi.lengths", "lengths")
    return prepared.finish(result)
