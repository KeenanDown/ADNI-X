from __future__ import annotations

from itertools import combinations
import json
import warnings
import numpy as np
import pandas as pd
from scipy import stats

from .constants import ATOMS, NETWORKS, TAXA_FOUR, TAXA_SIX
from .numerics import (compare, fdr, omnibus, hierarchy, correlation, cn_standardise,
                       interval, permutation_p, connection_mask)
from .provenance import rng_for
from .results import Result, InsufficientData


def context(prepared, minimum=1):
    df, atlas = prepared.table("subjects"), prepared.table("atlas")
    if len(df) < minimum:
        raise InsufficientData(f"Need at least {minimum} adjusted participants; found {len(df)}")
    return df, atlas, prepared.config.statistics


def correct(table, family_columns, settings):
    table = table.copy()
    table["q"], table["significant"] = np.nan, False
    groups = table.groupby(family_columns, dropna=False).groups if family_columns else {"all": table.index}
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


def group_matrices(prepared):
    df, atlas, settings = context(prepared, 2)
    data = prepared.array("adjusted")
    p = len(atlas)
    means = np.full((3, p, p, 16), np.nan)
    for gi, group in enumerate(prepared.groups):
        indices = np.flatnonzero(df.analysis_group == group)
        if len(indices):
            total, count = np.zeros((p, p, 16)), np.zeros((p, p, 16), int)
            for i in indices:
                valid = np.isfinite(data[i])
                total += np.where(valid, data[i], 0)
                count += valid
            np.divide(total, count, out=means[gi], where=count > 0)
    mask = connection_mask(p, prepared.config.transform.exclude_diagonal)
    ii, jj = np.where(mask)
    flat = data.reshape(len(data), -1)
    rows = []
    for atom in (0, 15):
        columns = (ii*p + jj)*16 + atom
        for a, b in prepared.pairs:
            ia, ib = np.flatnonzero(df.analysis_group == a), np.flatnonzero(df.analysis_group == b)
            for start in range(0, len(columns), prepared.config.transform.block_features):
                end = start + prepared.config.transform.block_features
                cols = columns[start:end]
                comparisons = _compare_columns(flat[np.ix_(ia, cols)], flat[np.ix_(ib, cols)], settings.welch)
                rows.append(pd.DataFrame({"atom": ATOMS[atom], "group_a": a, "group_b": b,
                    "parcel_i": ii[start:end], "parcel_j": jj[start:end], **comparisons}))
    table = correct(pd.concat(rows, ignore_index=True), ["atom", "group_a", "group_b"], settings)
    result = Result(tables={"edges": table, "atlas": atlas,
                           "group_counts": df.groupby("analysis_group").size().rename_axis("group").reset_index(name="n")},
                    arrays={"group_means": means}, metadata={"groups": prepared.groups, "cohort": "adjusted", "family": "edges per atom/contrast"})
    for atom in (0, 15):
        result.figure(f"group-matrices.{ATOMS[atom]}", "group_matrices", atom=atom, comparison=f"{prepared.groups[-1]}-vs-{prepared.groups[0]}")
    return prepared.finish(result)


def coefficients(prepared):
    _, atlas, _ = context(prepared, 2)
    model = json.loads((prepared.model_path / "design.json").read_text())
    beta = np.asarray(prepared.array("coefficients"))
    names = ["Intercept", *model["columns"]]
    result = Result(tables={"atlas": atlas, "design": pd.DataFrame(prepared.array("design_matrix"), columns=names)},
                    arrays={"coefficients": beta.reshape(len(beta), len(atlas), len(atlas), 16),
                            "nobs": np.asarray(prepared.array("coefficient_nobs"))},
                    metadata={"design": model, "names": names, "cohort": "adjusted"})
    for atom in (0, 15):
        result.figure(f"coefficients.{ATOMS[atom]}", "coefficients", atom=atom)
    return prepared.finish(result)


def network_effects(prepared):
    df, atlas, settings = context(prepared, 2)
    values = prepared.array("network")
    nets = sorted(atlas.network_id.unique())
    rows, observations = [], []
    for ni, net in enumerate(nets):
        for atom, atom_name in enumerate(ATOMS):
            for i, row in df.iterrows():
                observations.append({"scan_id": row.scan_id, "subject": row.subject, "group": row.analysis_group,
                                     "network_id": net, "network": NETWORKS[net], "atom": atom_name, "value": values[i, ni, atom]})
            for a, b in prepared.pairs:
                rows.append({"network_id": net, "network": NETWORKS[net], "atom": atom_name, "group_a": a, "group_b": b,
                             **compare(values[df.analysis_group == a, ni, atom], values[df.analysis_group == b, ni, atom], settings.welch)})
    table = correct(pd.DataFrame(rows), ["group_a", "group_b"], settings)
    result = Result(tables={"effects": table, "observations": pd.DataFrame(observations)},
                    metadata={"cohort": "adjusted", "family": "all 16 atoms x all networks per contrast"})
    for a, b in prepared.pairs:
        for effect in ("cohens_d", "hedges_g"):
            for subset in ("all-atoms", "rtr-sts"):
                result.figure(f"network-effects.{b}-vs-{a}.{effect}.{subset}", "effects_grid",
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


def hierarchies(prepared):
    df, atlas, settings = context(prepared, 2)
    values = np.asarray(prepared.array("regional"))
    names = list(prepared.groups)
    labels = df.analysis_group
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
        centres = profile.groupby("parcel_index")[atom].transform("mean")
        profile[f"{atom}_deviation"] = profile[atom] - centres
        endpoints = profile.loc[profile.group.isin([names[0], names[-1]])].groupby("parcel_index")[atom].mean()
        profile[f"{atom}_midpoint"] = profile[atom] - profile.parcel_index.map(endpoints)
    result = Result(tables={"profiles": profile, "correlations": pd.DataFrame(correlations), "atlas": atlas},
                    arrays={"bootstrap_correlations": samples}, metadata={"groups": names, "pairs": pairs, "cohort": "adjusted", "grouping": prepared.grouping})
    for measure in ("hierarchy", "rtr", "sts", "rtr_deviation", "sts_deviation", "rtr_midpoint", "sts_midpoint"):
        result.figure(f"{name}.{measure}", "profiles", measure=measure, anchor=names[0],
                      order="difference" if "_" in measure else "reference")
    result.figure(f"{name}.correlations", "correlation_intervals")
    return prepared.finish(result)


def taxa(prepared, six=False):
    df, atlas, settings = context(prepared, 2)
    values = np.asarray(prepared.array("network"))
    taxonomy = TAXA_SIX if six else TAXA_FOUR
    groups = list(prepared.groups)
    nets = sorted(atlas.network_id.unique())
    observations, tests, posthoc = [], [], []
    for taxon, atoms in taxonomy.items():
        vals = values[:, :, list(atoms)].sum(axis=2)
        for ni, net in enumerate(nets):
            samples = [vals[df.analysis_group == g, ni] for g in groups]
            for i, subject in df.iterrows():
                observations.append({"scan_id": subject.scan_id, "subject": subject.subject,
                                     "group": subject.analysis_group, "taxon": taxon, "network_id": net,
                                     "network": NETWORKS[net], "value": vals[i, ni]})
            tests.append({"taxon": taxon, "network_id": net, **omnibus(samples, settings.welch)})
            for ai, bi in combinations(range(len(groups)), 2):
                posthoc.append({"taxon": taxon, "network_id": net, "group_a": groups[ai], "group_b": groups[bi],
                                **compare(samples[ai], samples[bi], settings.welch)})
    tests = correct(pd.DataFrame(tests), [], settings)
    posthoc = correct(pd.DataFrame(posthoc), [], settings)
    posthoc = posthoc.merge(tests[["taxon", "network_id", "q"]].rename(columns={"q": "omnibus_q"}), on=["taxon", "network_id"], validate="many_to_one")
    name = f"taxa-{'six' if six else 'four'}-dx"
    result = Result(tables={"observations": pd.DataFrame(observations), "omnibus": tests, "posthoc": posthoc},
                    metadata={"groups": groups, "cohort": "adjusted", "taxonomy": taxonomy,
                              "family": "all taxon/network panels; all panel/contrast posthoc tests separately"})
    result.figure(f"{name}.panels", "taxa_panels")
    return prepared.finish(result)


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
