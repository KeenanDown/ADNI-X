from __future__ import annotations

import numpy as np
import pandas as pd

from .constants import ATOMS, GROUPS, NETWORKS
from .numerics import normalise_tdmi, network_means
from .results import Result, InsufficientData


def longitudinal(prepared, ratios=False):
    df = prepared.table("visits")
    atlas = prepared.table("atlas")
    t = prepared.config.transform
    networks = sorted(atlas.network_id.unique())
    raw_values = prepared.array("raw_network")
    values = raw_values
    invalid_edges = []
    if ratios:
        raw = prepared.array("raw")
        values = np.full(raw_values.shape, np.nan)
        for i in range(len(df)):
            normalised, invalid = normalise_tdmi(raw[i], t.epsilon, t.exclude_diagonal)
            values[i] = network_means(normalised, atlas.network_id, t.exclude_diagonal, t.network_edges)
            invalid_edges.append({"scan_id": df.scan_id.iloc[i], "invalid_tdmi_edges": invalid})
    dates = pd.to_datetime(df.EXAMDATE, errors="coerce")
    first = dates.groupby(df.subject).transform("min")
    elapsed = (dates - first).dt.days / 365.25
    rows = []
    for i, row in df.iterrows():
        if row.group not in GROUPS:
            continue
        for ni, net in enumerate(networks):
            measures = {"rtr": values[i, ni, 0], "sts": values[i, ni, 15]}
            if not ratios:
                measures["TDMI"] = raw_values[i, ni].sum()
            for measure, value in measures.items():
                rows.append({"scan_id": row.scan_id, "subject": row.subject, "session": row.session,
                             "group": row.group, "network_id": net, "network": NETWORKS[net],
                             "measure": measure, "value": value, "date": row.EXAMDATE,
                             "calendar_year": dates.iloc[i].year if pd.notna(dates.iloc[i]) else np.nan,
                             "years_since_first_eligible_scan": elapsed.iloc[i], "years_since_ADNI_baseline": row.Years_bl,
                             "visit_age": row.visit_age})
    if not rows:
        raise InsufficientData("No eligible longitudinal visits with diagnosis")
    visits = pd.DataFrame(rows)
    bins = []
    axes = {"calendar_year": 1., "years_since_first_eligible_scan": 1., "years_since_ADNI_baseline": 1., "visit_age": 5.}
    for axis, width in axes.items():
        data = visits.loc[np.isfinite(visits[axis])].copy()
        data["bin_start"] = np.floor(data[axis] / width) * width
        # Exact maximum values and a cohort containing only baseline visits are included.
        for keys, sub in data.groupby(["measure", "network_id", "group", "bin_start"], observed=True):
            measure, net, group, start = keys
            finite = sub.loc[np.isfinite(sub.value)]
            # Each participant contributes once per bin, avoiding extra weight from repeated scans.
            person_values = finite.groupby("subject").value.mean()
            bins.append({"axis": axis, "measure": measure, "network_id": net, "group": group,
                         "bin_start": start, "bin_end": start + width, "x": start + width / 2,
                         "mean": person_values.mean() if len(person_values) else np.nan,
                         "n_visits": len(finite), "n_subjects": len(person_values),
                         "n_missing_visits": len(sub) - len(finite),
                         "status": "ok" if len(person_values) else "no_finite_values"})
    columns = ["axis", "measure", "network_id", "group", "bin_start", "bin_end", "x", "mean", "n_visits", "n_subjects", "n_missing_visits", "status"]
    binned = pd.DataFrame(bins, columns=columns)
    # Preserve the notebook's exploratory 5-point rolling curve, on participant/bin means.
    rolling = binned.loc[binned.axis == "years_since_first_eligible_scan"].copy()
    rolling = rolling.sort_values("x", kind="stable")
    rolling["mean"] = rolling.groupby(["measure", "network_id", "group"], observed=True)["mean"].transform(lambda x: x.rolling(5, min_periods=2).mean())
    result = Result(tables={"visits": visits, "binned": binned, "rolling": rolling,
                           "normalisation_audit": pd.DataFrame(invalid_edges, columns=["scan_id", "invalid_tdmi_edges"])},
                    metadata={"cohort": "all eligible visits, true visit diagnosis", "descriptive": True,
                              "adjustment": "none", "tdmi_normalisation": "tdmi" if ratios else "none",
                              "aggregation": "mean within participant/bin, then mean across participants; no inferential CI",
                              "rolling": "5 adjacent elapsed-time bin means, minimum 2; descriptive only"})
    name = "longitudinal-tdmi-ratio" if ratios else "longitudinal"
    for measure in (["rtr", "sts"] if ratios else ["rtr", "sts", "TDMI"]):
        for axis in axes:
            result.figure(f"{name}.{measure}.{axis}", "longitudinal", table="binned", axis=axis, measure=measure)
        result.figure(f"{name}.{measure}.rolling", "longitudinal", table="rolling", axis="years_since_first_eligible_scan", measure=measure)
    return result
