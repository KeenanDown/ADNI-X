"""Scientific invariants for the HPC rerun, using synthetic participants only."""
from dataclasses import replace
from itertools import combinations
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from scipy import stats

from adni_x.analyses import coefficients, group_matrices, network_effects, taxa
from adni_x.constants import ATOMS, MOCA_GROUPS, TAXA_PRINCIPAL
from adni_x.inference import hc3_coefficients, validate_symmetric_atoms
from adni_x.inputs import clinical_table
from adni_x.numerics import connection_mask, fdr
from adni_x.plotting import Bundle, render, plot
from adni_x.preparation import prepare
from adni_x.provenance import file_hash
from adni_x.runner import run


class Proxy:
    def __init__(self, prepared, subjects=None, adjusted=None):
        self.prepared, self.subjects, self.adjusted = prepared, subjects, adjusted

    def __getattr__(self, name):
        return getattr(self.prepared, name)

    def table(self, name):
        return self.subjects if name == "subjects" and self.subjects is not None else self.prepared.table(name)

    def array(self, name):
        return self.adjusted if name == "adjusted" and self.adjusted is not None else self.prepared.array(name)


@pytest.mark.parametrize("p,count", [(232, 26796), (116, 6670)])
def test_unique_pair_counts(p, count):
    assert connection_mask(p, upper=True).sum() == count
    assert connection_mask(p, exclude_diagonal=False, upper=True).sum() == count+p


def test_symmetry_checks_participants_and_missingness(prepared):
    data = np.array(prepared.array("adjusted"))
    check = validate_symmetric_atoms(data, np.arange(len(data)), block_features=2)
    assert check["atoms"]["sts"]["unique_pairs"] == 6
    # Opposite participant errors cancel in the group mean but must still fail.
    data[0, 0, 1, 15] += .2
    data[1, 0, 1, 15] -= .2
    with pytest.raises(ValueError, match="sts symmetry failed"):
        validate_symmetric_atoms(data, np.arange(len(data)))
    data = np.array(prepared.array("adjusted")); data[0, 0, 1, 0] = np.nan
    with pytest.raises(ValueError, match="finite-value masks"):
        validate_symmetric_atoms(data, np.arange(len(data)))


@pytest.mark.parametrize("grouping", ["group", "MOCA_category"])
def test_unique_matrix_families_and_values(prepared, grouping):
    result = group_matrices(prepared, grouping)
    edges = result.tables["edges"]
    assert len(edges) == 2*3*6 and (edges.parcel_i < edges.parcel_j).all()
    df = prepared.table("subjects"); data = prepared.array("adjusted")
    for _, sub in edges.groupby(["atom", "group_a", "group_b"]):
        np.testing.assert_allclose(sub.q, fdr(sub.p)[0])
    row = edges.iloc[0]
    a, b = df.index[df[grouping] == row.group_a], df.index[df[grouping] == row.group_b]
    expected = stats.ttest_ind(data[b, row.parcel_i, row.parcel_j, ATOMS.index(row.atom)],
                              data[a, row.parcel_i, row.parcel_j, ATOMS.index(row.atom)], equal_var=False)
    np.testing.assert_allclose(row.p, expected.pvalue)


@pytest.mark.parametrize("grouping", ["group", "MOCA_category"])
def test_principal_modes_preserve_each_atom_and_families(prepared, grouping):
    result = taxa(prepared, grouping=grouping, principal=True)
    assert result.metadata["display_order"] == list(TAXA_PRINCIPAL)
    observations = result.tables["observations"]
    df = prepared.table("subjects"); values = prepared.array("network")
    for name, atoms in TAXA_PRINCIPAL.items():
        sub = observations.query("taxon == @name and network_id == 1")
        expected = values[:, 0, list(atoms)].sum(axis=1)
        np.testing.assert_allclose(sub.value, expected)
    assert len(result.tables["omnibus"]) == 4*2
    assert len(result.tables["posthoc"]) == 4*2*3
    for key in ("omnibus", "posthoc"):
        table = result.tables[key]
        np.testing.assert_allclose(table.q, fdr(table.p)[0])


def test_moca_missingness_keeps_original_participant_indices(prepared):
    df = prepared.table("subjects"); df.loc[[0, 3, 15], "MOCA_category"] = np.nan
    proxy = Proxy(prepared, subjects=df)
    assert len(group_matrices(proxy, "MOCA_category").tables["group_counts"]) == 3
    result = taxa(proxy, grouping="MOCA_category", principal=True)
    assert result.tables["observations"].subject.nunique() == len(df)-3
    net = network_effects(proxy, "MOCA_category")
    assert net.tables["observations"].subject.nunique() == len(df)-3
    assert net.metadata["groups"] == MOCA_GROUPS
    row = result.tables["observations"].iloc[0]
    participant = df.index[df.subject == row.subject][0]
    np.testing.assert_allclose(row.value, prepared.array("network")[participant, 0, [7, 11, 3]].sum())
    assert taxa(proxy, principal=True).tables["observations"].subject.nunique() == len(df)


def test_moca_boundaries_and_noninteger_scores(tmp_path):
    scores = [0, 19, 20, 28, 29, 30, -1, 31, 20.5, np.nan]
    pd.DataFrame({"PTID": [f"001_S_{i:04d}" for i in range(len(scores))],
                  "VISCODE": "bl", "MOCA": scores}).to_csv(tmp_path/"clinical.csv", index=False)
    table, _ = clinical_table(tmp_path/"clinical.csv")
    assert list(table.MOCA_category.iloc[:6]) == ["0-19", "0-19", "20-28", "20-28", "29-30", "29-30"]
    assert table.MOCA_category.iloc[6:].isna().all()


def test_hc3_matches_independent_full_sandwich():
    rng = np.random.default_rng(12)
    X = np.column_stack([np.ones(40), rng.normal(size=(40, 3))])
    Y = X@rng.normal(size=(4, 3))+rng.normal(size=(40, 3))*(1+abs(X[:, 1, None]))
    Y[2:6, 1] = np.nan
    beta = np.column_stack([np.linalg.lstsq(X[np.isfinite(y)], y[np.isfinite(y)], rcond=None)[0] for y in Y.T])
    result = hc3_coefficients(X, Y, beta)
    for j, y in enumerate(Y.T):
        valid = np.isfinite(y); x, y = X[valid], y[valid]
        bread = np.linalg.inv(x.T@x)
        hat = x@bread@x.T
        e = y-x@beta[:, j]
        meat = x.T@np.diag(e**2/(1-np.diag(hat))**2)@x
        se = np.sqrt(np.diag(bread@meat@bread))
        df = len(y)-x.shape[1]
        np.testing.assert_allclose(result["standard_error"][:, j], se, rtol=1e-10)
        np.testing.assert_allclose(result["p"][:, j], 2*stats.t.sf(abs(beta[:, j]/se), df))
        np.testing.assert_allclose(result["ci_lower"][:, j], beta[:, j]-stats.t.ppf(.975, df)*se)
        assert result["df_resid"][j] == df and result["nobs"][j] == len(y)


@pytest.mark.parametrize("kind,status", [("rank", "rank_deficient_design"), ("small", "insufficient_observations"),
                                        ("leverage", "unit_leverage"), ("zero", "zero_standard_error")])
def test_hc3_failure_states(kind, status):
    X = np.column_stack([np.ones(6), np.arange(6)])
    Y = np.arange(6)[:, None].astype(float)
    if kind == "rank": X[:, 1] = 1
    if kind == "small": X, Y = X[:2], Y[:2]
    if kind == "leverage": X[:, 1] = [1, 0, 0, 0, 0, 0]
    if kind == "zero": Y[:] = 0
    beta = np.linalg.lstsq(X, Y, rcond=None)[0]
    result = hc3_coefficients(X, Y, beta)
    assert (result["status"] == status).all()
    assert np.isnan(result["p"]).all()


@pytest.mark.parametrize("normalisation", ["none", "tdmi"])
def test_coefficient_inference_uses_original_responses(config, normalisation):
    config = replace(config, transform=replace(config.transform, tdmi=normalisation))
    prepared = prepare(config); result = coefficients(prepared)
    table = result.tables["inference"]
    assert len(table) == 2*6*4 and "Intercept" not in set(table.coefficient)
    assert set(table.status) == {"ok"}
    for _, sub in table.groupby("atom"):
        np.testing.assert_allclose(sub.q, fdr(sub.p)[0])
    df = prepared.table("subjects"); X = prepared.array("design_matrix")
    raw = prepared.array("raw")[df.tensor_row.to_numpy(int), 0, 1, :]
    Y = raw[:, [0]] if normalisation == "none" else raw[:, [0]]/raw.sum(axis=1, keepdims=True)
    beta = prepared.array("coefficients")[:, [(0*4+1)*16]]
    expected = hc3_coefficients(X, Y, beta)
    row = table.query("atom == 'rtr' and parcel_i == 0 and parcel_j == 1").iloc[0]
    np.testing.assert_allclose(row.standard_error, expected["standard_error"][1, 0])
    np.testing.assert_allclose(row.p, expected["p"][1, 0])


def test_intercept_only_coefficients_have_empty_inference(config):
    prepared = prepare(replace(config, transform=replace(config.transform, confounds=())))
    assert coefficients(prepared).tables["inference"].empty


def test_default_targets_are_moca_and_mmse_is_opt_in(prepared, tmp_path):
    # Existing runs can already contain explicitly requested optional targets.
    # Verify default selection in a fresh result root, without deleting history.
    config = replace(prepared.config, results=str(tmp_path/"default-results"))
    report = run("all", config, render=False)
    assert all(v["status"] == "complete" for v in report["targets"].values())
    assert "taxa-principal-moca" in report["targets"] and "group-matrices-moca" in report["targets"]
    assert "taxa-six-mmse" not in report["targets"]
    assert run("taxa-six-mmse", config, render=False)["targets"]["taxa-six-mmse"]["status"] == "complete"


@pytest.mark.parametrize("grouping", ["group", "MOCA_category"])
def test_profile_baseline_order_and_old_bundle_refresh(bundles, grouping):
    name = "hierarchies" if grouping == "group" else "hierarchies-moca"
    bundle = Bundle(bundles/name)
    recipe = next(r for r in json.loads((bundle.path/"figures.json").read_text()) if r["measure"] == "sts_deviation")
    baseline = bundle.details["groups"][0]
    df = bundle.table("profiles"); pivot = df.pivot(index="parcel_index", columns="group", values="sts")
    order = (pivot[bundle.details["groups"][-1]]-pivot[baseline]).sort_values(ascending=False, kind="stable").index
    # Model an old bundle whose stored deviation used the three-group mean.
    df.sts_deviation = df.sts-df.groupby("parcel_index").sts.transform("mean")
    original_table = bundle.table
    bundle.table = lambda name: df.copy() if name == "profiles" else original_table(name)
    fig = render(bundle, recipe)
    try:
        main = fig.axes[0]
        np.testing.assert_allclose(main.collections[0].get_offsets()[:, 1], 0)
        expected = pivot.loc[order, bundle.details["groups"][-1]]-pivot.loc[order, baseline]
        np.testing.assert_allclose(main.collections[2].get_offsets()[:, 1], expected)
        assert all(ax.get_ylim() == main.get_ylim() for ax in fig.axes[2:])
        assert f"difference from {baseline}" in main.get_ylabel()
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        note = next(t for t in fig.texts if t.get_text().startswith("Descriptive profiles"))
        assert not note.get_window_extent(renderer).overlaps(fig.legends[0].get_window_extent(renderer))
    finally: plt.close(fig)


@pytest.mark.parametrize("subset", ["all-atoms", "rtr-sts"])
def test_mci_cn_cohens_d_export_and_full_family(bundles, subset):
    bundle = Bundle(bundles/"network-effects")
    figure_id = f"network-effects.MCI-vs-CN.cohens_d.{subset}"
    recipe = next(r for r in json.loads((bundle.path/"figures.json").read_text()) if r["id"] == figure_id)
    df = bundle.table("effects").query("group_a == 'CN' and group_b == 'MCI'")
    np.testing.assert_allclose(df.q, fdr(df.p)[0])
    expected = df.pivot(index="atom", columns="network_id", values="cohens_d").reindex(recipe["atoms"])
    flags = df.pivot(index="atom", columns="network_id", values="significant").reindex_like(expected)
    fig = render(bundle, recipe)
    try:
        shown = np.where(flags, expected, np.nan)
        if np.isfinite(shown).any():
            np.testing.assert_allclose(fig.axes[0].images[0].get_array().filled(np.nan), shown, equal_nan=True)
        assert "MCI − CN" in fig.axes[0].get_title()
    finally: plt.close(fig)
    exported = plot(figure_id, bundle.path)[0]
    assert {Path(p).suffix for p in exported["files"]} == {".svg", ".png", ".pdf"}


def test_matrix_and_coefficient_significance_masks_are_mirrored(bundles):
    for target, figure in [("group-matrices", "group-matrices.rtr"), ("coefficients", "coefficients.rtr.significant")]:
        bundle = Bundle(bundles/target)
        recipe = next(r for r in json.loads((bundle.path/"figures.json").read_text()) if r["id"] == figure)
        fig = render(bundle, recipe)
        try:
            images = [im for ax in fig.axes for im in ax.images]
            for im in images:
                a = im.get_array().filled(np.nan)
                np.testing.assert_allclose(a, a.T, equal_nan=True)
        finally: plt.close(fig)


def test_refreshed_profiles_preserve_bundle_files(bundles):
    bundle = bundles/"hierarchies"
    before = {p.name: file_hash(p) for p in bundle.iterdir() if p.is_file()}
    plot("hierarchies.sts_deviation", bundle)
    assert before == {p.name: file_hash(p) for p in bundle.iterdir() if p.is_file()}
