import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

from adni_x.catalogue import CATALOGUE
from adni_x.provenance import rng_for, verify
from adni_x.runner import run


@pytest.mark.parametrize("spec", CATALOGUE, ids=lambda spec: spec.id)
def test_every_analysis_has_independent_numerical_outputs(prepared, spec):
    result = spec.function(prepared)
    assert result.tables or result.arrays
    assert {r["id"] for r in result.figures} == set(spec.figures)
    for array in result.arrays.values():
        assert np.asarray(array).dtype.kind != "O"


def test_prediction_participant_splits(prepared):
    from adni_x.prediction import lda, heldout_hierarchy
    for function in (lda, heldout_hierarchy):
        result = function(prepared)
        assert not set(result.metadata["fit_subjects"]) & set(result.metadata["test_subjects"])
        table = result.tables["predictions" if function is lda else "scores"]
        assert not table.subject.duplicated().any()


def test_fdr_does_not_depend_on_figure_subset(prepared):
    from adni_x.analyses import network_effects
    from adni_x.numerics import fdr
    result = network_effects(prepared)
    table = result.tables["effects"]
    for _, sub in table.groupby(["group_a", "group_b"]):
        expected, _ = fdr(sub.p)
        np.testing.assert_allclose(sub.q, expected)
        assert len(sub) == 16*2


def test_deterministic_bootstrap_and_permutation(prepared):
    from adni_x.analyses import hierarchies, alterations
    for function in (hierarchies, alterations):
        one, two = function(prepared), function(prepared)
        for key in one.arrays:
            np.testing.assert_array_equal(one.arrays[key], two.arrays[key])


def test_alone_and_all_identical_and_complete(prepared):
    one = run("hierarchies", prepared.config, render=False)
    bundle = Path(one["path"]) / "hierarchies"
    old = np.load(bundle / "bootstrap_correlations.npy")
    all_results = run("all", prepared.config, render=False, include_sensitivity=True)
    assert all(v["status"] == "complete" for v in all_results["targets"].values()), all_results["targets"]
    np.testing.assert_array_equal(old, np.load(bundle / "bootstrap_correlations.npy"))
    assert len(all_results["targets"]) == len(CATALOGUE)
    for target in all_results["targets"].values():
        assert verify(target["path"])


def test_longitudinal_includes_maximum_time_and_counts_people(prepared):
    from adni_x.longitudinal import longitudinal
    result = longitudinal(prepared)
    visits, binned = result.tables["visits"], result.tables["binned"]
    elapsed = binned.query("axis == 'years_since_first_eligible_scan'")
    assert elapsed.bin_start.max() >= 1
    assert (binned.n_subjects <= binned.n_visits).all()
    assert visits.scan_id.nunique() == 48
    assert "years_since_ADNI_baseline" in visits


def test_render_without_sources_or_preparation(prepared, tmp_path):
    import shutil
    from adni_x.plotting import plot
    report = run("group-matrices", prepared.config, render=False)
    copied = tmp_path / "portable-result"
    shutil.copytree(Path(report["path"]) / "group-matrices", copied)
    # No config is loaded and source/preparation locations need not exist.
    manifest_path = copied / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["configuration"]["inputs"] = {"tensors": "/does/not/exist"}
    manifest["sources"] = {"/does/not/exist": {"sha256": "absent"}}
    manifest_path.write_text(json.dumps(manifest))
    output = plot("group-matrices.rtr", copied)
    assert len(output) == 1
    for file in output[0]["files"]:
        assert Path(file).stat().st_size > 1000
        assert "tdmi-none" in Path(file).name and "latest" in Path(file).name
