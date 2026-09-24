"""Diagnosis and cognition are independent groupings, never relabelled aliases."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import pytest
from scipy.io import loadmat, savemat

from adni_x.analyses import network_effects, cn_reference
from adni_x.catalogue import ANALYSES, CATALOGUE
from adni_x.inputs import clinical_table
from adni_x.numerics import compare, Design, fit_coefficients, fdr
from adni_x.preparation import prepare, target_id
from adni_x.prediction import lda, heldout_hierarchy
from adni_x.runner import run


def test_catalogue_covers_both_groupings_without_mmse():
    dx = [a for a in CATALOGUE if a.grouping == "diagnosis"]
    assert len(dx) == 14
    assert {target_id(a.id, "moca") for a in dx} == {a.id for a in CATALOGUE if a.grouping == "moca"}
    assert not any("mmse" in a.id for a in CATALOGUE)
    assert {a.id for a in CATALOGUE if a.opt_in} == {"longitudinal-tdmi-ratio", "longitudinal-tdmi-ratio-moca"}


def test_moca_boundaries_and_mmse_not_categorised(config):
    df = pd.read_csv(config.inputs.clinical)
    df["MOCA"] = df.MOCA.astype(float)
    scores = [-1, 0, 19, 20, 28, 29, 30, 31, 28.5, np.nan]
    df.loc[:len(scores)-1, "MOCA"] = scores
    df.to_csv(config.inputs.clinical, index=False)
    clinical, _ = clinical_table(config.inputs.clinical)
    assert "MMSE_category" not in clinical
    actual = clinical.MOCA_category.iloc[:len(scores)].fillna("missing").tolist()
    assert actual == ["missing", "0-19", "0-19", "20-28", "20-28", "29-30", "29-30", "missing", "missing", "missing"]


def test_independent_cohorts_and_fits_with_discordant_labels(config):
    df = pd.read_csv(config.inputs.clinical)
    df["MOCA"] = np.tile(np.repeat([15, 30, 25], 2), 8)  # Each category crosses diagnoses.
    df.loc[1, "MOCA"] = np.nan  # baseline remains valid; must not substitute it
    df.loc[3, "DX"] = np.nan  # usable MoCA despite missing diagnosis
    df.to_csv(config.inputs.clinical, index=False)
    dataset = prepare(config)
    dx, moca = dataset.table("subjects"), dataset.for_group("moca").table("subjects")
    assert "ADNI001S0001" in set(dx.subject) - set(moca.subject)
    assert "ADNI001S0002" in set(moca.subject) - set(dx.subject)
    assert (moca.month == 12).all()
    assert (moca.group != moca.analysis_group).all()
    assert dataset.array("raw").filename == dataset.for_group("moca").array("raw").filename
    for grouping in ("diagnosis", "moca"):
        view = dataset.for_group(grouping)
        subjects = view.table("subjects")
        raw = dataset.array("raw")[subjects.tensor_row.to_numpy(int)]
        np.testing.assert_allclose(view.array("adjusted").mean(axis=0), raw.mean(axis=0), atol=1e-9)
        design = Design.fit(subjects, config.transform)
        beta, _ = fit_coefficients(design.matrix(subjects), raw.reshape(len(subjects), -1))
        np.testing.assert_allclose(view.array("coefficients"), beta, atol=1e-9)
    result = network_effects(dataset.for_group("moca"))
    observations = result.tables["observations"]
    assert set(observations.group) == {"29-30", "20-28", "0-19"}
    assert set(observations.diagnosis.dropna()) == {"CN", "MCI", "AD"}
    values = dataset.for_group("moca").array("network")
    row = result.tables["effects"].iloc[0]
    expected = compare(values[moca.analysis_group == row.group_a, 0, 0],
                       values[moca.analysis_group == row.group_b, 0, 0], True)
    assert row.difference == pytest.approx(expected["difference"])
    for _, panel in result.tables["effects"].groupby(["group_a", "group_b"]):
        np.testing.assert_allclose(panel.q, fdr(panel.p)[0])


def test_missing_moca_only_disables_moca_analyses(config):
    df = pd.read_csv(config.inputs.clinical)
    df.MOCA = np.nan
    df.to_csv(config.inputs.clinical, index=False)
    report = run("all", config, render=False)
    for spec in CATALOGUE:
        if spec.opt_in:
            continue
        assert report["targets"][spec.id]["status"] == ("complete" if spec.grouping == "diagnosis" else "insufficient_data")
    with pytest.raises(ValueError, match="Unknown analysis"):
        run("taxa-four-mmse", config, render=False)


def test_moca_reference_and_matlab_labels(prepared):
    view = prepared.for_group("moca")
    result = cn_reference(view)
    subjects = view.table("subjects")
    expected = view.array("regional")[subjects.analysis_group == "29-30"].mean(axis=0)
    np.testing.assert_allclose(result.arrays["reference_mean"], expected)
    vectors = result.matlab["regional_brain_exports"]
    assert "MOCA_20_28_rtr" in vectors and "MOCA_0_19_PC1" in vectors
    assert result.metadata["reference_group"] == "29-30"


@pytest.mark.parametrize("function", [lda, heldout_hierarchy])
def test_moca_prediction_training_fit_unaffected_by_heldout_changes(config, function):
    dataset = prepare(config).for_group("moca")
    original = function(dataset)
    assert not set(original.metadata["fit_subjects"]) & set(original.metadata["test_subjects"])
    person = original.metadata["test_subjects"][0]
    subject = dataset.table("subjects").set_index("subject").loc[person]
    tensor = loadmat(subject.tensor_file)["PhiIDTensor"]
    savemat(subject.tensor_file, {"PhiIDTensor": tensor + 123})
    df = pd.read_csv(config.inputs.clinical)
    df.loc[(df.PTID.str.replace("_", "").radd("ADNI") == person) & (df.VISCODE == "m12"), "ICV"] = 1e9
    df.to_csv(config.inputs.clinical, index=False)
    changed = function(prepare(config).for_group("moca"))
    assert original.metadata["confound_design"] == changed.metadata["confound_design"]
    np.testing.assert_array_equal(original.arrays["confound_coefficients"], changed.arrays["confound_coefficients"])
    if function is heldout_hierarchy:
        np.testing.assert_array_equal(original.arrays["reference_training_hierarchy"], changed.arrays["reference_training_hierarchy"])
    else:
        assert original.metadata["confusion_labels"] == ["29-30", "0-19"]
        assert original.metadata["probability_class"] == "0-19"
        np.testing.assert_array_equal(original.arrays["lda_coefficients"], changed.arrays["lda_coefficients"])


def test_moca_alone_and_all_numeric_outputs_identical(prepared, tmp_path):
    from dataclasses import replace
    standalone = run("hierarchies-moca", replace(prepared.config, results=str(tmp_path / "alone")), render=False)
    combined = run("all", replace(prepared.config, results=str(tmp_path / "all")), render=False)
    a = Path(standalone["path"]) / "hierarchies-moca"
    b = Path(combined["path"]) / "hierarchies-moca"
    for file in a.glob("*.npy"):
        np.testing.assert_array_equal(np.load(file), np.load(b / file.name))
    for name in ("profiles", "correlations"):
        pd.testing.assert_frame_equal(pd.read_csv(a / f"{name}.csv"), pd.read_csv(b / f"{name}.csv"))
    cohort = pd.read_csv(a / "cohort_adjusted.csv")
    assert set(cohort.group) == {"CN", "MCI", "AD"}
    assert set(cohort.analysis_group) == {"29-30", "20-28", "0-19"}
    assert json.loads((a / "manifest.json").read_text())["grouping"] == "moca"


def test_longitudinal_moca_uses_visit_score_and_original_time_origin(config):
    from adni_x.longitudinal import longitudinal
    clinical = pd.read_csv(config.inputs.clinical)
    clinical.loc[0, "MOCA"] = np.nan
    clinical.loc[1, "MOCA"] = 15  # CN diagnosis, low MoCA at the later visit
    clinical.to_csv(config.inputs.clinical, index=False)
    result = longitudinal(prepare(config).for_group("moca"))
    visits = result.tables["visits"]
    person = visits.loc[visits.subject == "ADNI001S0001"]
    assert set(person.session) == {"M012"}
    assert set(person.group) == {"0-19"}
    assert set(person.diagnosis) == {"CN"}
    assert (person.years_since_first_eligible_scan > .99).all()
    assert (result.tables["binned"].n_subjects <= result.tables["binned"].n_visits).all()


def test_source_free_moca_export_and_plot(prepared, tmp_path):
    import shutil
    from adni_x.plotting import plot
    report = run("moca-reference", prepared.config, render=False)
    portable = tmp_path / "portable"
    shutil.copytree(Path(report["path"]) / "moca-reference", portable)
    manifest_file = portable / "manifest.json"
    manifest = json.loads(manifest_file.read_text())
    manifest["configuration"]["inputs"] = {"tensors": "/unavailable"}
    manifest["configuration"]["cache"] = "/unavailable"
    manifest["sources"] = {}
    manifest_file.write_text(json.dumps(manifest))
    exports = loadmat(portable / "regional_brain_exports.mat")
    assert exports["MOCA_20_28_rtr"].shape == (4, 1)
    np.testing.assert_array_equal(exports["parcel_index"].ravel(), np.arange(4))
    output = plot("moca-reference.regions", portable)[0]
    assert {Path(f).suffix for f in output["files"]} == {".svg", ".png", ".pdf"}
    svg = Path(next(f for f in output["files"] if f.endswith(".svg"))).read_text()
    assert "29-30 z-score" in svg and "CN z-score" not in svg


def test_missing_moca_reference_reports_insufficient_data(config):
    clinical = pd.read_csv(config.inputs.clinical)
    clinical.MOCA = 25
    clinical.to_csv(config.inputs.clinical, index=False)
    report = run("moca-reference", config, render=False)
    assert report["targets"]["moca-reference"]["status"] == "insufficient_data"
    bundle = Path(report["path"]) / "moca-reference"
    assert "29-30" in json.loads((bundle / "manifest.json").read_text())["reason"]
