from dataclasses import replace
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
import h5py
from scipy.io import savemat, loadmat

from adni_x.inputs import inventory, read_tensor, series_metadata, entities, read_atlas


def test_boundaries_and_filter_before_latest(config):
    files = sorted(Path(config.inputs.series).glob("*ses-M012*"))
    for path, count in zip(files[:5], [149, 150, 250, 251, 1000]):
        savemat(path, {"activity": np.ones((4, count))})
    audit = inventory(config).scans
    for path, expected_month, reason in zip(files[:5], [0, 12, 12, 0, 0],
                                           ["below_length_window", "eligible", "eligible", "above_length_window", "above_length_window"]):
        subject = entities(path)["sub"]
        rows = audit.loc[audit.subject == subject]
        assert rows.loc[rows.selected, "month"].item() == expected_month
        assert rows.loc[rows.month == 12, "reason"].item() == reason


def test_missing_and_ambiguous_series(config):
    paths = sorted(Path(config.inputs.series).glob("*ses-M012*"))
    paths[0].unlink()
    original = loadmat(paths[1])["activity"]
    savemat(paths[1].with_name(paths[1].name.replace("task-rest", "task-rest_run-2")), {"activity": original})
    audit = inventory(config).scans
    assert "missing_series" in audit.reason.values
    assert "ambiguous_series" in audit.reason.values


def test_conflicting_clinical_and_resolution(config):
    clinical = pd.read_csv(config.inputs.clinical)
    duplicate = clinical.iloc[1].copy()
    duplicate.DX = "Dementia"
    pd.concat([clinical, duplicate.to_frame().T], ignore_index=True).to_csv(config.inputs.clinical, index=False)
    inv = inventory(config)
    bad = inv.scans.loc[inv.scans.reason == "conflicting_clinical_records"]
    assert len(bad) == 1
    resolution = Path(config.inputs.clinical).parent / "resolutions.csv"
    pd.DataFrame({"tensor_file": bad.tensor_file, "clinical_row": [1]}).to_csv(resolution, index=False)
    repaired = inventory(replace(config, inputs=replace(config.inputs, resolutions=str(resolution))))
    assert "conflicting_clinical_records" not in repaired.scans.reason.values


def test_duplicate_visit_requires_explicit_choice(config):
    path = sorted(Path(config.inputs.tensors).glob("*ses-M012*"))[0]
    duplicate = path.with_name(path.name.replace("task-rest", "task-rest_run-2"))
    savemat(duplicate, {"PhiIDTensor": loadmat(path)["PhiIDTensor"]})
    inv = inventory(config)
    assert (inv.scans.reason == "ambiguous_same_visit_scans").sum() == 2
    resolution = Path(config.inputs.clinical).parent / "resolutions.csv"
    pd.DataFrame({"tensor_file": [str(duplicate)], "use_scan": [False]}).to_csv(resolution, index=False)
    fixed = inventory(replace(config, inputs=replace(config.inputs, resolutions=str(resolution))))
    assert (fixed.scans.reason == "ambiguous_same_visit_scans").sum() == 0


def test_matlab_storage_and_explicit_axes(config, tmp_path):
    a = np.arange(4*4*16, dtype=float).reshape(4, 4, 16)
    matlab5, matlab73, transposed = tmp_path / "v5.mat", tmp_path / "v73.mat", tmp_path / "transposed.mat"
    savemat(matlab5, {"PhiIDTensor": a})
    with h5py.File(matlab73, "w") as f:
        ds = f.create_dataset("PhiIDTensor", data=a.transpose(2, 1, 0))
        ds.attrs["MATLAB_class"] = np.bytes_("double")
    savemat(transposed, {"PhiIDTensor": a.transpose(2, 0, 1)})
    np.testing.assert_array_equal(read_tensor(matlab5, config), a)
    np.testing.assert_array_equal(read_tensor(matlab73, config), a)
    with pytest.raises(ValueError, match="Expected"):
        read_tensor(transposed, config)
    explicit = replace(config, inputs=replace(config.inputs, tensor_axes=(1, 2, 0)))
    np.testing.assert_array_equal(read_tensor(transposed, explicit), a)


def test_series_ambiguity_and_atlas_order(config, tmp_path):
    path = tmp_path / "ambiguous.mat"
    savemat(path, {"a": np.ones((4, 200)), "b": np.ones((4, 197))})
    with pytest.raises(ValueError, match="unambiguous"):
        series_metadata(path, 4)
    assert series_metadata(path, 4, "b")[1] == 197
    atlas = pd.read_csv(config.inputs.atlas).iloc[::-1]
    atlas.to_csv(config.inputs.atlas, index=False)
    with pytest.raises(ValueError, match="tensor order"):
        read_atlas(config.inputs.atlas, 4)


def test_estimator_discovery_and_entity_conflicts(config, tmp_path):
    ccs = replace(config, cohort=replace(config.cohort, estimator="ccs"))
    with pytest.raises(ValueError, match="No tensor"):
        inventory(ccs)
    path = sorted(Path(config.inputs.tensors).glob("*.mat"))[0]
    savemat(path.with_name(path.name.replace("desc-phiid", "desc-CCS-phiid")), {"PhiIDTensor": loadmat(path)["PhiIDTensor"]})
    assert len(inventory(ccs).scans) == 1
    with pytest.raises(ValueError, match="conflicting"):
        entities(tmp_path / "sub-ADNI001S9999" / path.name)

