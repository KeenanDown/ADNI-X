from dataclasses import replace
from pathlib import Path
import json
import numpy as np
import pandas as pd

from adni_x.preparation import prepare
from adni_x.provenance import verify, transaction


def test_no_visit_substitution_for_missing_confounds(config):
    clinical = pd.read_csv(config.inputs.clinical)
    clinical.loc[1, "ICV"] = np.nan
    clinical.to_csv(config.inputs.clinical, index=False)
    prepared = prepare(config)
    selected, subjects = prepared.table("selected"), prepared.table("subjects")
    subject = selected.iloc[0].subject
    assert selected.loc[selected.subject == subject, "month"].item() == 12
    assert subject not in subjects.subject.values
    assert "missing_ICV" in selected.adjustment_exclusion.values
    assert len(prepared.table("visits")) == 48


def test_freesurfer_recovery(config):
    clinical = pd.read_csv(config.inputs.clinical)
    clinical.loc[1, "ICV"] = np.nan
    clinical.to_csv(config.inputs.clinical, index=False)
    root = Path(config.cache).parent / "freesurfer"
    stats = root / "sub-ADNI001S0001" / "stats" / "aseg.stats"
    stats.parent.mkdir(parents=True)
    stats.write_text("# Measure EstimatedTotalIntraCranialVol, eTIV, Estimated Total Intracranial Volume, 1234000.0, mm^3\n")
    p = prepare(replace(config, inputs=replace(config.inputs, freesurfer=str(root))))
    row = p.table("subjects").query("subject == 'ADNI001S0001'").iloc[0]
    assert row.ICV == 1234000 and row.ICV_source == "FreeSurfer_eTIV"


def test_cache_identity_changes_with_inputs_and_settings(config):
    a = prepare(config)
    b = prepare(config)
    assert a.path == b.path and verify(a.path)
    c = prepare(replace(config, statistics=replace(config.statistics, seed=3)))
    assert c.path == a.path  # statistics are not preparation dependencies
    d = prepare(replace(config, cohort=replace(config.cohort, visit="baseline")))
    assert d.path != a.path
    clinical = pd.read_csv(config.inputs.clinical)
    clinical.loc[0, "AGE"] += 1
    clinical.to_csv(config.inputs.clinical, index=False)
    e = prepare(config)
    assert e.path != a.path


def test_interrupted_transaction_is_not_a_completed_cache(tmp_path):
    target = tmp_path / "bundle"
    try:
        with transaction(target) as temp:
            (temp / "half.npy").write_bytes(b"partial")
            raise RuntimeError("simulated interruption")
    except RuntimeError:
        pass
    assert not target.exists()
    assert not list(tmp_path.iterdir())


def test_adjusted_means_and_raw_inputs_preserved(prepared):
    subjects = prepared.table("subjects")
    raw = prepared.array("raw")[subjects.tensor_row.to_numpy(int)]
    adjusted = prepared.array("adjusted")
    np.testing.assert_allclose(adjusted.mean(axis=0), raw.mean(axis=0), atol=1e-9)
    assert prepared.manifest["configuration"]["transform"]["tdmi"] == "none"

