from dataclasses import replace
from pathlib import Path
import json
import numpy as np
import pandas as pd
import pytest

from adni_x.catalogue import CATALOGUE
from adni_x.config import Config, Transform
from adni_x.numerics import Design, fdr
from adni_x.preparation import prepare
from adni_x.provenance import verify
from adni_x.runner import run


def test_all_missing_confounds_leave_raw_targets_usable(config):
    clinical = pd.read_csv(config.inputs.clinical)
    clinical.ICV = np.nan
    clinical.to_csv(config.inputs.clinical, index=False)
    dataset = prepare(config)
    assert len(dataset.table("subjects")) == 0
    report = run("all", config, render=False)
    assert report["targets"]["tdmi"]["status"] == "complete"
    assert report["targets"]["longitudinal"]["status"] == "complete"
    assert report["targets"]["group-matrices"]["status"] == "insufficient_data"
    assert not any(v["status"] == "failed" for v in report["targets"].values())


def test_numeric_imputation_uses_each_training_column(config):
    clinical = pd.read_csv(config.inputs.clinical)
    clinical.loc[1, ["ICV", "PTEDUCAT", "PTGENDER"]] = [np.nan, np.nan, None]
    clinical.to_csv(config.inputs.clinical, index=False)
    dataset = prepare(replace(config, transform=replace(config.transform, missing="median")))
    design = json.loads((dataset.path / "design.json").read_text())
    assert len(dataset.table("subjects")) == 24
    assert design["medians"]["ICV"] > 1e6
    assert design["medians"]["education"] < 30
    assert {i["confound"] for i in design["imputations"]} == {"ICV", "education", "sex"}


def test_tdmi_mode_and_no_confounds_have_separate_cache(config):
    dataset = prepare(config)
    normalised_config = replace(config, transform=replace(config.transform, tdmi="tdmi", confounds=()))
    normalised = prepare(normalised_config)
    assert dataset.path != normalised.path
    values = normalised.array("adjusted")
    mask = ~np.eye(4, dtype=bool)
    np.testing.assert_allclose(values[:, mask, :].sum(axis=-1), 1)


def test_corrupt_cache_and_failed_target_do_not_hide_failures(config, monkeypatch):
    import adni_x.runner as runner
    from adni_x.catalogue import Analysis
    def broken(_):
        raise RuntimeError("simulated target failure")
    monkeypatch.setitem(runner.ANALYSES, "broken", Analysis("broken", "test", "", broken, ()))
    report = run("all", config, render=False)
    assert report["targets"]["broken"]["status"] == "failed"
    assert report["targets"]["tdmi"]["status"] == "complete"
    assert (Path(report["path"]) / "failures" / "broken.json").exists()
    dataset = prepare(config)
    with (dataset.path / "raw_tdmi.npy").open("ab") as f:
        f.write(b"corruption")
    assert not verify(dataset.path)
    with pytest.raises(ValueError, match="corrupt"):
        prepare(config)


def test_config_validation_and_relative_paths(tmp_path):
    config_file = tmp_path / "config.yaml"
    config_file.write_text('inputs:\n  tensors: scans\n  series: series\n  atlas: atlas.csv\n  clinical: clinical.csv\n')
    c = Config.load(config_file)
    assert c.inputs.tensors == str(tmp_path / "scans")
    with pytest.raises(ValueError, match="Unknown configuration"):
        config_file.write_text('typo: true\n')
        Config.load(config_file)
    with pytest.raises(ValueError):
        Config(transform=Transform(confounds=("wrong",)))


def test_fdr_inclusive_boundary():
    q, reject = fdr([.05])
    assert q.item() == .05 and reject.item()


def test_source_notebook_coverage_is_exhaustive():
    root = Path(__file__).resolve().parents[1]
    data = json.loads((root / "docs" / "notebook_cells.json").read_text())
    assert [r["cell"] for r in data["cells"]] == list(range(297))
    assert all(r["targets"] or r["category"] in {"preparation", "unimplemented TODO (outside v1)", "empty cell"}
               for r in data["cells"])


def test_concurrent_independent_status_updates_are_not_lost(tmp_path):
    from concurrent.futures import ThreadPoolExecutor
    from adni_x.runner import publish_status
    def update(i):
        name = f"analysis-{i}"
        publish_status(tmp_path, {"run_id": "one", "targets": {name: {"status": "complete"}}}, name)
    with ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(update, range(20)))
    report = json.loads((tmp_path / "run.json").read_text())
    assert set(report["targets"]) == {f"analysis-{i}" for i in range(20)}
