from dataclasses import replace
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from scipy.io import savemat

from adni_x.config import Config, Inputs, Cohort, Transform, Statistics
from adni_x.constants import ATOMS


def make_dataset(root, per_group=8):
    root = Path(root)
    tensor_dir, series_dir = root / "tensors", root / "series"
    tensor_dir.mkdir(parents=True); series_dir.mkdir()
    rng = np.random.default_rng(2026)
    p = 4
    pd.DataFrame({"parcel_index": range(p), "parcel_id": [f"ROI{i+1}" for i in range(p)],
                  "name": [f"Parcel {i+1}" for i in range(p)], "network_id": [1, 1, 2, 2]}).to_csv(root / "atlas.csv", index=False)
    records = []
    for gi, group in enumerate(["CN", "MCI", "Dementia"]):
        for j in range(per_group):
            n = gi*per_group + j + 1
            subject = f"ADNI001S{n:04d}"
            age, icv = 60 + rng.uniform(0, 20), 1.3e6 + rng.normal(0, 2e4)
            education = 12 + j % 7
            for month in (0, 12):
                name = f"sub-{subject}_ses-M{month:03d}_task-rest_atlas-ExtSchaefer{p}"
                tensor = rng.gamma(3, .05, (p, p, 16))
                tensor += .0001*age + .0002*education + 1e-8*icv
                tensor[:, :, 0] += gi*.03*np.arange(1, p + 1)[:, None]
                tensor[:, :, 15] -= gi*.015*np.arange(1, p + 1)[:, None]
                # Physical ROI exchange swaps unique-X/unique-Y atoms; rtr/sts
                # and total TDMI are symmetric, while directed atoms may differ.
                swap = {"r": "r", "x": "y", "y": "x", "s": "s"}
                permutation = [ATOMS.index(swap[a[0]]+"t"+swap[a[2]]) for a in ATOMS]
                tensor = (tensor+tensor.transpose(1, 0, 2)[..., permutation])/2
                tensor[np.arange(p), np.arange(p)] = 100  # detect accidental self-edge inclusion
                savemat(tensor_dir / f"{name}_desc-phiidtensor.mat", {"PhiIDTensor": tensor})
                savemat(series_dir / f"{name}_ROI_activity.mat", {"activity": rng.normal(size=(p, 200))})
                records.append({"PTID": f"001_S_{n:04d}", "VISCODE": "bl" if not month else "m12",
                                "DX": group, "PTGENDER": "Female" if j % 2 else "Male", "PTEDUCAT": education,
                                "AGE": age, "Years_bl": month/12, "ICV": icv, "MMSE": [29, 25, 15][gi],
                                "MOCA": [30, 25, 15][gi], "EXAMDATE": f"{2020+month//12}-01-{j+1:02d}"})
    pd.DataFrame(records).to_csv(root / "clinical.csv", index=False)
    return Config(inputs=Inputs(tensors=str(tensor_dir), series=str(series_dir), clinical=str(root / "clinical.csv"), atlas=str(root / "atlas.csv")),
                  cohort=Cohort(parcels=p), transform=Transform(block_features=37), statistics=Statistics(resamples=20, permutations=20),
                  cache=str(root / "cache"), results=str(root / "results"))


@pytest.fixture
def config(tmp_path):
    return make_dataset(tmp_path)


@pytest.fixture(scope="session")
def prepared(tmp_path_factory):
    from adni_x.preparation import prepare
    return prepare(make_dataset(tmp_path_factory.mktemp("prepared")))


@pytest.fixture(scope="session")
def bundles(prepared):
    from adni_x.runner import run
    report = run("all", prepared.config, render=False, include_sensitivity=True)
    assert all(item["status"] == "complete" for item in report["targets"].values()), report["targets"]
    return Path(report["path"])

