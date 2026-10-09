#!/usr/bin/env python3
"""Synthetic integration check at both real atlas sizes; no ADNI inputs used.

Run from the source root after installation:
    python scripts/verify_rerun_release.py --output /path/to/verification
Large synthetic tensors/caches live in a temporary directory and are removed
on completion. Only aggregate checks and selected editable figures are kept.
"""
import argparse
import ast
import json
from pathlib import Path
import shutil
import tempfile

import numpy as np
import pandas as pd
from scipy.io import savemat

from adni_x import __version__
from adni_x.catalogue import CATALOGUE, FIGURES
from adni_x.config import Config, Inputs, Cohort, Statistics
from adni_x.constants import ATOMS
from adni_x.plotting import plot
from adni_x.provenance import environment, verify
from adni_x.runner import run


def synthetic_inputs(root, parcels):
    root.mkdir()
    tensors, series = root/"tensors", root/"series"
    tensors.mkdir(); series.mkdir()
    network = np.arange(parcels)*8//parcels+1
    pd.DataFrame({"parcel_index": np.arange(parcels), "parcel_id": [f"SYNTH{i}" for i in range(parcels)],
                  "name": [f"Synthetic parcel {i}" for i in range(parcels)], "network_id": network}).to_csv(root/"atlas.csv", index=False)
    rng = np.random.default_rng(20261009+parcels)
    swap = {"r": "r", "x": "y", "y": "x", "s": "s"}
    permutation = [ATOMS.index(swap[a[0]]+"t"+swap[a[2]]) for a in ATOMS]
    records = []
    for gi, group in enumerate(("CN", "MCI", "Dementia")):
        for j in range(6):
            subject = gi*6+j+1
            age, education, icv = rng.uniform(60, 85), rng.uniform(10, 20), rng.uniform(1.2e6, 1.6e6)
            male = j % 2
            prefix = f"sub-ADNI001S{subject:04d}_ses-M000_task-rest_atlas-ExtSchaefer{parcels}"
            for estimator in ("mmi", "ccs"):
                tensor = rng.normal(.5, .025, size=(parcels, parcels, 16))
                tensor += .005*(age-70)+.004*(education-15)+3e-7*(icv-1.4e6)+.04*male
                tensor[:, :, 0] += gi*.04*(network[:, None]+network[None, :])/16
                tensor[:, :, 15] += gi*.06*(network[:, None]+network[None, :])/16
                tensor = (tensor+tensor.transpose(1, 0, 2)[..., permutation])/2
                tensor[np.arange(parcels), np.arange(parcels)] = 0
                suffix = "-CCS" if estimator == "ccs" else ""
                savemat(tensors/f"{prefix}_desc{suffix}-phiidtensor.mat", {"PhiIDTensor": tensor})
            savemat(series/f"{prefix}_ROI_activity.mat", {"activity": rng.normal(size=(parcels, 200))})
            records.append({"PTID": f"001_S_{subject:04d}", "VISCODE": "bl", "DX": group,
                            "PTGENDER": "Male" if male else "Female", "PTEDUCAT": education,
                            "AGE": age, "Years_bl": 0, "ICV": icv, "MMSE": [29, 25, 15][gi],
                            "MOCA": [30, 25, 15][gi], "EXAMDATE": "2020-01-01"})
    pd.DataFrame(records).to_csv(root/"clinical.csv", index=False)
    return Inputs(tensors=str(tensors), series=str(series), clinical=str(root/"clinical.csv"), atlas=str(root/"atlas.csv"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    notebook = json.loads((Path(__file__).resolve().parents[1]/"notebooks/ADNI_X.ipynb").read_text())
    for cell in notebook["cells"]:
        if cell["cell_type"] == "code":
            ast.parse("".join(cell["source"]))
            assert cell["outputs"] == []
    assert len({c["id"] for c in notebook["cells"]}) == len(notebook["cells"])
    summary = {"synthetic_only": True, "source_package_version": __version__, "environment": environment(), "analyses": len(CATALOGUE),
               "registered_figures": len(FIGURES), "notebook_cells": len(notebook["cells"]), "configurations": {}}
    with tempfile.TemporaryDirectory(prefix="adnix-release-") as temporary:
        root = Path(temporary)
        for parcels in (232, 116):
            inputs = synthetic_inputs(root/f"inputs-{parcels}", parcels)
            for estimator in ("mmi", "ccs"):
                name = f"{parcels}-{estimator}"
                config = Config(inputs=inputs, cohort=Cohort(parcels=parcels, estimator=estimator),
                                statistics=Statistics(resamples=5, permutations=5),
                                cache=str(root/"cache"), results=str(root/"results"))
                print(f"Verifying {name}", flush=True)
                report = run("all", config, render=False)
                assert all(t["status"] == "complete" for t in report["targets"].values()), report["targets"]
                assert not any(t.endswith("mmse") for t in report["targets"])
                expected = parcels*(parcels-1)//2
                checks = {"participants": 18, "networks": 8, "unique_pairs": expected,
                          "completed_targets": len(report["targets"]), "figures": []}
                for target in ("group-matrices", "group-matrices-moca"):
                    bundle = Path(report["targets"][target]["path"])
                    assert verify(bundle)
                    edges = pd.read_csv(bundle/"edges.csv")
                    assert (edges.parcel_i < edges.parcel_j).all()
                    assert all(len(sub) == expected for _, sub in edges.groupby(["atom", "group_a", "group_b"]))
                    checks[target+"_symmetry"] = json.loads((bundle/"details.json").read_text())["symmetry"]
                inference = pd.read_csv(Path(report["targets"]["coefficients"]["path"])/"inference.csv")
                assert len(inference) == 2*expected*4
                assert set(inference.status) == {"ok"}
                checks["coefficient_tests"] = len(inference)
                selected = ["taxa-principal-dx.panels", "taxa-principal-moca.panels",
                            "network-effects.MCI-vs-CN.cohens_d.all-atoms", "network-effects.MCI-vs-CN.cohens_d.rtr-sts",
                            "hierarchies.sts_deviation", "hierarchies-moca.sts_deviation",
                            "group-matrices.sts", "group-matrices-moca.rtr",
                            "coefficients.rtr.significant", "coefficients.combined.significant"]
                # Inspect representative exports from both atlas sizes and the
                # CCS route, without keeping the large temporary datasets.
                if name in ("232-mmi", "116-ccs"):
                    destination = output/name
                    destination.mkdir(exist_ok=True)
                    for figure in selected:
                        target = "coefficients" if figure.startswith("coefficients.") else FIGURES[figure]
                        exported = plot(figure, report["targets"][target]["path"])[0]
                        for file in exported["files"]:
                            suffix = Path(file).suffix
                            target_file = destination/f"{figure}{suffix}"
                            shutil.copy2(file, target_file)
                            checks["figures"].append(str(target_file.relative_to(output)))
                        sidecar = Path(exported["files"][0]).with_suffix(".json")
                        record = json.loads(sidecar.read_text())
                        record["original_export_files"] = record["files"]
                        record["files"] = {f"{figure}{Path(filename).suffix}": checksum
                                           for filename, checksum in record["files"].items()}
                        (destination/f"{figure}.json").write_text(json.dumps(record, indent=2)+"\n")
                summary["configurations"][name] = checks
                (output/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps({"report": str(output/"summary.json"), "verified": list(summary["configurations"])}), flush=True)


if __name__ == "__main__":
    main()
