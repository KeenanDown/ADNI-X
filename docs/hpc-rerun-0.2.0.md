# ADNI-X 0.2.0: HPC rerun handoff

This release prepares the package for a new cluster run. Local checks use synthetic data; no real-data cluster jobs have been submitted. Keep source tensors read-only and place caches/results in the established scratch/project locations. Read the existing [HPC setup guide](hpc.md) for module, Slurm and input-path details.

## Install the release

Commit and push the source changes, then `git pull` in the HPC checkout before reinstalling. The repository includes the notebook, tests, locked requirements, configurations and Slurm scripts. Keep an existing running job's checkout untouched until it finishes. Release archives are optional alternatives: unpack `dist/adni_x-0.2.0.tar.gz` into a new versioned code directory, or install the wheel. Do not copy the Windows/WSL virtual environment.

Load the same Python 3.12 module used for the existing cluster virtual environment, then install from the unpacked source root:

```sh
source /absolute/path/to/cluster/venv/bin/activate
python -m pip install -r requirements.lock
python -m pip install --no-deps --no-build-isolation .
python -m pip check
adni-x --version
python -m pytest
```

The version command must report `ADNI-X 0.2.0`. To use a wheel, replace the source-install command with `python -m pip install --no-deps --force-reinstall /absolute/path/adni_x-0.2.0-py3-none-any.whl`. Preserve the previous code version and results for comparison. New source hashes create distinct preparation caches and result IDs.

The package adds no runtime dependency. Plotting uses the Agg backend; GUI/Jupyter sessions are unnecessary. Symmetry and HC3 calculations use bounded feature blocks over memory-mapped data. Keep the Slurm script's single-process/BLAS-thread settings and inspect pilot memory/time before choosing production resources. Disk space must accommodate fresh raw and adjusted caches for each configuration.

## Create four configurations

For the Git workflow, review the four committed `config.232.mmi.local.yaml`, `config.232.ccs.local.yaml`, `config.116.mmi.local.yaml` and `config.116.ccs.local.yaml` files and use their existing paths. Alternatively, copy `configs/hpc.example.yaml` into four verified local files:

| File | `cohort.parcels` | `cohort.estimator` | Atlas mapping |
|---|---:|---|---|
| `config.232-mmi.yaml` | 232 | `mmi` | Full 232-region CSV |
| `config.232-ccs.yaml` | 232 | `ccs` | Same 232-region CSV |
| `config.116-mmi.yaml` | 116 | `mmi` | Full 116-region CSV |
| `config.116-ccs.yaml` | 116 | `ccs` | Same 116-region CSV |

Verify `inputs.atlas` and all dataset paths separately. Configured custom tensor globs must match the selected estimator/atlas; remove stale custom globs to use the package's defaults. Preserve the inclusive 150–250 retained-time-point and 2.875–3.150-second TR filters before visit selection, plus the established covariates and transformation settings. A common writable cache/results root is supported because configuration and source hashes separate runs. Give the pilot separate roots.

## Pilot and submissions

Use a small audited input subset with all three diagnostic and MoCA groups for the pilot. Do not alter the full dataset to create it. From the new source root, after exporting the established module and virtual-environment settings:

```sh
export ADNI_X_PYTHON_MODULE=YOUR_VERIFIED_PYTHON_MODULE
export ADNI_X_VENV=/absolute/path/to/cluster/venv
sbatch scripts/qmul_slurm.sh pipeline /absolute/path/config.pilot.yaml
```

Inspect the pilot's audit, `run.json`, symmetry discrepancies, coefficient statuses and figures. Production submissions use the same script, one configuration per job:

```sh
sbatch scripts/qmul_slurm.sh pipeline /absolute/path/config.232-mmi.yaml
sbatch scripts/qmul_slurm.sh pipeline /absolute/path/config.232-ccs.yaml
sbatch scripts/qmul_slurm.sh pipeline /absolute/path/config.116-mmi.yaml
sbatch scripts/qmul_slurm.sh pipeline /absolute/path/config.116-ccs.yaml
```

Run a first production configuration before submitting the others if scratch, memory or runtime remain uncertain. The default pipeline validates, prepares and runs all default analyses with headless figures. Diagnosis and MoCA are included; MMSE is omitted. To run optional MMSE explicitly:

```sh
sbatch scripts/qmul_slurm.sh run taxa-six-mmse --config /absolute/path/config.232-mmi.yaml
```

`run all` now selects 28 diagnosis/MoCA targets. `--include-sensitivity` requests both grouping-specific longitudinal TDMI-ratio branches and both optional MMSE taxonomy targets, so use named MMSE targets when the ratio branches are not wanted.

## Check the outputs

- `taxa-principal-dx` and `taxa-principal-moca` contain upward, downward, transfer and storage panels; original six-mode diagnosis/MoCA bundles are also retained.
- `network-effects` includes `MCI-vs-CN.cohens_d.all-atoms` and `.rtr-sts`. Positive values mean MCI minus CN. Both subsets reuse correction across all 16 atoms × networks.
- `group-matrices` and `group-matrices-moca` export one row per unique pair in `edges.csv`: 26,796 per atom/contrast for 232 parcels, 6,670 for 116. Each of the two atoms has three contrasts. Significant full-matrix plots mirror these results.
- `coefficients/inference.csv` and `coefficients-moca/inference.csv` contain pooled non-intercept rtr/sts covariate tests for their separately fitted cohorts, with HC3 standard errors, approximate t tests, 95% CIs and corrected q-values. The original unthresholded maps remain; new `.significant` maps show corrected results. Inspect unavailable statuses instead of treating missing p-values as nonsignificance.
- `hierarchies.*_deviation` subtracts each parcel's CN mean, so CN is zero. MoCA counterparts subtract 29–30. The endpoint-minus-reference ordering is descriptive. Axes identify the baseline.
- MoCA uses 29–30, 20–28 and 0–19; missing scores affect only cognitive outputs. Check real cohort counts and available group sizes.
- A rtr/sts symmetry failure is an input/convention diagnostic, not a reason to bypass validation or silently retain only one direction.

`adni-x plot` can refresh saved profiles and other existing figures without input tensors. New coefficient inference and unique-pair tables require a fresh analysis run; refreshing old figures cannot create those results. Historical result files remain unchanged.
