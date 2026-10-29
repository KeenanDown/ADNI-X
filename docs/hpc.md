# HPC pilot and complete rerun

For the version 0.2.0 rerun, start with [the release handoff](hpc-rerun-0.2.0.md), which adds the new output checks and four atlas/estimator submissions to this existing setup guide.

## Readiness boundary

Local implementation and tests use synthetic data. No cluster job has been submitted and no HPC files have been changed. The real-data pilot is **not completed locally**. Do not interpret fixture cohort counts as real ADNI counts.

Copy/clone the package to your normal code directory on the HPC. Install the locked environment using Python 3.12, then make a local configuration from `configs/hpc.example.yaml`. Inputs must be read-only in normal operation. Cache/output directories must be separate writable locations with adequate scratch space.

## QMUL Apocrita setup

QMUL's current documentation uses Slurm (`sbatch`), replacing Grid Engine
(`qsub`); see the [migration guide](https://docs.hpc.qmul.ac.uk/using/uge_to_slurm/)
and [Python guide](https://docs.hpc.qmul.ac.uk/apps/languages/python/).
Confirm the scheduler in your connected terminal with `command -v sbatch`.
If only `qsub` is available, do not submit the Slurm script: first check which
cluster/session you are using.

Copy the source to `~/code/ADNI-X` (for example, using MobaXterm's SFTP pane).
Do not transfer a Windows virtual environment. Keep code and the Linux virtual
environment on persistent storage, with caches in scratch and results in the
project storage paths configured in YAML.

The example uses `/path/to/project/adni-x/results` as a placeholder for results
and validation reports. Replace it with your own writable project directory and
create that directory before submitting jobs. Give a pilot its own `results`
subdirectory and separate cache.

From the HPC terminal:

```sh
cd ~/code/ADNI-X
module avail python
# Example module; check availability in your HPC session.
export ADNI_X_PYTHON_MODULE=python/3.12.1-gcc-12.2.0
module load "$ADNI_X_PYTHON_MODULE"
python --version
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.lock
python -m pip install --no-deps --no-build-isolation -e .
python -m pip check
adni-x --version
cp configs/hpc.example.yaml config.local.yaml
```

Do not continue past a failed install. The lock was verified locally and still
needs Linux verification; do not silently substitute unpinned versions if a
wheel/version is unavailable. Load this same module whenever using the venv,
and export `ADNI_X_PYTHON_MODULE` again after reconnecting.

Edit `config.local.yaml` to supply the verified data and atlas paths. The atlas
example is a required input, not a file shipped with ADNI-X. Run the synthetic
tests in a compute allocation before the real-data pilot. For example, request
an interactive allocation with `srun --pty --mem=8G --time=01:00:00 bash`, then
load the same module, activate the venv and run `python -m pytest`. Exit the
allocation when finished.

Submit from the repository root, with the module variable exported:

```sh
sbatch scripts/qmul_slurm.sh validate --config config.local.yaml
squeue -u "$USER"
# Replace JOBID with the number returned by sbatch.
tail -f adni-x-JOBID.log
sacct -j JOBID --format=JobID,State,ExitCode,Elapsed,MaxRSS
```

The launcher requests one CPU, 32 GB total RAM and four hours as **unmeasured
pilot starting values**; it selects no account or partition. Override resources
with `sbatch --mem=... --time=...` based on the pilot and your allocation. It
activates `.venv` (override with exported `ADNI_X_VENV`), limits numerical library
threads to the allocated CPUs and logs output immediately. More CPUs do not
automatically parallelise the package's analysis targets.

After validation and the input checks below, create the separate pilot described
in the next section and submit:

```sh
sbatch scripts/qmul_slurm.sh run all --config config.pilot.yaml --include-sensitivity
```

Once the pilot passes, submit the full workflow using measured resource limits:

```sh
sbatch scripts/qmul_slurm.sh pipeline "$PWD/config.local.yaml"
```

Do not submit simultaneous preparation jobs against the same cache. The
launcher also accepts individual commands such as `prepare --config ...`,
`run group-matrices --config ...`, or `plot all --results ...`.

## Input checks

1. Check the S8b estimator, atom order, source/target axes, and complete cortical/subcortical parcel order against the upstream implementation. Supply the full atlas CSV.
2. Confirm the original time-series variable and matching filename entities. The defaults require 150–250 points and TR 2.875–3.150 seconds, both inclusive. Each tensor JSON must link via `source_bold` to an accessible BOLD JSON containing `RepetitionTime` in seconds. Review `TR_seconds`, `tr_status` and exclusion reasons; metadata failures are exclusions, not passes.
3. Run `adni-x validate --config config.local.yaml`. Inspect `inventory.csv`, reasons, partial-identity matches, unknown diagnoses, missing/invalid MoCA, conflicting clinical records and both grouping counts. Resolve ambiguous files through the optional CSV rather than editing raw data.
4. Confirm available disk space. 2,907 × 232 × 232 × 16 float64 values are about 20 GB before the two grouping-specific adjusted stores, temporary work and result files. The raw store is shared by both groupings. Blocks limit computation memory, not on-disk cache size. Start with one process; avoid multiplying independent full preparation jobs.

## Small real-data pilot

Use a separate pilot input directory containing symlinks to a deterministic, documented selection of approximately 10–20 participants per diagnosis and all their candidate visits. Include boundary lengths, an advanced-length scan, and any known ambiguous runs. Never copy a misleading latest-only subset: filtering must be exercised before selecting the latest eligible scan. Supply matching series symlinks and the same clinical/atlas sources; point a separate pilot YAML at these directories and separate cache/results paths.

Tensor symlinks use the JSON beside the resolved tensor target. If copying tensors instead, include their JSON sidecars and preserve accessible `source_bold` links. Include the known 197-point/6.015-second scan in the pilot and confirm `above_tr_window`, then check whether an earlier eligible visit is selected. TR filtering also removes ineligible visits from longitudinal outputs.

Keep default methods and use `resamples: 100`, `permutations: 100` only for this pilot. These small counts are for runtime checking, not inference. Run:

```sh
adni-x validate --config config.pilot.yaml
adni-x run all --config config.pilot.yaml --include-sensitivity
```

Check `run.json` for failed and insufficient-data targets. Check source identities, selected visits, group counts, ICV sources, adjustment exclusions, feature validity, design rank and figure labels. Inspect several tensors and reductions against independently calculated values. Expected changes from the old notebook include the length rule, filter-before-latest selection, complete cases, visit age, removed diagonals, mean-preserving adjustment and corrected inference. Do not require corrected p-values to match old displayed output.

Move a copied result bundle to a directory with no access to its original inputs and run `adni-x plot all --results /path/to/bundle`. Confirm SVG/PNG rendering and the parcel-indexed MATLAB exports.

## Full run

For existing configurations, the new TR defaults apply even when the YAML omits those fields. Update/reinstall the package, run validation and rerun analyses with the existing configuration. New source/code/configuration fingerprints produce separate preparation and result directories; existing numerical bundles are not edited. Replotting an old bundle alone does **not** apply the new exclusion.

Restore the full inputs and 10,000 resamples/permutations. Use the scheduler-neutral wrapper inside the HPC's approved compute-job mechanism:

```sh
sh scripts/run_hpc.sh /absolute/path/config.local.yaml
```

The wrapper does not select a scheduler, queue or account. It sets numerical-library thread counts before starting Python, logs validate/prepare/run phases, and returns a nonzero exit code on failure. Choose HPC resource limits based on the pilot's measured peak memory/time. Do not run the complete analysis on a login node.

For independent jobs after preparation, invoke `adni-x run ANALYSIS --config … --no-plots`; use separate job logs. Raw scans are never modified. Plotting can run later without inputs.

## Resume and diagnose

Repeat the same command to reuse valid preparation and completed numerical targets. Failed targets are recorded under `failures/` and retried; other targets continue. Insufficient-data targets are completed reports with explicit reasons, not fabricated figures. Changes to source contents, method settings or package source create a new cache/result identity.

Incomplete temporary directories are never treated as completed caches. A checksum failure is a hard error: move the affected bundle aside for inspection, then rerun. Do not edit a sealed result in place. Preserve the run manifest, dependency lock and code revision with the final analysis archive.

## Diagnosis and MoCA update

Existing configuration files automatically run both groupings with `adni-x run all --config ...`: 28 default targets, or 32 with `--include-sensitivity` (both TDMI-ratio targets and both MMSE taxonomy targets). Named MMSE targets remain optional and use the diagnosis fitting cohort. Reinstall the updated checkout and start a new run; replotting old results cannot add these analyses. The implementation fingerprint creates new cache/result directories without overwriting earlier runs. Do not update an editable installation while existing jobs are still using it.

Check both `selected_groups` and `selected_moca_categories` in validation. Preparation's `audit.json` includes separate `groupings.diagnosis` and `groupings.moca` cohort counts and adjustment exclusions. Each result's manifest identifies its grouping and reference; `cohort_adjusted.csv` retains diagnostic `group`, `MOCA_category`, and `analysis_group`. Inspect missing-score exclusions and any `insufficient_data` targets before interpreting MoCA output. The highest MoCA band can be small; a missing reference is reported, never replaced with a diagnostic group.
