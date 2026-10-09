# ADNI-X

Reproducible ADNI analyses starting from **completed S8b PhiID tensors**. ADNI-X does not recompute PhiID or preprocess MRI scans.

The workflow is **prepare once → run independent analyses → replot saved numerical results**. A Jupyter kernel, MATLAB, and the original notebook's intermediate pickles are not required.

Version **0.2.0** adds separate upward/downward, transfer and storage panels, default MoCA analyses, unique symmetric rtr/sts pair tests, CN-centred profiles and HC3 coefficient inference. See the [0.2.0 HPC rerun handoff](docs/hpc-rerun-0.2.0.md) for installation, pilot checks and the four atlas/estimator submissions.

## Installation

Use Python 3.12 (Python ≥3.11 is supported by the package). From this repository:

```sh
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.lock
python -m pip install --no-deps --no-build-isolation -e .
```

The lock records the exact versions used for local verification; the same versions must also be installed and tested on the HPC. No job submission or remote changes are performed by installation.

## First run

Copy `configs/hpc.example.yaml` to `config.local.yaml`, verify its paths, and supply a **full, ordered atlas CSV**. See [input contracts](docs/inputs.md) before the first run. The example atlas path is intentionally a required user-supplied mapping; subcortical parcel identities are not guessed.

```sh
adni-x validate --config config.local.yaml
adni-x prepare --config config.local.yaml
adni-x list
adni-x run group-matrices --config config.local.yaml
adni-x run all --config config.local.yaml
```

`run` produces numerical bundles and figures. Use `--no-plots` to compute only. MMSE taxonomy targets and `longitudinal-tdmi-ratio` are omitted by `all` unless `--include-sensitivity` is supplied; requesting a target by name is also explicit opt-in. Diagnosis and MoCA outputs are included by default, including exploratory ones.

To regenerate just one figure, supply its saved analysis directory or the parent run directory:

```sh
adni-x plot group-matrices.rtr --results /path/to/results/RUN_ID/group-matrices
adni-x plot all --results /path/to/results/RUN_ID
```

Replotting needs **only the result bundle**, not raw data or preparation caches. Numerical outputs have stable names inside immutable bundles. New editable SVG, PDF and 300-dpi PNG files have informative names and their actual UTC creation timestamp. Every plot has a JSON provenance sidecar. See the [manuscript plotting guide](docs/plotting.md) for the original-notebook styling and how to refresh an existing run.

## Defaults

- Every analysis family runs for **diagnosis (CN/MCI/AD)** and **MoCA category (29–30 / 20–28 / 0–19)** under `run all`: 28 default targets, plus two opt-in TDMI-ratio and two opt-in MMSE taxonomy targets. MMSE categories are calculated only for explicitly requested MMSE analyses.
- 232 parcels, MMI, no TDMI division.
- Inclusive **150–250 retained time points** and **2.875–3.150-second repetition time**, then latest eligible scan per participant. TR comes from the BOLD JSON linked by the tensor sidecar; missing or invalid acquisition metadata excludes the scan. These checks do not establish a Basic/Advanced protocol classification.
- Confounds: sex, education, visit age (`AGE + Years_bl`) and ICV. Missing ICV may be recovered from FreeSurfer; otherwise use complete cases after visit selection, with no older-visit substitution.
- Adjustment preserves the fitting cohort's mean. Self-connections are excluded; directed regional/network reductions retain both directions and regional sums use axis 1. Symmetric rtr/sts matrix tests use each unique pair once after validation.
- Welch tests, two-sided α=.05, BH-FDR, seed 42, 10,000 resamples/permutations.
- Prediction uses a stratified 75/25 participant split with training-only preprocessing.

Each grouping has its own cohort and confound fit. Missing MoCA excludes a participant only from MoCA analyses; missing diagnosis excludes them only from diagnostic analyses. Missing grouping scores do not trigger selection of an older visit. Raw TDMI requires the relevant group label but not complete confounds; longitudinal plots use all eligible visits with that visit's label. Cohort exports preserve diagnosis in `group` and the active label in `analysis_group`. Observation tables use `group` for the active label and retain `diagnosis` and `MOCA_category` alongside it.

Diagnostic target IDs remain unchanged. MoCA targets generally add `-moca`, for example `network-effects-moca`; the reference/PCA target is `moca-reference`, and taxa targets include `taxa-four-moca`, `taxa-six-moca` and `taxa-principal-moca`. A named target runs only that grouping; `run all` runs both. MoCA reference analyses use 29–30 as their reference, and MoCA prediction compares 29–30 against 0–19 with 20–28 projection-only. These score bands are not clinical diagnoses. `taxa-four-mmse` and `taxa-six-mmse` remain optional sensitivity targets and use the diagnosis fitting cohort.

## Python and notebook use

```python
from adni_x import prepare, run, plot, validate
from adni_x.config import Config

config = Config.load("config.local.yaml")
audit = validate(config)
prepared = prepare(config)
report = run("all", config)  # diagnosis and MoCA; numerical results and figures
figures = plot("hierarchies.hierarchy", report["path"])
```

`notebooks/ADNI_X.ipynb` is a thin entry point to this API. Numerical analysis functions are also available in `adni_x.analyses`, `adni_x.prediction`, and `adni_x.longitudinal`; they accept a prepared dataset (`prepared.for_group("moca")` for MoCA) and return a `Result` containing tables, numeric arrays, figure recipes and metadata.

## Documentation and verification

- [Analysis catalogue and notebook coverage](docs/analysis_catalogue.md)
- [Methods, correction families and legacy options](docs/methods.md)
- [Input and atlas contracts](docs/inputs.md)
- [Bug ledger](docs/bug_ledger.md)
- [HPC pilot and rerun guide](docs/hpc.md)
- [QMUL Apocrita installation and Slurm submission](docs/hpc.md#qmul-apocrita-setup)
- [Local verification results and remaining HPC checks](docs/verification.md)

```sh
python -m pytest
python -m build --no-isolation
```

Tests use synthetic participants only. Real HPC tensor conventions, atlas identity and cohort counts must be checked in the pilot; passing local tests does not certify the real dataset. Dataset files, caches and participant-level results are ignored by Git. Do not commit them.

Updating to this grouping workflow creates a new cache/run identity. Existing YAML files need no extra settings, but saved results must be recomputed to add the MoCA targets. Old result directories remain untouched.
