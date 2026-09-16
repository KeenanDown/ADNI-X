# ADNI-X

Reproducible ADNI analyses starting from **completed S8b PhiID tensors**. ADNI-X does not recompute PhiID or preprocess MRI scans.

The workflow is **prepare once → run independent analyses → replot saved numerical results**. A Jupyter kernel, MATLAB, and the original notebook's intermediate pickles are not required.

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

`run` produces numerical bundles and figures. Use `--no-plots` to compute only. The opt-in `longitudinal-tdmi-ratio` analysis is omitted by `all` unless `--include-sensitivity` is supplied; requesting that analysis by name is also explicit opt-in. All other implemented analyses are included, including exploratory ones.

To regenerate just one figure, supply its saved analysis directory or the parent run directory:

```sh
adni-x plot group-matrices.rtr --results /path/to/results/RUN_ID/group-matrices
adni-x plot all --results /path/to/results/RUN_ID
```

Replotting needs **only the result bundle**, not raw data or preparation caches. Numerical outputs have stable names inside immutable bundles. New SVG and 300-dpi PNG files have informative names and their actual UTC creation timestamp. Every plot has a JSON provenance sidecar.

## Defaults

- 232 parcels, MMI, no TDMI division.
- Inclusive **150–250 retained time points**, then latest eligible scan per participant. This excludes shorter scans too and is **not a verified Basic/Advanced protocol classifier**.
- Confounds: sex, education, visit age (`AGE + Years_bl`) and ICV. Missing ICV may be recovered from FreeSurfer; otherwise use complete cases after visit selection, with no older-visit substitution.
- Adjustment preserves the fitting cohort's mean. Self-connections are excluded; both off-diagonal directions are retained; regional sums use axis 1.
- Welch tests, two-sided α=.05, BH-FDR, seed 42, 10,000 resamples/permutations.
- Prediction uses a stratified 75/25 participant split with training-only preprocessing.

Raw TDMI uses all selected eligible participants with a diagnosis, including those lacking confounds. Adjusted analyses use their documented complete-case cohort. Longitudinal plots use all eligible visits. Each bundle contains cohort tables so these differences are visible.

## Python and notebook use

```python
from adni_x import prepare, run, plot, validate
from adni_x.config import Config

config = Config.load("config.local.yaml")
audit = validate(config)
prepared = prepare(config)
report = run("hierarchies", config)
figures = plot("hierarchies.hierarchy", report["path"])
```

`notebooks/ADNI_X.ipynb` is a thin entry point to this API. Numerical analysis functions are also available in `adni_x.analyses`, `adni_x.prediction`, and `adni_x.longitudinal`; they accept a prepared dataset and return a `Result` containing tables, numeric arrays, figure recipes and metadata.

## Documentation and verification

- [Analysis catalogue and notebook coverage](docs/analysis_catalogue.md)
- [Methods, correction families and legacy options](docs/methods.md)
- [Input and atlas contracts](docs/inputs.md)
- [Bug ledger](docs/bug_ledger.md)
- [HPC pilot and rerun guide](docs/hpc.md)
- [Local verification results and remaining HPC checks](docs/verification.md)

```sh
python -m pytest
python -m build --no-isolation
```

Tests use synthetic participants only. Real HPC tensor conventions, atlas identity and cohort counts must be checked in the pilot; passing local tests does not certify the real dataset. Dataset files, caches and participant-level results are ignored by Git. Do not commit them.
