# Local verification — 2026-09-16

## Release 0.2.0 verification — 9 October 2026

- Full WSL test suite: **185 passed**, using Python 3.12.3 and the existing locked dependency versions.
- **22 analysis targets / 105 registered figures**; default runs include 19 targets, with MMSE and the longitudinal ratio branch optional.
- Synthetic integration runs at actual 232/116 matrix sizes completed all default targets for **232-MMI, 232-CCS, 116-MMI and 116-CCS**, using 18 synthetic participants and eight networks per configuration.
- Unique matrix families contained **26,796** or **6,670** pairs per atom/contrast. Coefficient tables contained **214,368** or **53,360** non-intercept tests per configuration (two atoms × four covariates × unique pairs).
- Batched HC3 standard errors, p-values and confidence intervals matched an independent full-sandwich calculation. Original and TDMI-normalised responses, missing feature rows, invalid regression designs, correction families and mirrored displays were checked.
- CN/MoCA-baseline centring, endpoint-difference ordering, saved-result compatibility, MCI–CN Cohen's d exports and explicit MoCA score boundaries/missingness were tested. Notebook code parses, IDs are unique and participant outputs are cleared.
- After visual review corrected a profile-note/legend collision, the affected rendering and rerun suites passed **132 checks**, including a non-overlap assertion.
- Aggregate integration and wheel-install checks are in `reports/hpc-rerun-0.2.0/`. The synthetic example SVG/PDF/PNG exports were inspected and removed after verification; `scripts/verify_rerun_release.py` can regenerate them. These tensors are format/integration fixtures, not estimates computed from neural signals; their group results have no scientific interpretation.

The older verification record below is retained as history. No real-data cluster jobs have been run for this release. See [the HPC handoff](hpc-rerun-0.2.0.md).

The implementation was verified on macOS arm64, Python 3.12.14, using `requirements.lock`.

| Check | Result |
|---|---|
| Automated tests | **134 passed** |
| Independent numerical analysis targets | **17 / 17 completed** on synthetic data |
| Figure recipes | **86 / 86 rendered** from saved numerical bundles |
| Full image export | **172 files:** 86 SVG + 86 PNG |
| PNG resolution | All verified at **300 dpi** |
| Package build | Source distribution and universal Python wheel built successfully |
| Installed command | `adni-x --version`, catalogue and dependency consistency checks passed |
| Notebook | Code cells parse successfully; stable cell IDs; no saved participant outputs |
| Notebook coverage | All **297 source cells** classified and mapped in `notebook_cells.json` |

Tests cover length boundaries and exclusion-before-latest selection, ambiguous/missing series, conflicting clinical records and explicit resolution, MATLAB storage formats and axis declarations, asymmetric tensors, exact synthetic regression coefficients, all confound-toggle combinations, mean preservation, training-only preprocessing, FDR families and boundary handling, deterministic randomisation, empty adjusted cohorts, ICV recovery, column-wise sensitivity imputation, cache invalidation/corruption, interrupted writes, failure isolation, standalone/all consistency, and result-only rendering.

Every registered renderer was drawn from its saved result bundle in the test suite. A separate end-to-end smoke run exported all SVGs and PNGs. Representative group-matrix and taxon layouts were visually inspected. Package distribution contents were checked for the expected code, fixtures and handoff files.

The synthetic cohort contains 24 participants across CN/MCI/AD, two visits each, and four explicitly labelled parcels. Its numerical results have no scientific interpretation. Small bootstrap/permutation counts are used in tests; production defaults remain 10,000.

## Still requires the HPC

The real-data pilot and full ADNI rerun have **not** been performed. The open HPC terminal was not accessible through the available computer-use tool. No cluster inputs were inspected or modified during implementation and no cluster jobs were submitted.

Before scientific use, verify actual S8b atom/axis conventions, the full atlas including subcortical identities, original-series matching, clinical joins, exclusion counts, confound design rank, and Linux environment compatibility. Follow `hpc.md`; these checks cannot be certified by synthetic tests.
