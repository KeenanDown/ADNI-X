# Local verification — 2026-09-16

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
