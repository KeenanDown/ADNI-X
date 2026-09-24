# Local verification

## Diagnosis and MoCA update — 2026-09-24

**258 tests passed**. All **28 numerical targets** (26 default plus two opt-in TDMI-ratio targets) completed on synthetic data, and all **152 registered figure recipes** rendered from saved numerical bundles. Tests include missing/invalid MoCA, exact score boundaries, discordant diagnosis/MoCA labels, independent cohort membership and confound coefficients, no older-visit substitution for missing scores, mean preservation, separate FDR families, training-only prediction preprocessing, longitudinal visit labels/time origins, MATLAB identifiers/parcel order, standalone/all equivalence, rejection of removed MMSE target IDs, and source-free MoCA SVG/PNG/PDF export.

The replacement notebook's 98 cells have unique IDs, parseable Python and no stored outputs; all 28 individual targets and the default paired run are represented. Source distribution and wheel builds succeeded. The earlier 297-cell notebook catalogue remains complete, with historical MMSE grouping explicitly superseded by MoCA.

This verification uses synthetic fixtures only. No cluster jobs or real-data analyses were run for this update. Existing configurations need no additional grouping fields; the implementation fingerprint creates new preparation/results directories. Diagnosis and MoCA have separate adjusted stores and a shared raw store.

## Acquisition TR update — 2026-09-22

**167 tests passed** after adding the default inclusive 2.875–3.150-second repetition-time filter. New coverage verifies exact boundaries (and immediately adjacent floating-point values), the 197-point/6.015-second outlier, exclusion before latest/baseline/specified-visit selection, exclusion from all retained visits, relative/absolute source links, `.nii`/`.nii.gz` metadata lookup, tensor symlinks, missing/invalid/unreadable metadata, conflicting source identities, invalid configuration and cache invalidation after metadata/settings changes. Synthetic fixtures now include S8b source links and BOLD acquisition JSONs. The full existing numerical and saved-result rendering suite also passed.

No new real-data run or HPC changes were performed for this update. Existing runs must be recomputed to apply the new cohort rule; source-free replotting preserves the cohort of the original saved result.

## Initial implementation — 2026-09-16

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
