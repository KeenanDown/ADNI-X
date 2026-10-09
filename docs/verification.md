# Local verification

## Combined HPC/TR/MoCA merge — 9 October 2026

- The full WSL suite passed **297 tests**. A final focused check passed **7 tests** for independent-model HC3 inference, optional MMSE rendering and MoCA coefficient maps.
- All **32 targets / 160 registered figures** are retained; `run all` selects **28** diagnosis/MoCA targets, with two MMSE and two TDMI-ratio targets optional.
- Fresh actual-size synthetic runs completed all default targets for 232-MMI, 232-CCS, 116-MMI and 116-CCS. Both diagnosis and MoCA HC3 tables contained the expected 214,368 or 53,360 tests, with valid statuses and validated unique matrix pairs.
- The incoming TR acquisition code, cohort filter, independent grouping preparation, runner and TR tests were preserved unchanged. Eligibility still requires inclusive TR 2.875–3.150 seconds before visit selection, with missing/invalid metadata excluded.
- The notebook has **110** unique cells, valid Python and no saved outputs. Principal-mode sections and optional MMSE sections coexist with all incoming diagnosis/MoCA sections.
- The synthetic rerun used `--no-figures`, keeping only an aggregate summary outside the code repository. Historical records below describe their respective earlier code states.

## Diagnosis and MoCA update — 2026-09-24

**258 tests passed**. All **28 numerical targets** (26 default plus two opt-in TDMI-ratio targets) completed on synthetic data, and all **152 registered figure recipes** rendered from saved numerical bundles. Tests include missing/invalid MoCA, exact score boundaries, discordant diagnosis/MoCA labels, independent cohort membership and confound coefficients, no older-visit substitution for missing scores, mean preservation, separate FDR families, training-only prediction preprocessing, longitudinal visit labels/time origins, MATLAB identifiers/parcel order, standalone/all equivalence, rejection of removed MMSE target IDs, and source-free MoCA SVG/PNG/PDF export.

The replacement notebook's 98 cells have unique IDs, parseable Python and no stored outputs; all 28 individual targets and the default paired run are represented. Source distribution and wheel builds succeeded. The earlier 297-cell notebook catalogue remains complete, with historical MMSE grouping explicitly superseded by MoCA.

This verification uses synthetic fixtures only. No cluster jobs or real-data analyses were run for this update. Existing configurations need no additional grouping fields; the implementation fingerprint creates new preparation/results directories. Diagnosis and MoCA have separate adjusted stores and a shared raw store.

## Acquisition TR update — 2026-09-22

**167 tests passed** after adding the default inclusive 2.875–3.150-second repetition-time filter. New coverage verifies exact boundaries (and immediately adjacent floating-point values), the 197-point/6.015-second outlier, exclusion before latest/baseline/specified-visit selection, exclusion from all retained visits, relative/absolute source links, `.nii`/`.nii.gz` metadata lookup, tensor symlinks, missing/invalid/unreadable metadata, conflicting source identities, invalid configuration and cache invalidation after metadata/settings changes. Synthetic fixtures now include S8b source links and BOLD acquisition JSONs. The full existing numerical and saved-result rendering suite also passed.

No new real-data run or HPC changes were performed for this update. Existing runs must be recomputed to apply the new cohort rule; source-free replotting preserves the cohort of the original saved result.

## Initial implementation — 2026-09-16

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
