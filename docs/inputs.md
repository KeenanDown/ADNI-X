# Input contracts

## Tensors and scan matching

The default recursive tensor patterns are `*232_desc-phiidtensor.mat` (MMI) and `*232_desc-CCS-phiidtensor.mat` (CCS), with the parcel count substituted from configuration. `tensor_glob` is an explicit override. A custom pattern must isolate one estimator and parcellation.

Each tensor contains a real numeric `PhiIDTensor` with logical shape `(P, P, 16)`. Atom order is:

`rtr, rtx, rty, rts, xtr, xtx, xty, xts, ytr, ytx, yty, yts, str, stx, sty, sts`.

MATLAB v5–v7.2 files use SciPy. MATLAB v7.3 files use HDF5 and reverse storage dimensions according to MATLAB's format, then apply the **explicit** `tensor_axes` permutation. Generic HDF5 arrays without MATLAB metadata are rejected. Shape cannot establish scientific source/target orientation, especially for square matrices: verify this once against S8b code or a known tensor. No symmetry or implicit scientific-axis reversal is applied.

Subject/session identifiers come from BIDS-style `sub-…` / `ses-…` entities. Directory and filename identifiers must agree. The matcher also compares every available `task`, `acq`, `run`, `echo`, `rec` and `dir` entity. Missing entities permit a match only when one compatible candidate exists; the audit marks a partial identity match. Distinct eligible tensors for the same subject/visit remain ambiguous until explicitly resolved.

Original series default to `*ExtSchaefer232_ROI_activity.mat`, excluding Gaussian-copula and correlation-index companions. Exactly one real numeric 2-D matrix must have exactly one dimension equal to P. Specify `series_variable` if multiple arrays qualify. A square P×P matrix cannot establish a time axis and is rejected.

## Acquisition repetition time

Every tensor requires its adjacent `.json` sidecar with a `source_bold` string. For a tensor symlink, the sidecar is read beside its resolved target. Resolve relative `source_bold` paths against that sidecar's directory; absolute paths are also accepted. Replace the linked `.nii.gz` or `.nii` suffix with `.json` and read its numeric `RepetitionTime` in **seconds**. BOLD image data are not loaded. There is no inferred TR, alternate-scan lookup, NIfTI-header fallback or BIDS inheritance search. If relocated metadata still point to inaccessible sources, fix the copied input layout before running.

Both `cohort.tr_min_seconds: 2.875` and `cohort.tr_max_seconds: 3.150` are inclusive. Comparisons use the recorded value without rounding or tolerance. Both length and TR eligibility precede latest/baseline/specified-visit selection and determine which visits enter longitudinal analyses. A scan with 197 time points and TR 6.015 seconds passes the length check but is excluded as `above_tr_window`.

The inventory preserves `TR_seconds`, `tensor_json`, `source_bold`, `bold_json`, `tr_status` and `tr_detail`; retained cohort tables carry these columns into result bundles. `tr_status: ok` means the metadata were readable and valid, **not** that TR passed the configured range. `reason` records the first eligibility failure, with separate `below_tr_window` and `above_tr_window` reasons. Missing/unreadable sidecars, missing/invalid TR and conflicting source identifiers have distinct reasons. Positive finite JSON numbers are required; strings and booleans are rejected. Missing metadata never silently pass the filter. Available subject/session/run/acquisition entities must not conflict between tensor and source BOLD.

Both JSON files are content-hashed as preparation dependencies; missing/unreadable metadata states are recorded too. Changing, removing or repairing sidecars, changing `source_bold`, or changing TR bounds invalidates affected preparation/results. No raw NIfTI hash or access is needed for this check. This enforces the specified TR requirement; it does not independently verify every acquisition-protocol parameter.

## Clinical records

Required join columns: `PTID`, `VISCODE`. `001_S_0001` and `ADNI001S0001` map to the same subject; `bl`, `m12`, `M012` map to canonical visits. Unsupported visit codes are reported. Exactly repeated clinical records are accepted and counted; conflicting records are excluded unless a resolution identifies the row.

Analysis fields are `DX`, `PTGENDER`, `PTEDUCAT`, `AGE`, `Years_bl`, `ICV`, `MMSE`, `MOCA`, `EXAMDATE`. `Years.bl` is accepted as an alias. Absent fields remain missing; they are not fabricated. `Dementia` is displayed as `AD`, while the original `DX` is retained. Unknown sex categories are missing. MoCA scores outside integer 0–30 are not binned. MoCA categories are 0–19, 20–28 and 29–30. MMSE is retained only as source metadata and is never categorised or used to group analyses.

FreeSurfer recovery searches `sub-*/stats/aseg.stats` for EstimatedTotalIntraCranialVol. Only a positive, unambiguous subject-level estimate is used; no arbitrary value is chosen from conflicting estimates. The source of every usable ICV is recorded.

## Atlas mapping

Supply a CSV with **one row per tensor parcel, in exact tensor order**:

| parcel_index | parcel_id | name | network_id |
|---|---|---|---|
| 0 | atlas-specific stable ID | human-readable name | integer 1–8 |

`parcel_index` must be exactly `0 … P−1`; parcel IDs must be unique. Networks 1–7 are Visual, Somatomotor, Dorsal Attention, Salience/Ventral Attention, Limbic, Control and Default; 8 is Subcortical. A cortical-only Schaefer resource cannot establish subcortical identities. Use the actual S8b parcellation lookup to supply those rows. The package checks consistency, not anatomical truth.

`scripts/build_atlas.py` converts a Schaefer centroid CSV and an explicitly ordered subcortical CSV to the canonical format. For cortical-only atlases omit `--subcortex`. This supports 100, 116, 200 and 232 parcels without assuming the unnamed final parcels are interchangeable. Synthetic tests use smaller explicit mappings.

## Manual resolution table

Optional CSV columns:

- `tensor_file` (required): exact tensor path; relative paths resolve against the resolution CSV.
- `series_file`: choose one compatible original series; conflicting scan identifiers cannot be overridden.
- `clinical_row`: zero-based data-row number from the original clinical CSV (excluding header).
- `use_scan`: `false` to exclude that tensor explicitly. Exclude unwanted same-visit runs so one remains.

Unknown tensor references and duplicate resolution rows are errors. Manual choices are part of the source fingerprint and audit. Do not edit S8b inputs to resolve an ambiguity.

## Preparation artifacts

The cache contains `inventory.csv`, `visits.csv`, `selected.csv`, `subjects.csv`, the ordered atlas, raw/adjusted disk-backed arrays, reduced regional/network arrays, design coefficients and references, normalisation counts and provenance. No pickle is required. Arrays are float64. Global and within-network means use the defined connection mask; a missing included edge propagates into its affected aggregate rather than quietly changing its denominator.

Source files are content-hashed. Cache validation reads checksums, so a cache hit still performs input/integrity I/O, but does not repeat regression or analysis. This deliberately prioritises reproducibility over trusting filename dates.
