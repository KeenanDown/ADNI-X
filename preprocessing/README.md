# Upstream preprocessing and PhiID tensor generation

This directory preserves custom scripts supplied from the author's ADNI-PID HPC
checkout on 29 October 2026. They precede ADNI-X's `prepare` and `run` commands.
They are not a new implementation or an independently validated end-to-end
reproduction of the manuscript. No participant data, tensors, atlas binaries,
credentials, or third-party toolbox code are included.

## Files and processing order

| Stage | Source | Role |
|---|---|---|
| S3 | `hpc/S3-ADNItoBIDS.sh` | Historical Clinica ADNI-to-BIDS invocation. |
| S4 | `hpc/S4-ADNIfMRIPrep.sh` | Historical fMRIPrep/Apptainer invocation. |
| S5 | `hpc/S5-RegressFilter.sh`, `RegressFilter.py` | Confound regression and band-pass filtering using Nilearn. |
| Parcellation | `parcellate.m` | Average voxels within each atlas label to obtain parcel-by-time series. |
| MMI tensors | `create_phiID_tensors_v3.m` | Call external `PhiIDFull` with MMI for each off-diagonal parcel pair. |
| MMI HPC jobs | `hpc/S8bx-CreatePhiIDTensor-116.sh`, `hpc/S8bx-CreatePhiIDTensor-232.sh` | Historical SGE wrappers for the two resolutions. |
| Earlier generation | `historical/create_phiID_tensors.m`, `historical/create_phiID_tensors_v2.m`, `hpc/S8b-CreatePhiIDTensor.sh` | Earlier default-estimator implementations; see the CCS limitation below. |

The supplied S8b wrapper copies its input directly from S5 output. This is
consistent with the author's report that S6 deconvolution was skipped. S8bx
subsequently operates on that copied BOLD directory. S6/S7 and experimental
Gaussian-copula scripts are not included in this workflow.

The `hpc/` files are historical SGE examples, not ready-to-submit jobs. Site paths
have been replaced with `/path/to/...` placeholders and the notification email
removed. Their original copy/overwrite and work-directory cleanup operations
remain. Use working copies of inputs, adapt the script locations and MATLAB
path, and verify inputs before submitting jobs. The old S3 project date and S4
container filename are retained as evidence, not asserted as the final cohort's
provenance.

`source-manifest.json` records original and distributed SHA-256 hashes and every
sanitization. Numerical logic is preserved; line endings and trailing whitespace
have been normalized where necessary. The
v2 documentation example uses an anonymous placeholder instead of a participant
identifier. These source snapshots do not establish historical execution by
themselves.

## External dependencies and inputs

- **elph / PhiID** is an external dependency, not vendored here. Supply the actual
  checkout containing `PhiID/PhiIDFull.m`, its `DoubleRedundancyMMI` and
  `DoubleRedundancyCCS` helpers, `octaveToJavaDoubleMatrix`, and the compatible
  JIDT `infodynamics.jar`. Its upstream URL, historical commit and any local
  changes still need recording; no unverified upstream version is prescribed.
- **MATLAB** must provide `niftiread`, `niftiinfo`, and the statistical functions
  used by elph, including `mvnpdf`, with Java support for JIDT. The historical
  MATLAB/toolbox and JIDT versions remain to be recorded.
- **Nilearn, nibabel, NumPy and pandas** are needed by `RegressFilter.py`.
  The historical NiPy environment was not version-pinned in the supplied setup
  script. ADNI-X's root `requirements.lock` covers downstream analysis only;
  it does not specify the preprocessing environment.
- **Clinica, fMRIPrep, Apptainer and a separately obtained FreeSurfer license**
  are needed if starting before S5. Use execution reports to identify the
  historical versions, not the current scripts alone.
- **Atlas resources:** supply the actual extended 116- and 232-label NIfTI maps
  used for processing, with matching voxel grids, label order and source
  attribution. The historical filenames are
  `Schaefer2018_116Parcels_7Networks_order_FSLMNI152_2mm.nii` and
  `Schaefer2018_232Parcels_7Networks_order_FSLMNI152_2mm.nii`.
  Their files/checksums and label mappings are not included yet. A standard
  cortical-only Schaefer atlas is not an interchangeable substitute for these
  extended atlases. `parcellate.m` assumes nonzero labels are consecutive 1..P.
- **Data:** obtain ADNI inputs through the appropriate access process. The
  release does not redistribute scans, confound tables or participant lists.

## What the available code does

`RegressFilter.py` selects available columns from a requested set of five CSF
CompCor and five WM CompCor components, six motion parameters and their first
derivatives, and `non_steady_state_outlier00`. It fills missing confound values
from the next row, then calls `nilearn.image.clean_img` with `low_pass=0.05`,
`high_pass=0.0025`, the supplied TR in seconds, and `standardize=False`.
Filtering and nuisance regression are delegated to that single call. No sample
mask or explicit temporal-resampling step is supplied. Filter order, detrending
defaults and internal operation order depend on the historical Nilearn version.
The S5 wrapper overwrites the copied `preproc_bold` file; its filename alone
does not distinguish cleaned from uncleaned data.

The MMI generator calls `parcellate.m` directly and then
`PhiIDFull(two_region_series, [], 'MMI')`. The supplied HPC `PhiIDFull.m`, inspected
separately, defaults an empty lag to one sample, scales series by their standard
deviations, and evaluates Gaussian log densities with natural logarithms. It
contains no rank-to-normal transformation of the input series. Its default
redundancy estimator is CCS. The separate Gaussian-copula branch writes `-gcop`
tensor files; the saved manuscript run manifests instead list ordinary MMI and
CCS tensor filenames. These observations support the ordinary branch, but are
not a complete audit of every scan or every historical toolbox revision.

## Single-scan examples

These examples illustrate the supplied script interfaces. Resolve dependency
versions first; neither command has been run here on participant data. From this
directory, S5 can be invoked with separate input and output files:

```sh
python RegressFilter.py /path/to/input_bold.nii.gz \
  /path/to/working-copy/scan_desc-preproc_bold.nii.gz \
  /path/to/confounds.tsv 3.0 32G
```

Replace `3.0` with that acquisition's actual `RepetitionTime`. The last argument
is retained for compatibility with the historical wrapper; setting `PYTHONMEM`
inside this script is not itself a memory allocator or limit. Preserve the
matching BOLD JSON beside the output for downstream TR lookup.

In MATLAB, for a working copy of the cleaned BOLD file:

```matlab
addpath('/path/to/ADNI-X/preprocessing');
create_phiID_tensors_v3( ...
    '/path/to/working-copy/scan_desc-preproc_bold.nii.gz', ...
    '/path/to/elph', ...
    '/path/to/Schaefer2018_232Parcels_7Networks_order_FSLMNI152_2mm.nii', ...
    'ExtSchaefer232');
```

Use the matching 116 atlas and `ExtSchaefer116` for the other resolution. Outputs
are written beside the BOLD file and existing outputs may be overwritten. The
historical MATLAB functions catch errors and print them; inspect the job log
and verify newly produced outputs rather than relying on the shell exit status.

## Connecting to ADNI-X

The MMI output is `*_parc-P_desc-phiidtensor.mat`, containing `PhiIDTensor` with
shape P-by-P-by-16, plus a JSON sidecar with `source_bold`, parcel count and atom
order. `parcellate.m` also writes `*_ExtSchaeferP_ROI_activity.mat`. Point the
ADNI-X configuration at the tensors, original series, clinical metadata and
ordered atlas CSV as described in [input contracts](../docs/inputs.md).
The tensor sidecar's `source_bold` must remain accessible with its corresponding
BOLD JSON containing `RepetitionTime`; correct paths if relocating inputs.
The v3 function's generic `parcellationMatrix.mat` companion is not the
resolution-specific series file selected by ADNI-X.

The earlier generators use `PhiIDFull`'s default (CCS in the supplied toolbox
copy). V2 writes an unsuffixed tensor name, whereas ADNI-X expects
`*_parc-P_desc-CCS-phiidtensor.mat`; the historical renaming step is not recorded
here. Keep CCS and MMI outputs separate to avoid overwriting either estimator.

## Version records to add

Record the elph upstream URL/commit and changes, JIDT and MATLAB versions, the
preprocessing environment and actual atlas checksums/label mappings, and verify
fMRIPrep versions from execution reports. In particular,
the supplied current S4 wrapper names a 24.1.1 container whereas the manuscript
previously specified 23.0.2. That difference requires historical records and
must not be resolved by silently changing the script or manuscript to match.

Only source integrity, Python parsing and shell syntax were checked when this
directory was added. No MRI preprocessing, MATLAB execution, numerical
equivalence test or complete cohort rerun was performed.
