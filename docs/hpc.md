# HPC pilot and complete rerun

## Readiness boundary

Local implementation and tests use synthetic data. No cluster job has been submitted and no HPC files have been changed. The real-data pilot is **not completed locally**. Do not interpret fixture cohort counts as real ADNI counts.

Copy/clone the package to your normal code directory on the HPC. Install the locked environment using Python 3.12, then make a local configuration from `configs/hpc.example.yaml`. Inputs must be read-only in normal operation. Cache/output directories must be separate writable locations with adequate scratch space.

## Before computing

1. Check the S8b estimator, atom order, source/target axes, and complete cortical/subcortical parcel order against the upstream implementation. Supply the full atlas CSV.
2. Confirm the original time-series variable and matching filename entities. Confirm that the 150–250-point rule is the intended operational filter.
3. Run `adni-x validate --config config.local.yaml`. Inspect `inventory.csv`, reasons, partial-identity matches, unknown diagnoses, conflicting clinical records and cohort counts. Resolve ambiguous files through the optional CSV rather than editing raw data.
4. Confirm available disk space. 2,907 × 232 × 232 × 16 float64 values are about 20 GB before adjusted copies, temporary work and result files. Blocks limit computation memory, not on-disk cache size. Start with one process; avoid multiplying independent full preparation jobs.

## Small real-data pilot

Use a separate pilot input directory containing symlinks to a deterministic, documented selection of approximately 10–20 participants per diagnosis and all their candidate visits. Include boundary lengths, an advanced-length scan, and any known ambiguous runs. Never copy a misleading latest-only subset: filtering must be exercised before selecting the latest eligible scan. Supply matching series symlinks and the same clinical/atlas sources; point a separate pilot YAML at these directories and separate cache/results paths.

Keep default methods and use `resamples: 100`, `permutations: 100` only for this pilot. These small counts are for runtime checking, not inference. Run:

```sh
adni-x validate --config config.pilot.yaml
adni-x run all --config config.pilot.yaml --include-sensitivity
```

Check `run.json` for failed and insufficient-data targets. Check source identities, selected visits, group counts, ICV sources, adjustment exclusions, feature validity, design rank and figure labels. Inspect several tensors and reductions against independently calculated values. Expected changes from the old notebook include the length rule, filter-before-latest selection, complete cases, visit age, removed diagonals, mean-preserving adjustment and corrected inference. Do not require corrected p-values to match old displayed output.

Move a copied result bundle to a directory with no access to its original inputs and run `adni-x plot all --results /path/to/bundle`. Confirm SVG/PNG rendering and the parcel-indexed MATLAB exports.

## Full run

Restore the full inputs and 10,000 resamples/permutations. Use the scheduler-neutral wrapper inside the HPC's approved compute-job mechanism:

```sh
sh scripts/run_hpc.sh /absolute/path/config.local.yaml
```

The wrapper does not select a scheduler, queue or account. It sets numerical-library thread counts before starting Python, logs validate/prepare/run phases, and returns a nonzero exit code on failure. Choose HPC resource limits based on the pilot's measured peak memory/time. Do not run the complete analysis on a login node.

For independent jobs after preparation, invoke `adni-x run ANALYSIS --config … --no-plots`; use separate job logs. Raw scans are never modified. Plotting can run later without inputs.

## Resume and diagnose

Repeat the same command to reuse valid preparation and completed numerical targets. Failed targets are recorded under `failures/` and retried; other targets continue. Insufficient-data targets are completed reports with explicit reasons, not fabricated figures. Changes to source contents, method settings or package source create a new cache/result identity.

Incomplete temporary directories are never treated as completed caches. A checksum failure is a hard error: move the affected bundle aside for inspection, then rerun. Do not edit a sealed result in place. Preserve the run manifest, dependency lock and code revision with the final analysis archive.

