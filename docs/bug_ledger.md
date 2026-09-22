# Notebook bug and method-change ledger

Cell numbers below are zero-based in the supplied 297-cell notebook. Historical displayed outputs may reflect earlier execution states and are not golden truth for corrected calculations.

| Location | Problem | Resolution and verification |
|---|---|---|
| 5,15,27,33 | Save helper ignores its argument; cache save/load filenames disagree and omit settings | Immutable content/config/code-keyed cache, manifest checksums, atomic publication; cache and interruption tests |
| 11,27,273 | Subject/session-only records discard source identity and arbitrarily drop duplicate scans | Preserve BIDS entities and exact paths; explicit ambiguous-run/clinical exclusions and resolution table tests |
| 39,295 | Length sensitivity follows latest selection rather than preceding it | Apply length eligibility before visit selection; boundary and fallback tests |
| 29,31 | Normalised tensors computed even when disabled; zero TDMI becomes a false zero | Default no division; opt-in before regression, mask/count undefined denominators |
| 45,49 | Mean filling across unlike confounds; chained assignment; categorical encoding fails when sex disabled | Complete-case default, unambiguous ICV recovery, optional column-wise fitted imputation; all confound combinations tested |
| 47,49 | Baseline age used for later visits | Default AGE + Years_bl; explicit baseline sensitivity |
| 49,51 | Object-array regression; implicit global model and unsafe intercept branch | Numeric multi-output block regression with explicit design; known-coefficient and mean-preservation tests |
| 51 | Adding intercept at zero age/ICV has an unintended reference level | Default mean-preserving centred adjustment; named historical intercept option |
| 53 | Mixed atlas sources and assumed subcortical order | Full ordered atlas with stable parcel IDs; reject mismatch, no inferred anatomy |
| 58,271 | Diagnosis overwritten with cognitive labels and restoration commented out | Separate immutable diagnosis and cognition columns |
| 63,75,97,122,157 | Diagonal, axis and triangle conventions disagree | Shared masks, axis-1 reduction, both directions; asymmetric-tensor tests |
| 77-115,135,163 | Different test/correction settings and signed effects; nonsignificant results lost | One comparison API, fixed FDR families, B−A convention, all results saved |
| 122-126 | Fresh unseeded generators despite global seed; five bootstraps | Stable per-target streams, 10,000 iterations, full distributions and repeatability tests |
| 179 | Zero CN SD replaced by 1 | Undefined standardisations flagged with validity mask |
| 195-202 | Shared PCA state, incomplete outputs, network IDs exported as regions | Group-specific persisted PCA, parcel IDs, deterministic sign and degeneracy handling |
| 215,224,229 | Unseeded splits, global preprocessing leakage and Series index misalignment | Participant splits before fitted transforms, explicit scan-indexed outputs and training-only tests |
| 242-287 | Mixed regional/label nulls, irrelevant ICV permutation, hidden in-function shuffling | Separate explicit nulls, observed statistic never shuffled, joint regional null and saved distributions |
| 273-282 | Longitudinal data rebuilt without consistent filtering; baseline age and first-visit time confused | Reuse eligible visits; retain both time origins and true visit age |
| 278-282 | Terminal bin edge can omit the latest observation; repeat visits overweight people | Explicit floor bins and participant-first aggregation; maximum-time and count tests |
| 291 | Excluding columns by the word 'tensor' still exports residual arrays/objects | Explicit scalar cohort and tidy numerical tables |
| Throughout | Results depend on prior cells, overwritten names and global plotting state; many figures unsaved | Independent registered analyses, saved numerical recipes, fresh figures and source-free replot tests |

Method changes are documented separately from coding bugs in `methods.md`. The agreed 150–250 range also excludes 140-point scans; it must not be described as verified Advanced-only exclusion.

On 2026-09-22, following the user's post hoc acquisition-metadata review, an inclusive 2.875–3.150-second TR requirement was added before visit selection. Length alone had admitted a 197-point scan with TR 6.015 seconds. TR now comes from the tensor sidecar's linked BOLD JSON; unresolved metadata excludes a scan. Both sidecars and configured bounds participate in cache identity. This is a cohort-method revision; existing saved results must be recomputed, not merely replotted.
