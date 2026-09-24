# Methods and explicit choices

## Cohorts and transformations

Length eligibility (150–250 retained time points inclusive) and acquisition repetition time (2.875–3.150 seconds inclusive) are checked before selecting the maximum canonical visit month. TR is read from the BOLD JSON identified by each S8b tensor's `source_bold` metadata. Missing or invalid metadata excludes the scan. The TR requirement was added following the user's post hoc acquisition review, which identified a 197-point scan with TR 6.015 seconds. This is an acquisition eligibility change, not a change to PhiID computation.

Confound completeness is checked afterwards. Thus a newer scan failing either length or TR eligibility can be replaced by an older scan satisfying both rules, but missing confounds in the selected eligible scan do not cause a silent older-visit substitution. All eligible visits are retained for longitudinal use. Multiple eligible same-visit runs are excluded until resolved. As a consequence, applying this rule to all candidate visits may differ from simply removing one row from an existing selected cohort; review the new inclusion audit and cohort counts.

Every default analysis is run independently for diagnosis (CN/MCI/AD) and MoCA category (29–30, 20–28, 0–19). Only integer MoCA scores in 0–30 are categorised; missing, fractional and out-of-range values are excluded from MoCA analyses. These are the existing notebook MoCA bands, not diagnoses. MMSE may be retained as source metadata, but no MMSE category or MMSE-based analysis is calculated.

Each grouping's adjusted cohort requires its own valid label and complete confounds; its confound model is fitted separately. A valid MoCA score does not require a known diagnosis, and a missing MoCA score does not affect diagnostic eligibility. Missing labels are checked after visit selection and never cause substitution of an older visit. Raw TDMI requires the relevant group label without requiring confound completeness. These distinctions, model sizes and exclusion reasons are exported with every result. Raw tensors are shared; adjusted arrays and model artifacts are stored separately for diagnosis and MoCA.

Visit age is `AGE + Years_bl`; ADNI confirms that `AGE` is baseline age: https://adni.loni.usc.edu/support/experts-knowledge-base/question/?QID=596 .

For numeric response features Y and confounds C, fit `Y = b0 + (C − mean(C)) B + error`. Default adjustment is `Y − (C − mean(C)) B`. All means, imputations and coefficients are fitted on the applicable fitting cohort. The centred intercept equals its mean response when a feature is fully observed. Binary sex coding uses Female as reference. Unidentifiable coefficients are flagged through design rank; the least-squares prediction remains defined, but individual collinear coefficients must not be interpreted.

Default complete-case handling performs no cohort imputation. The optional `missing: median` sensitivity setting fills each numeric variable using its own fitting-cohort median and uses an explicit missing-sex indicator. It records every preparation imputation. Predictive fitting uses training-only medians and references.

For optional TDMI division, compute the sum of all 16 **raw** atoms per connection first, divide each atom by that sum, then adjust confounds. Denominators with absolute value below `epsilon` or nonfinite values are undefined; the invalid-edge counts are saved. Regression uses feature-specific valid rows where necessary. Downstream sums propagate invalid values.

Regional sums use axis 1, excluding self-connections. Intranetwork means include both off-diagonal directions; the program never assumes atom symmetry. Reference scaling uses per-region/per-atom mean and population SD across CN participants for diagnosis and across MoCA 29–30 participants for MoCA. Zero reference variance produces NaN with a saved validity mask. CN scaling does not enable TDMI division.

PCA is two-dimensional on group-mean reference-standardised RTR and STS, for MCI/AD in diagnosis and 20–28/0–19 in MoCA. MATLAB MoCA fields use valid identifiers such as `MOCA_20_28_rtr`, preserving parcel order. PC1 is oriented toward increasing STS, with RTR as tie-breaker; scores retain original parcel order. Degenerate eigenspaces are reported instead of assigning an arbitrary axis. MATLAB exports include **parcel IDs**, not network labels masquerading as parcel IDs.

## Inference

Welch tests are the default; signed t statistics and effect sizes always describe **B minus A**. Cohen's d uses pooled within-group SD and Hedges' g uses `1 − 3/(4df − 1)`. A constant or undersized sample receives a reason and undefined statistics, not an invented zero effect. Welch ANOVA is the primary omnibus test; raw TDMI also saves ordinary ANOVA as sensitivity. Eta squared is the descriptive sums-of-squares effect size, not a Welch-specific estimator.

BH-FDR is applied to finite p-values within fixed families:

| Result | Correction family |
|---|---|
| Edge matrices | All included directed edges, per atom and contrast |
| Regional effects | All parcels, per atom and contrast |
| Intranetwork atoms | All 16 atoms × all available networks, per contrast |
| Network means of regional sums | Networks, per atom and contrast |
| Taxon omnibus | Every taxon × network panel, per taxonomy and grouping |
| Taxon post-hoc | All panel × pairwise contrasts, independently of omnibus selection |
| Network regional-shuffle tests | Networks, per atom and contrast |
| Alteration correlations | All atom/contrast correlations in that analysis |

Correction families are separate for diagnosis and MoCA; no correction pools the two groupings. RTR/STS-only figures reuse the full intranetwork correction family. All effects and p/q values are exported, whether significant or not. Undefined tests do not contaminate finite FDR calculations. Display-only ordering never changes the test family.

Bootstrap hierarchy intervals and alteration intervals resample participants within the active analysis groups, conditional on the already fitted descriptive confound adjustment. Regional network nulls shuffle regional changes jointly across networks, centred at the global regional mean. They assume exchangeable regions and **do not preserve spatial autocorrelation**; treat them as exploratory. Analysis-group-label permutations for alteration/baseline correlations preserve their shared-baseline coupling, using the empirical null mean for a two-sided test. These null hypotheses are separate and never implemented by shuffling ICV.

Random streams derive from seed 42 and a stable analysis/contrast identifier. Full null/bootstrap arrays, seed and iteration count are saved. Five-resample and ten-permutation notebook experiments are replaced by 10,000 by default; smaller values are for smoke tests only.

## Predictive and longitudinal outputs

CN/AD LDA uses the notebook's SVD solver and stratified 75/25 participant split. Confound adjustment and finite/varying-feature selection are fitted on training participants only. Saved scores include training, held-out and MCI projection-only labels. Accuracy, balanced accuracy, AD ROC-AUC and confusion matrices use held-out participants only. The LDA decision score includes its fitted intercept; its positive class is recorded. Held-out hierarchy analysis builds its CN reference from training CN participants only. The MoCA counterpart compares 29–30 against 0–19; 20–28 participants are projection-only. Its probability and ROC-AUC refer to 0–19, while the decision score's positive class is recorded separately. Its held-out hierarchy reference uses training MoCA 29–30 participants only. Both groupings use separate deterministic random streams and fitted models. See https://scikit-learn.org/stable/common_pitfalls.html .

Longitudinal outputs are descriptive, using raw atoms and raw TDMI, plus an opt-in atom/TDMI branch. Diagnosis and MoCA labels are taken from each visit. Time since first eligible scan uses the first acquisition-eligible visit, even if its active grouping label is missing. They distinguish calendar year, time since first eligible scan, ADNI `Years_bl`, and visit age. Bins are left-closed with widths 1 year or 5 years of age, including the maximum value and baseline-only datasets. Average within participant/bin before averaging across participants, and export visit/participant counts. No independence-based confidence intervals are imposed on repeated scans. The exploratory five-bin rolling mean is labelled separately.

## Explicit historical options

Set these only in a separately named sensitivity configuration: `transform.age: baseline`, `adjustment: intercept` (legacy uncentred intercept addback) or `residuals`, `exclude_diagonal: false`, `regional_axis: 0`, `network_edges: upper`, `statistics.welch: false`, or `statistics.correction: none`. These settings remain visible in provenance. `network_edges: upper` changes directed information summaries and is not a claim of symmetry. No single legacy switch reproduces arbitrary stale notebook state or restores known coding errors.
