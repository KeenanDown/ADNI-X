# Methods and explicit choices

## Cohorts and transformations

Length eligibility is checked before selecting the maximum canonical visit month. Confound completeness is checked afterwards. Thus a newer long scan can be replaced by an older length-eligible scan, but missing confounds in the selected eligible scan do not cause a silent older-visit substitution. All eligible visits are retained for longitudinal use. Multiple eligible same-visit runs are excluded until resolved.

The primary adjusted cohort includes CN/MCI/AD and complete confounds. Raw TDMI retains selected diagnosed participants without requiring confound completeness. Missing MMSE/MoCA affects only the corresponding cognitive analysis. These cohort distinctions are exported with every result.

Visit age is `AGE + Years_bl`; ADNI confirms that `AGE` is baseline age: https://adni.loni.usc.edu/support/experts-knowledge-base/question/?QID=596 .

For numeric response features Y and confounds C, fit `Y = b0 + (C − mean(C)) B + error`. Default adjustment is `Y − (C − mean(C)) B`. All means, imputations and coefficients are fitted on the applicable fitting cohort. The centred intercept equals its mean response when a feature is fully observed. Binary sex coding uses Female as reference. Unidentifiable coefficients are flagged through design rank; the least-squares prediction remains defined, but individual collinear coefficients must not be interpreted.

Default complete-case handling performs no cohort imputation. The optional `missing: median` sensitivity setting fills each numeric variable using its own fitting-cohort median and uses an explicit missing-sex indicator. It records every preparation imputation. Predictive fitting uses training-only medians and references.

For optional TDMI division, compute the sum of all 16 **raw** atoms per connection first, divide each atom by that sum, then adjust confounds. Denominators with absolute value below `epsilon` or nonfinite values are undefined; the invalid-edge counts are saved. Regression uses feature-specific valid rows where necessary. Downstream sums propagate invalid values.

Regional sums use axis 1, excluding self-connections. Intranetwork means include both off-diagonal directions. Unique-pair matrix and coefficient tests validate rtr/sts symmetry for every included participant (matching finite masks; rtol=1e-7, atol=1e-10). A failed check fails the analysis rather than dropping a direction. Other atom reductions retain their directed definitions. CN reference scaling uses per-region/per-atom mean and population SD across CN participants. Zero reference variance produces NaN with a saved validity mask. CN scaling does not enable TDMI division.

PCA is two-dimensional on group-mean CN-standardised RTR and STS. PC1 is oriented toward increasing STS, with RTR as tie-breaker; scores retain original parcel order. Degenerate eigenspaces are reported instead of assigning an arbitrary axis. MATLAB exports include **parcel IDs**, not network labels masquerading as parcel IDs.

## Inference

Welch tests are the default; signed t statistics and effect sizes always describe **B minus A**. Cohen's d uses pooled within-group SD and Hedges' g uses `1 − 3/(4df − 1)`. A constant or undersized sample receives a reason and undefined statistics, not an invented zero effect. Welch ANOVA is the primary omnibus test; raw TDMI also saves ordinary ANOVA as sensitivity. Eta squared is the descriptive sums-of-squares effect size, not a Welch-specific estimator.

BH-FDR is applied to finite p-values within fixed families:

| Result | Correction family |
|---|---|
| Edge matrices | All unique included rtr/sts ROI pairs, per atom and contrast |
| Confound coefficients | All non-intercept coefficients × unique ROI pairs, separately for rtr and sts |
| Regional effects | All parcels, per atom and contrast |
| Intranetwork atoms | All 16 atoms × all available networks, per contrast |
| Network means of regional sums | Networks, per atom and contrast |
| Taxon omnibus | Every taxon × network panel, per taxonomy and grouping |
| Taxon post-hoc | All panel × pairwise contrasts, independently of omnibus selection |
| Network regional-shuffle tests | Networks, per atom and contrast |
| Alteration correlations | All atom/contrast correlations in that analysis |

RTR/STS-only figures reuse the full intranetwork correction family. All effects and p/q values are exported, whether significant or not. Undefined tests do not contaminate finite FDR calculations. Display-only ordering never changes the test family.

For 232 regions, unique off-diagonal matrix families contain 26,796 pairs per atom/contrast; for 116 regions they contain 6,670. Diagonal-inclusive sensitivity analyses add each self-pair once. Full-matrix significance plots mirror unique-pair results. Removing exact duplicate hypotheses does not automatically strengthen BH significance: duplicated identical p-values can produce the same q-values.

Coefficient tests use the original (optionally TDMI-normalised) responses and existing centred confound-only regression. For each feature-valid participant subset, HC3 uses covariance `X+ diag((e/(1-h))²) X+ᵀ`, where `X+` is the design pseudoinverse and `h` is leverage. Two-sided approximate Student-t tests use df=n−rank(X); 95% confidence intervals use the same reference distribution. The intercept remains descriptive and is not tested. Rank-deficient designs, nonpositive residual degrees of freedom, unit leverage and zero/nonfinite standard errors produce explicit unavailable statuses. These are pooled covariate associations conditional on the other covariates; no diagnosis terms are added. Existing group tests remain conditional on fitted adjustment.

The principal taxonomy contains upward causation, downward causation, transfer and storage as separate sums of their original atoms. It does not impose U=D or any copy/erasure/transfer identity. All six original modes remain available. Principal omnibus correction covers all four modes × networks; post-hoc correction covers all modes × networks × contrasts independently of omnibus selection, separately within each grouping/configuration.

Diagnosis and MoCA analyses run by default. MoCA categories are 29–30, 20–28 and 0–19. Score missingness excludes a participant only from that cognitive analysis; scores are not imputed. MMSE targets are optional. Deviations in regional profile plots are group mean minus the per-parcel CN mean (29–30 reference for MoCA), without SD division. Parcels are ordered by AD−CN (0–19 minus 29–30 for MoCA). Network panels retain that relative order and share the overview's scale. These profiles and their ordering are descriptive.

Bootstrap hierarchy intervals and alteration intervals resample participants within diagnosis groups, conditional on the already fitted descriptive confound adjustment. Regional network nulls shuffle regional changes jointly across networks, centred at the global regional mean. They assume exchangeable regions and **do not preserve spatial autocorrelation**; treat them as exploratory. Diagnosis-label permutations for alteration/baseline correlations preserve their shared-baseline coupling, using the empirical null mean for a two-sided test. These null hypotheses are separate and never implemented by shuffling ICV.

Random streams derive from seed 42 and a stable analysis/contrast identifier. Full null/bootstrap arrays, seed and iteration count are saved. Five-resample and ten-permutation notebook experiments are replaced by 10,000 by default; smaller values are for smoke tests only.

## Predictive and longitudinal outputs

CN/AD LDA uses the notebook's SVD solver and stratified 75/25 participant split. Confound adjustment and finite/varying-feature selection are fitted on training participants only. Saved scores include training, held-out and MCI projection-only labels. Accuracy, balanced accuracy, AD ROC-AUC and confusion matrices use held-out participants only. The LDA decision score includes its fitted intercept; its positive class is recorded. Held-out hierarchy analysis builds its CN reference from training CN participants only. See https://scikit-learn.org/stable/common_pitfalls.html .

Longitudinal outputs are descriptive, using raw atoms and raw TDMI, plus an opt-in atom/TDMI branch. They distinguish calendar year, time since first eligible scan, ADNI `Years_bl`, and visit age. Bins are left-closed with widths 1 year or 5 years of age, including the maximum value and baseline-only datasets. Average within participant/bin before averaging across participants, and export visit/participant counts. No independence-based confidence intervals are imposed on repeated scans. The exploratory five-bin rolling mean is labelled separately.

## Explicit historical options

Set these only in a separately named sensitivity configuration: `transform.age: baseline`, `adjustment: intercept` (legacy uncentred intercept addback) or `residuals`, `exclude_diagonal: false`, `regional_axis: 0`, `network_edges: upper`, `statistics.welch: false`, or `statistics.correction: none`. These settings remain visible in provenance. `network_edges: upper` changes directed information summaries and is not a claim of symmetry. No single legacy switch reproduces arbitrary stale notebook state or restores known coding errors.

