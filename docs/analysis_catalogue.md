# Analysis catalogue

`run all` includes 26 independent diagnosis/MoCA targets. Two additional TDMI-ratio targets remain opt-in. MMSE targets have been removed. MoCA uses bands 29–30, 20–28, 0–19 and the 29–30 band as reference. Every listed figure can be regenerated from its saved numerical bundle.

| Analysis ID | Grouping | Notebook cells (zero-based) | Outputs |
|---|---|---|---|
| `group-matrices` | diagnosis | 59-65 | Group atomic matrices and edge contrasts |
| `group-matrices-moca` | moca | 59-65 | Group atomic matrices and edge contrasts (MoCA categories) |
| `coefficients` | diagnosis | 66-71 | Confound coefficient grids |
| `coefficients-moca` | moca | 66-71 | Confound coefficient grids (MoCA categories) |
| `network-effects` | diagnosis | 73-94 | Intranetwork atom comparisons and effect sizes |
| `network-effects-moca` | moca | 73-94 | Intranetwork atom comparisons and effect sizes (MoCA categories) |
| `regional-effects` | diagnosis | 95-115,133-150 | Regional comparisons, significant counts and network summaries |
| `regional-effects-moca` | moca | 95-115,133-150 | Regional comparisons, significant counts and network summaries (MoCA categories) |
| `hierarchies` | diagnosis | 116-129,151-152,167-176 | Rank hierarchies, profiles and bootstrap correlations |
| `hierarchies-moca` | moca | 116-129,151-152,167-176,130-132 | Rank hierarchies, profiles and bootstrap correlations (MoCA categories) |
| `taxa-four-dx` | diagnosis | 153-166 | Four taxa by DX |
| `taxa-four-moca` | moca | 153-166 | Four taxa by MoCA (MoCA categories) |
| `taxa-six-dx` | diagnosis | 153-166 | Six taxa by DX |
| `taxa-six-moca` | moca | 153-166 | Six taxa by MoCA (MoCA categories) |
| `cn-reference` | diagnosis | 177-209 | CN-standardised regional comparisons, PCA and brain exports |
| `moca-reference` | moca | 177-209 | MoCA 29–30-standardised regional comparisons, PCA and brain exports (MoCA categories) |
| `lda` | diagnosis | 210-226 | Held-out CN/AD LDA and labelled projections |
| `lda-moca` | moca | 210-226 | Held-out MoCA 29–30/0–19 LDA and labelled projections (MoCA categories) |
| `heldout-hierarchy` | diagnosis | 227-238 | Held-out participant hierarchy correlations |
| `heldout-hierarchy-moca` | moca | 227-238 | Held-out participant hierarchy correlations (MoCA categories) |
| `alterations` | diagnosis | 239-271,285-290 | Alteration correlations and exploratory network permutations |
| `alterations-moca` | moca | 239-271,285-290 | Alteration correlations and exploratory network permutations (MoCA categories) |
| `longitudinal` | diagnosis | 272-281 | Raw all-visit atom and TDMI trajectories |
| `longitudinal-moca` | moca | 272-281 | Raw all-visit atom and TDMI trajectories (MoCA categories) |
| `longitudinal-tdmi-ratio` | diagnosis | 282-284 | Optional raw atom/TDMI longitudinal sensitivity |
| `longitudinal-tdmi-ratio-moca` | moca | 282-284 | Optional raw atom/TDMI longitudinal sensitivity (MoCA categories) |
| `tdmi` | diagnosis | 291-296 | Raw global TDMI, diagnosis comparisons and length audit |
| `tdmi-moca` | moca | 291-296 | Raw global TDMI, MoCA comparisons and length audit (MoCA categories) |

## Coverage rules

Cells 0–58 become preparation, validation and configuration. Cell 72 is an unimplemented coefficient-significance TODO and remains out of v1. Repeated function definitions, inspection expressions and commented duplicate versions are consolidated; they are not extra independent scientific results. Cell 291 becomes explicit scalar metadata exports. External MATLAB surface rendering remains outside the package, but its CSV/MAT input exports are produced by `cn-reference` and `moca-reference`.

The original notebook is not executed or shipped with participant outputs. `notebook_cells.json` records its content hash and assigns every cell to a migrated analysis, preparation, or a stated non-analysis category. Numerical outputs are not used as correctness targets where the notebook contains bugs.

Historical MMSE-based taxa cells map to the replacement MoCA target and are explicitly superseded; no MMSE-category calculation is retained. Other MoCA targets extend the corresponding diagnostic calculations.

## Figure IDs

### `group-matrices`

Dependency: adjusted selected cohort. Grouping: diagnosis.

- `group-matrices.rtr`
- `group-matrices.sts`

### `group-matrices-moca`

Dependency: adjusted selected cohort. Grouping: moca.

- `group-matrices-moca.rtr`
- `group-matrices-moca.sts`

### `coefficients`

Dependency: adjusted selected cohort. Grouping: diagnosis.

- `coefficients.rtr`
- `coefficients.sts`

### `coefficients-moca`

Dependency: adjusted selected cohort. Grouping: moca.

- `coefficients-moca.rtr`
- `coefficients-moca.sts`

### `network-effects`

Dependency: adjusted selected cohort. Grouping: diagnosis.

- `network-effects.MCI-vs-CN.cohens_d.all-atoms`
- `network-effects.MCI-vs-CN.cohens_d.rtr-sts`
- `network-effects.MCI-vs-CN.hedges_g.all-atoms`
- `network-effects.MCI-vs-CN.hedges_g.rtr-sts`
- `network-effects.AD-vs-CN.cohens_d.all-atoms`
- `network-effects.AD-vs-CN.cohens_d.rtr-sts`
- `network-effects.AD-vs-CN.hedges_g.all-atoms`
- `network-effects.AD-vs-CN.hedges_g.rtr-sts`
- `network-effects.AD-vs-MCI.cohens_d.all-atoms`
- `network-effects.AD-vs-MCI.cohens_d.rtr-sts`
- `network-effects.AD-vs-MCI.hedges_g.all-atoms`
- `network-effects.AD-vs-MCI.hedges_g.rtr-sts`

### `network-effects-moca`

Dependency: adjusted selected cohort. Grouping: moca.

- `network-effects-moca.20-28-vs-29-30.cohens_d.all-atoms`
- `network-effects-moca.20-28-vs-29-30.cohens_d.rtr-sts`
- `network-effects-moca.20-28-vs-29-30.hedges_g.all-atoms`
- `network-effects-moca.20-28-vs-29-30.hedges_g.rtr-sts`
- `network-effects-moca.0-19-vs-29-30.cohens_d.all-atoms`
- `network-effects-moca.0-19-vs-29-30.cohens_d.rtr-sts`
- `network-effects-moca.0-19-vs-29-30.hedges_g.all-atoms`
- `network-effects-moca.0-19-vs-29-30.hedges_g.rtr-sts`
- `network-effects-moca.0-19-vs-20-28.cohens_d.all-atoms`
- `network-effects-moca.0-19-vs-20-28.cohens_d.rtr-sts`
- `network-effects-moca.0-19-vs-20-28.hedges_g.all-atoms`
- `network-effects-moca.0-19-vs-20-28.hedges_g.rtr-sts`

### `regional-effects`

Dependency: adjusted selected cohort. Grouping: diagnosis.

- `regional-effects.rtr.MCI-vs-CN`
- `regional-effects.rtr.AD-vs-CN`
- `regional-effects.rtr.AD-vs-MCI`
- `regional-effects.sts.MCI-vs-CN`
- `regional-effects.sts.AD-vs-CN`
- `regional-effects.sts.AD-vs-MCI`
- `regional-effects.rtr.network-counts`
- `regional-effects.sts.network-counts`

### `regional-effects-moca`

Dependency: adjusted selected cohort. Grouping: moca.

- `regional-effects-moca.rtr.20-28-vs-29-30`
- `regional-effects-moca.rtr.0-19-vs-29-30`
- `regional-effects-moca.rtr.0-19-vs-20-28`
- `regional-effects-moca.sts.20-28-vs-29-30`
- `regional-effects-moca.sts.0-19-vs-29-30`
- `regional-effects-moca.sts.0-19-vs-20-28`
- `regional-effects-moca.rtr.network-counts`
- `regional-effects-moca.sts.network-counts`

### `hierarchies`

Dependency: adjusted selected cohort. Grouping: diagnosis.

- `hierarchies.hierarchy`
- `hierarchies.rtr`
- `hierarchies.sts`
- `hierarchies.rtr_deviation`
- `hierarchies.sts_deviation`
- `hierarchies.rtr_midpoint`
- `hierarchies.sts_midpoint`
- `hierarchies.correlations`

### `hierarchies-moca`

Dependency: adjusted selected cohort. Grouping: moca.

- `hierarchies-moca.hierarchy`
- `hierarchies-moca.rtr`
- `hierarchies-moca.sts`
- `hierarchies-moca.rtr_deviation`
- `hierarchies-moca.sts_deviation`
- `hierarchies-moca.rtr_midpoint`
- `hierarchies-moca.sts_midpoint`
- `hierarchies-moca.correlations`

### `taxa-four-dx`

Dependency: adjusted selected cohort. Grouping: diagnosis.

- `taxa-four-dx.panels`

### `taxa-four-moca`

Dependency: adjusted selected cohort. Grouping: moca.

- `taxa-four-moca.panels`

### `taxa-six-dx`

Dependency: adjusted selected cohort. Grouping: diagnosis.

- `taxa-six-dx.panels`

### `taxa-six-moca`

Dependency: adjusted selected cohort. Grouping: moca.

- `taxa-six-moca.panels`

### `cn-reference`

Dependency: adjusted selected cohort. Grouping: diagnosis.

- `cn-reference.regions`
- `cn-reference.subjects`
- `cn-reference.MCI.networks`
- `cn-reference.AD.networks`
- `cn-reference.MCI.pca`
- `cn-reference.AD.pca`

### `moca-reference`

Dependency: adjusted selected cohort. Grouping: moca.

- `moca-reference.regions`
- `moca-reference.subjects`
- `moca-reference.20-28.networks`
- `moca-reference.0-19.networks`
- `moca-reference.20-28.pca`
- `moca-reference.0-19.pca`

### `lda`

Dependency: selected raw regional features; training-only adjustment. Grouping: diagnosis.

- `lda.confusion`
- `lda.projections`

### `lda-moca`

Dependency: selected raw regional features; training-only adjustment. Grouping: moca.

- `lda-moca.confusion`
- `lda-moca.projections`

### `heldout-hierarchy`

Dependency: selected raw regional features; training-only adjustment. Grouping: diagnosis.

- `heldout-hierarchy.correlations`

### `heldout-hierarchy-moca`

Dependency: selected raw regional features; training-only adjustment. Grouping: moca.

- `heldout-hierarchy-moca.correlations`

### `alterations`

Dependency: adjusted selected cohort. Grouping: diagnosis.

- `alterations.rtr.MCI-vs-CN`
- `alterations.rtr.AD-vs-CN`
- `alterations.rtr.AD-vs-MCI`
- `alterations.sts.MCI-vs-CN`
- `alterations.sts.AD-vs-CN`
- `alterations.sts.AD-vs-MCI`

### `alterations-moca`

Dependency: adjusted selected cohort. Grouping: moca.

- `alterations-moca.rtr.20-28-vs-29-30`
- `alterations-moca.rtr.0-19-vs-29-30`
- `alterations-moca.rtr.0-19-vs-20-28`
- `alterations-moca.sts.20-28-vs-29-30`
- `alterations-moca.sts.0-19-vs-29-30`
- `alterations-moca.sts.0-19-vs-20-28`

### `longitudinal`

Dependency: all eligible raw visits. Grouping: diagnosis.

- `longitudinal.rtr.calendar_year`
- `longitudinal.rtr.years_since_first_eligible_scan`
- `longitudinal.rtr.years_since_ADNI_baseline`
- `longitudinal.rtr.visit_age`
- `longitudinal.rtr.rolling`
- `longitudinal.sts.calendar_year`
- `longitudinal.sts.years_since_first_eligible_scan`
- `longitudinal.sts.years_since_ADNI_baseline`
- `longitudinal.sts.visit_age`
- `longitudinal.sts.rolling`
- `longitudinal.TDMI.calendar_year`
- `longitudinal.TDMI.years_since_first_eligible_scan`
- `longitudinal.TDMI.years_since_ADNI_baseline`
- `longitudinal.TDMI.visit_age`
- `longitudinal.TDMI.rolling`

### `longitudinal-moca`

Dependency: all eligible raw visits. Grouping: moca.

- `longitudinal-moca.rtr.calendar_year`
- `longitudinal-moca.rtr.years_since_first_eligible_scan`
- `longitudinal-moca.rtr.years_since_ADNI_baseline`
- `longitudinal-moca.rtr.visit_age`
- `longitudinal-moca.rtr.rolling`
- `longitudinal-moca.sts.calendar_year`
- `longitudinal-moca.sts.years_since_first_eligible_scan`
- `longitudinal-moca.sts.years_since_ADNI_baseline`
- `longitudinal-moca.sts.visit_age`
- `longitudinal-moca.sts.rolling`
- `longitudinal-moca.TDMI.calendar_year`
- `longitudinal-moca.TDMI.years_since_first_eligible_scan`
- `longitudinal-moca.TDMI.years_since_ADNI_baseline`
- `longitudinal-moca.TDMI.visit_age`
- `longitudinal-moca.TDMI.rolling`

### `longitudinal-tdmi-ratio`

Dependency: all eligible raw visits; explicit ratio normalisation. Grouping: diagnosis. Opt-in sensitivity.

- `longitudinal-tdmi-ratio.rtr.calendar_year`
- `longitudinal-tdmi-ratio.rtr.years_since_first_eligible_scan`
- `longitudinal-tdmi-ratio.rtr.years_since_ADNI_baseline`
- `longitudinal-tdmi-ratio.rtr.visit_age`
- `longitudinal-tdmi-ratio.rtr.rolling`
- `longitudinal-tdmi-ratio.sts.calendar_year`
- `longitudinal-tdmi-ratio.sts.years_since_first_eligible_scan`
- `longitudinal-tdmi-ratio.sts.years_since_ADNI_baseline`
- `longitudinal-tdmi-ratio.sts.visit_age`
- `longitudinal-tdmi-ratio.sts.rolling`

### `longitudinal-tdmi-ratio-moca`

Dependency: all eligible raw visits; explicit ratio normalisation. Grouping: moca. Opt-in sensitivity.

- `longitudinal-tdmi-ratio-moca.rtr.calendar_year`
- `longitudinal-tdmi-ratio-moca.rtr.years_since_first_eligible_scan`
- `longitudinal-tdmi-ratio-moca.rtr.years_since_ADNI_baseline`
- `longitudinal-tdmi-ratio-moca.rtr.visit_age`
- `longitudinal-tdmi-ratio-moca.rtr.rolling`
- `longitudinal-tdmi-ratio-moca.sts.calendar_year`
- `longitudinal-tdmi-ratio-moca.sts.years_since_first_eligible_scan`
- `longitudinal-tdmi-ratio-moca.sts.years_since_ADNI_baseline`
- `longitudinal-tdmi-ratio-moca.sts.visit_age`
- `longitudinal-tdmi-ratio-moca.sts.rolling`

### `tdmi`

Dependency: selected raw scans with analysis-group labels. Grouping: diagnosis.

- `tdmi.groups`
- `tdmi.lengths`

### `tdmi-moca`

Dependency: selected raw scans with analysis-group labels. Grouping: moca.

- `tdmi-moca.groups`
- `tdmi-moca.lengths`

Additional presentation-only coefficient layouts: `coefficients.combined` and `coefficients-moca.combined`.
