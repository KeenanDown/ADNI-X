# Analysis catalogue

Every analysis can run directly from preparation. Every listed figure can be regenerated from its numerical bundle. Cell numbers refer to the supplied notebook, starting at zero.

## Version 0.2.0 additions

| Analysis ID | Outputs |
|---|---|
| `group-matrices-moca` | MoCA group matrices; validated unique rtr/sts pair contrasts |
| `network-effects-moca` | All 16 atom × network effects for three MoCA contrasts |
| `taxa-principal-dx` | Upward/downward causation, transfer and storage by diagnosis |
| `taxa-principal-moca` | The same four modes by MoCA |
| `taxa-six-moca` | All six original modes by MoCA |

`coefficients` now includes `inference.csv` and per-atom `.significant` figure IDs. `taxa-four-mmse` and `taxa-six-mmse` remain callable but are opt-in for `run all`. `adni-x list` is the authoritative complete catalogue. The original cell mapping below records migration history; cell 72's coefficient-significance TODO is implemented in 0.2.0.

| Analysis ID | Notebook cells | Outputs |
|---|---|---|
| `group-matrices` | 59-65 | Group atomic matrices and edge contrasts |
| `coefficients` | 66-71 | Confound coefficient grids |
| `network-effects` | 73-94 | Intranetwork atom comparisons and effect sizes |
| `regional-effects` | 95-115,133-150 | Regional comparisons, significant counts and network summaries |
| `hierarchies` | 116-129,151-152,167-176 | Rank hierarchies, profiles and bootstrap correlations |
| `hierarchies-moca` | 130-132 | MoCA-group regional hierarchies |
| `taxa-four-dx` | 153-166 | Four taxa by DX |
| `taxa-four-mmse` | 153-166 | Four taxa by MMSE |
| `taxa-six-dx` | 153-166 | Six taxa by DX |
| `taxa-six-mmse` | 153-166 | Six taxa by MMSE |
| `cn-reference` | 177-209 | CN-standardised regional comparisons, PCA and brain exports |
| `lda` | 210-226 | Held-out CN/AD LDA and labelled projections |
| `heldout-hierarchy` | 227-238 | Held-out participant hierarchy correlations |
| `alterations` | 239-271,285-290 | Alteration correlations and exploratory network permutations |
| `longitudinal` | 272-281 | Raw all-visit atom and TDMI trajectories |
| `longitudinal-tdmi-ratio` | 282-284 | Optional raw atom/TDMI longitudinal sensitivity |
| `tdmi` | 291-296 | Raw global TDMI, diagnosis comparisons and length audit |

## Coverage rules

Cells 0–58 become preparation, validation and configuration. Cell 72 was an unimplemented coefficient-significance TODO in v1 and is implemented in 0.2.0. Repeated function definitions, inspection expressions and commented duplicate versions are consolidated; they are not extra independent scientific results. Cell 291 becomes explicit scalar metadata exports. External MATLAB surface rendering remains outside the package, but its CSV/MAT input exports are produced by `cn-reference`.

The original notebook is not executed or shipped with participant outputs. `notebook_cells.json` records its content hash and assigns every cell to a migrated analysis, preparation, or a stated non-analysis category. Numerical outputs are not used as correctness targets where the notebook contains bugs.

## Figure IDs

### `group-matrices`

Dependency: adjusted selected cohort. 

- `group-matrices.rtr`
- `group-matrices.sts`

### `coefficients`

Dependency: adjusted selected cohort. 

- `coefficients.rtr`
- `coefficients.sts`

### `network-effects`

Dependency: adjusted selected cohort. 

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

### `regional-effects`

Dependency: adjusted selected cohort. 

- `regional-effects.rtr.MCI-vs-CN`
- `regional-effects.rtr.AD-vs-CN`
- `regional-effects.rtr.AD-vs-MCI`
- `regional-effects.sts.MCI-vs-CN`
- `regional-effects.sts.AD-vs-CN`
- `regional-effects.sts.AD-vs-MCI`
- `regional-effects.rtr.network-counts`
- `regional-effects.sts.network-counts`

### `hierarchies`

Dependency: adjusted selected cohort. 

- `hierarchies.hierarchy`
- `hierarchies.rtr`
- `hierarchies.sts`
- `hierarchies.rtr_deviation`
- `hierarchies.sts_deviation`
- `hierarchies.rtr_midpoint`
- `hierarchies.sts_midpoint`
- `hierarchies.correlations`

### `hierarchies-moca`

Dependency: adjusted selected cohort. Exploratory. 

- `hierarchies-moca.hierarchy`
- `hierarchies-moca.rtr`
- `hierarchies-moca.sts`
- `hierarchies-moca.rtr_deviation`
- `hierarchies-moca.sts_deviation`
- `hierarchies-moca.rtr_midpoint`
- `hierarchies-moca.sts_midpoint`
- `hierarchies-moca.correlations`

### `taxa-four-dx`

Dependency: adjusted selected cohort. 

- `taxa-four-dx.panels`

### `taxa-four-mmse`

Dependency: adjusted selected cohort. 

- `taxa-four-mmse.panels`

### `taxa-six-dx`

Dependency: adjusted selected cohort. 

- `taxa-six-dx.panels`

### `taxa-six-mmse`

Dependency: adjusted selected cohort. 

- `taxa-six-mmse.panels`

### `cn-reference`

Dependency: adjusted selected cohort. 

- `cn-reference.regions`
- `cn-reference.subjects`
- `cn-reference.MCI.networks`
- `cn-reference.AD.networks`
- `cn-reference.MCI.pca`
- `cn-reference.AD.pca`

### `lda`

Dependency: selected raw regional features; training-only adjustment. Exploratory. 

- `lda.confusion`
- `lda.projections`

### `heldout-hierarchy`

Dependency: selected raw regional features; training-only adjustment. Exploratory. 

- `heldout-hierarchy.correlations`

### `alterations`

Dependency: adjusted selected cohort. Exploratory. 

- `alterations.rtr.MCI-vs-CN`
- `alterations.rtr.AD-vs-CN`
- `alterations.rtr.AD-vs-MCI`
- `alterations.sts.MCI-vs-CN`
- `alterations.sts.AD-vs-CN`
- `alterations.sts.AD-vs-MCI`

### `longitudinal`

Dependency: all eligible raw visits. Exploratory. 

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

### `longitudinal-tdmi-ratio`

Dependency: all eligible raw visits; explicit ratio normalisation. Exploratory. Explicit opt-in. 

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

### `tdmi`

Dependency: selected raw scans with diagnosis. 

- `tdmi.groups`
- `tdmi.lengths`

