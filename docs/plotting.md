# Manuscript figures

The default renderer follows the presentation in `ADNI_Epilogue (2).ipynb`:

- Magma group matrices with a shared scale, blue–white–red significant changes,
  top/left network strips and a shared network legend.
- Icefire coefficient maps, including a combined redundancy/synergy grid.
- Blue–white–orange intranetwork effect grids with sign annotations and horizontal
  network labels. Numerical effect sizes remain in `effects.csv`.
- Compact network-coloured regional significance strips, ordered by AD–CN.
- Hierarchy/profile scatter plots with a large overview and a network mini-grid.
- Grouped, notched taxon boxplots: one panel per taxon, networks on the x-axis,
  the original group palette, and stars from the saved omnibus tests. Outliers
  are hidden in these boxplots for visual correspondence, not removed from tests.
- CN-reference scatter/PCA panels and horizontal prediction/hierarchy violins.

The style changes do not change the numerical analyses, selected cohorts,
confound adjustment, edge masks, tests, correction families, bootstrap intervals,
train/test split or saved results. Plotters consume saved significance flags.
Thresholds, exploratory status and held-out versus training panels remain explicit.
Effect colour limits are symmetric and are not clipped at the old arbitrary cap.
Coefficient colour limits are also symmetric, so equal magnitudes have equal
colour intensity. Named anatomical labels are never inferred from generic IDs.

## Restyle an existing run

After updating the package, restart the notebook kernel. Set `RESULTS_RUN` in the
helper cell to the existing run directory, then run:

```python
refresh_figures("network-effects")  # one analysis
show_figures("network-effects")
# refresh_figures()                # every completed analysis in that saved run
```

Or use the Python API directly:

```python
from adni_x import plot
plot("all", "/path/to/results/RUN_ID")
plot("coefficients.combined", "/path/to/results/RUN_ID/coefficients")
```

The combined coefficient figure is derived from the arrays already present in old
bundles; no migration or recomputation is needed. Existing per-atom figure IDs
remain available. The notebook displays the combined grid when present.

Do not call `run` or `prepare` just to change the appearance. Package source changes
affect the computation cache identity, but `plot` works directly with an existing
bundle and does not need the source tensors or preparation cache.

## Editable exports and provenance

Every rendering creates new timestamped SVG, PDF and 300-dpi PNG files. SVG text
remains editable; PDF uses TrueType fonts. Fonts follow the original preference:
Helvetica, Arial, Liberation Sans, Nimbus Sans, then DejaVu Sans. Actual appearance
can vary with fonts installed on the rendering machine.

Each figure sidecar records the style version, plotting source hash, original
manifest hash and checksums of all three exports. Numerical manifests and old
figures are left intact. `show_figures` chooses the newest rendering per figure ID,
avoiding duplicate old/new versions in the notebook. Styles apply in a local
Matplotlib context rather than overriding subsequent user plots.

The supplied notebook did not contain rendered examples for every plot. Where
there was no directly comparable original panel, the renderer uses the same
fonts and palettes while preserving the corrected analysis's distinct meaning.
