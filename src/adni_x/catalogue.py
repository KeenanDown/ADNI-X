from dataclasses import dataclass
from functools import partial
from . import analyses
from .longitudinal import longitudinal
from .prediction import lda, heldout_hierarchy
from .constants import PAIRS


@dataclass(frozen=True)
class Analysis:
    id: str
    description: str
    cells: str
    function: object
    figures: tuple[str, ...]
    exploratory: bool = False
    opt_in: bool = False
    dependency: str = "adjusted selected cohort"


CATALOGUE = [
    Analysis("group-matrices", "Group atomic matrices and edge contrasts", "59-65", analyses.group_matrices,
             tuple(f"group-matrices.{a}" for a in ("rtr", "sts"))),
    Analysis("coefficients", "Confound coefficient grids", "66-71", analyses.coefficients,
             tuple(f"coefficients.{a}" for a in ("rtr", "sts"))),
    Analysis("network-effects", "Intranetwork atom comparisons and effect sizes", "73-94", analyses.network_effects,
             tuple(f"network-effects.{b}-vs-{a}.{e}.{s}" for a, b in PAIRS for e in ("cohens_d", "hedges_g") for s in ("all-atoms", "rtr-sts"))),
    Analysis("regional-effects", "Regional comparisons, significant counts and network summaries", "95-115,133-150", analyses.regional_effects,
             tuple([f"regional-effects.{atom}.{b}-vs-{a}" for atom in ("rtr", "sts") for a, b in PAIRS] +
                   [f"regional-effects.{atom}.network-counts" for atom in ("rtr", "sts")])),
    Analysis("hierarchies", "Rank hierarchies, profiles and bootstrap correlations", "116-129,151-152,167-176", analyses.hierarchies,
             tuple(f"hierarchies.{m}" for m in ("hierarchy", "rtr", "sts", "rtr_deviation", "sts_deviation", "rtr_midpoint", "sts_midpoint", "correlations"))),
    Analysis("hierarchies-moca", "MoCA-group regional hierarchies", "130-132", partial(analyses.hierarchies, grouping="MOCA_category"),
             tuple(f"hierarchies-moca.{m}" for m in ("hierarchy", "rtr", "sts", "rtr_deviation", "sts_deviation", "rtr_midpoint", "sts_midpoint", "correlations")), True),
    *[Analysis(f"taxa-{size}-{group}", f"{size.title()} taxa by {group.upper()}", "153-166",
               partial(analyses.taxa, six=size == "six", grouping="group" if group == "dx" else "MMSE_category"),
               (f"taxa-{size}-{group}.panels",)) for size in ("four", "six") for group in ("dx", "mmse")],
    Analysis("cn-reference", "CN-standardised regional comparisons, PCA and brain exports", "177-209", analyses.cn_reference,
             ("cn-reference.regions", "cn-reference.subjects", "cn-reference.MCI.networks", "cn-reference.AD.networks",
              "cn-reference.MCI.pca", "cn-reference.AD.pca")),
    Analysis("lda", "Held-out CN/AD LDA and labelled projections", "210-226", lda,
             ("lda.confusion", "lda.projections"), True, dependency="selected raw regional features; training-only adjustment"),
    Analysis("heldout-hierarchy", "Held-out participant hierarchy correlations", "227-238", heldout_hierarchy,
             ("heldout-hierarchy.correlations",), True, dependency="selected raw regional features; training-only adjustment"),
    Analysis("alterations", "Alteration correlations and exploratory network permutations", "239-271,285-290", analyses.alterations,
             tuple(f"alterations.{atom}.{b}-vs-{a}" for atom in ("rtr", "sts") for a, b in PAIRS), True),
    Analysis("longitudinal", "Raw all-visit atom and TDMI trajectories", "272-281", longitudinal,
             tuple(f"longitudinal.{m}.{a}" for m in ("rtr", "sts", "TDMI") for a in
                   ("calendar_year", "years_since_first_eligible_scan", "years_since_ADNI_baseline", "visit_age", "rolling")),
             True, dependency="all eligible raw visits"),
    Analysis("longitudinal-tdmi-ratio", "Optional raw atom/TDMI longitudinal sensitivity", "282-284", partial(longitudinal, ratios=True),
             tuple(f"longitudinal-tdmi-ratio.{m}.{a}" for m in ("rtr", "sts") for a in
                   ("calendar_year", "years_since_first_eligible_scan", "years_since_ADNI_baseline", "visit_age", "rolling")),
             True, True, "all eligible raw visits; explicit ratio normalisation"),
    Analysis("tdmi", "Raw global TDMI, diagnosis comparisons and length audit", "291-296", analyses.tdmi,
             ("tdmi.groups", "tdmi.lengths"), dependency="selected raw scans with diagnosis"),
]
ANALYSES = {item.id: item for item in CATALOGUE}
FIGURES = {figure: item.id for item in CATALOGUE for figure in item.figures}


def catalogue():
    return [{"id": a.id, "description": a.description, "notebook_cells_zero_based": a.cells,
             "figures": list(a.figures), "exploratory": a.exploratory, "opt_in": a.opt_in,
             "preparation_dependency": a.dependency} for a in CATALOGUE]
