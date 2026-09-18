"""Rendering depends exclusively on immutable, saved numerical result bundles."""
from __future__ import annotations

from pathlib import Path
import json
import math
import os
import re
import tempfile
from inspect import signature
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import seaborn as sns
from matplotlib.patches import Patch

from .constants import GROUPS, GROUP_COLORS, NETWORK_COLORS, NETWORKS, ATOMS
from .provenance import atomic_json, file_hash, stamp, utc_now, verify

# Reference: ADNI_Epilogue notebook. Styling only; inference is read from bundles.
STYLE_VERSION = "epilogue-1"
STYLE = {"font.family": "sans-serif",
         "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"],
         "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
         "axes.spines.top": False, "axes.spines.right": False,
         "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
         "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white"}
SHORT_NETWORKS = {1: "Visual", 2: "SomMot", 3: "DorsAttn", 4: "Salience/VA",
                  5: "Limbic", 6: "Control", 7: "Default", 8: "Subcortex"}
EFFECT_CMAP = LinearSegmentedColormap.from_list(
    "blue_white_orange", ["#008deb", "#ffffff", "#eb8d00"], N=256)


def clean_axes(ax):
    for spine in ax.spines.values():
        spine.set_visible(False)


def stars(value):
    if not np.isfinite(value):
        return ""
    return "***" if value < .001 else "**" if value < .01 else "*" if value < .05 else ""


def measure_label(measure):
    labels = {"hierarchy": "Rank difference (synergy − redundancy)",
              "rtr": "Regional redundancy", "sts": "Regional synergy"}
    if measure in labels:
        return labels[measure]
    atom, suffix = measure.split("_", 1)
    return f"{labels.get(atom, atom)}: deviation from " + (
        "three-group mean" if suffix == "deviation" else "endpoint-group midpoint")


def render(bundle, recipe):
    # Do not change the user's notebook-wide matplotlib settings.
    with plt.rc_context(STYLE):
        return _render(bundle, recipe)



class Bundle:
    def __init__(self, path):
        self.path = Path(path)
        self.manifest = json.loads((self.path / "manifest.json").read_text())
        self.details = json.loads((self.path / "details.json").read_text())

    def table(self, name):
        return pd.read_csv(self.path / f"{name}.csv")

    def array(self, name):
        return np.load(self.path / f"{name}.npy", mmap_mode="r", allow_pickle=False)


def no_data(ax, text="Insufficient data"):
    ax.text(.5, .5, text, ha="center", va="center", transform=ax.transAxes, wrap=True)
    ax.set_xticks([]); ax.set_yticks([])


def grid(n, width=4.5, height=3.5, columns=3):
    n = max(1, n)
    columns = min(columns, n)
    rows = math.ceil(n / columns)
    fig, axes = plt.subplots(rows, columns, figsize=(width*columns, height*rows), squeeze=False, layout="constrained")
    for ax in axes.flat[n:]:
        ax.set_visible(False)
    return fig, list(axes.flat[:n])


def color(group, index=0):
    return GROUP_COLORS.get(str(group), plt.get_cmap("tab10")(index % 10))


def matrix(ax, array, title, *, diverging=False, limits=None, cmap=None, colorbar=True, aspect="equal"):
    arr = np.asarray(array, float)
    valid = arr[np.isfinite(arr)]
    ax.set_title(title)
    if not len(valid):
        no_data(ax, "No finite values / no significant effects")
        return
    if limits is None:
        if diverging:
            bound = max(float(np.abs(valid).max()), 1e-12)
            limits = (-bound, bound)
        else:
            limits = (float(valid.min()), float(valid.max()))
    im = ax.imshow(arr, cmap=cmap if cmap is not None else ("bwr" if diverging else "magma"),
                   vmin=limits[0], vmax=limits[1], interpolation="nearest", aspect=aspect)
    ax.set_xticks([]); ax.set_yticks([])
    clean_axes(ax)
    if colorbar:
        ax.figure.colorbar(im, ax=ax, shrink=.8).outline.set_visible(False)
    return im


def bars(ax, networks, *, sides="bottom", offset=-.065):
    rgba = np.array([matplotlib.colors.to_rgba(NETWORK_COLORS[int(n)]) for n in networks])
    if sides == "matrix":
        top = ax.inset_axes([0, 1.006, 1, .018])
        top.imshow(rgba[None, :, :], aspect="auto", interpolation="nearest")
        top.axis("off")
        left = ax.inset_axes([-.024, 0, .018, 1])
        left.imshow(rgba[:, None, :], aspect="auto", interpolation="nearest")
        left.axis("off")
    else:
        bottom = ax.inset_axes([0, offset, 1, .025])
        bottom.imshow(rgba[None, :, :], aspect="auto", interpolation="nearest")
        bottom.axis("off")


def network_legend(fig, networks):
    ids = sorted(set(map(int, networks)))
    fig.legend(handles=[Patch(color=NETWORK_COLORS[n], label=SHORT_NETWORKS[n]) for n in ids],
               loc="outside lower center", ncol=min(8, len(ids)), fontsize=8, frameon=False)


def distributions(ax, df, value, groups=None):
    groups = groups or [g for g in GROUPS if g in df.group.unique()]
    if not groups:
        groups = sorted(df.group.dropna().unique())
    for gi, group in enumerate(groups):
        vals = df.loc[df.group == group, value].dropna().to_numpy(float)
        vals = vals[np.isfinite(vals)]
        if not len(vals):
            continue
        if len(vals) >= 2 and np.ptp(vals) > 0:
            violin = ax.violinplot(vals, positions=[gi], showmedians=True, showextrema=False)
            for body in violin["bodies"]:
                body.set_facecolor(color(group, gi)); body.set_alpha(.5)
        jitter = np.random.default_rng(42).uniform(-.12, .12, len(vals))
        ax.scatter(gi + jitter, vals, s=9, alpha=.4, color=color(group, gi))
        ax.scatter([gi], [np.mean(vals)], marker="D", s=35, c="white", edgecolors="black", zorder=4)
    ax.set_xticks(range(len(groups)), [f"{g}\n(n={df.loc[df.group == g, value].notna().sum()})" for g in groups])
    ax.set_ylabel(value.replace("_", " "))


def _render(bundle, recipe):
    kind = recipe["kind"]
    correction = bundle.manifest["configuration"]["statistics"]["correction"]
    threshold_name = "BH-FDR" if correction == "fdr_bh" else "uncorrected p threshold"
    if kind == "group_matrices":
        atom = recipe["atom"]
        atlas = bundle.table("atlas")
        order = np.argsort(atlas.network_id, kind="stable")
        means = np.array(bundle.array("group_means")[:, :, :, atom])
        if bundle.manifest["configuration"]["transform"]["exclude_diagonal"]:
            means[:, np.arange(len(atlas)), np.arange(len(atlas))] = np.nan
        edges = bundle.table("edges")
        sub = edges.loc[(edges.atom == ATOMS[atom]) & (edges.group_a == "CN") & (edges.group_b == "AD")]
        difference = np.full((len(atlas), len(atlas)), np.nan)
        sig = sub.loc[sub.significant]
        difference[sig.parcel_i.to_numpy(int), sig.parcel_j.to_numpy(int)] = sig.difference
        fig = plt.figure(figsize=(12, 10), layout="constrained")
        layout = fig.add_gridspec(2, 4, width_ratios=[.045, 1, 1, .045], wspace=.10)
        axes = [fig.add_subplot(layout[r, c]) for r, c in ((0, 1), (0, 2), (1, 1), (1, 2))]
        finite = means[np.isfinite(means)]
        limits = (finite.min(), finite.max()) if len(finite) else None
        counts = bundle.table("group_counts").set_index("group")["n"].to_dict()
        first = None
        for gi, group in enumerate(GROUPS):
            im = matrix(axes[gi], means[gi][np.ix_(order, order)],
                        f"{group} (n={counts.get(group, 0)})", limits=limits, colorbar=False)
            if first is None and im is not None:
                first = im
            axes[gi].set_title(axes[gi].get_title(), pad=15)
            bars(axes[gi], atlas.network_id.to_numpy()[order], sides="matrix")
        im = matrix(axes[3], difference[np.ix_(order, order)],
                    f"Significant changes (AD − CN)\n{threshold_name}", diverging=True, colorbar=False)
        axes[3].set_title(axes[3].get_title(), pad=15)
        bars(axes[3], atlas.network_id.to_numpy()[order], sides="matrix")
        if first is not None:
            fig.colorbar(first, cax=fig.add_subplot(layout[:, 0])).outline.set_visible(False)
        if im is not None:
            fig.colorbar(im, cax=fig.add_subplot(layout[1, 3])).outline.set_visible(False)
        network_legend(fig, atlas.network_id)
        fig.suptitle("Redundancy" if atom == 0 else "Synergy", fontsize=14)
    elif kind == "coefficients":
        # Both atom slices are already saved in every coefficients bundle.
        # Keep individual atom IDs, and offer the original side-by-side layout.
        values = bundle.array("coefficients")
        names = bundle.details["names"]
        atoms = (0, 15) if recipe.get("combined") else (recipe["atom"],)
        fig, axes = plt.subplots(len(names), len(atoms), figsize=(4*len(atoms), 2.65*len(names)),
                                 squeeze=False, layout="constrained")
        atlas = bundle.table("atlas")
        order = np.argsort(atlas.network_id, kind="stable")
        labels = {"sex_Male": "Male − female", "education": "Education (per year)",
                  "visit_age": "Visit age (per year)", "ICV": "ICV (per volume unit)"}
        for row, name in enumerate(names):
            for col, atom in enumerate(atoms):
                ax = axes[row, col]
                a = np.array(values[row, :, :, atom])[np.ix_(order, order)]
                if bundle.manifest["configuration"]["transform"]["exclude_diagonal"]:
                    np.fill_diagonal(a, np.nan)
                matrix(ax, a, ("Redundancy (rtr)" if atom == 0 else "Synergy (sts)") if row == 0 else "",
                       diverging=True, cmap=sns.color_palette("icefire", as_cmap=True))
                ax.set_title(ax.get_title(), pad=12)
                bars(ax, atlas.network_id.to_numpy()[order], sides="matrix")
                if col == 0:
                    ax.set_ylabel(labels.get(name, name), labelpad=14)
        network_legend(fig, atlas.network_id)
        fig.suptitle("Confound coefficients (unthresholded)", fontsize=13)
    elif kind == "effects_grid":
        df = bundle.table(recipe["table"])
        df = df.loc[(df.group_a == recipe["group_a"]) & (df.group_b == recipe["group_b"]) & df.atom.isin(recipe["atoms"])]
        nets = sorted(df.network_id.unique())
        pivot = df.pivot(index="atom", columns="network_id", values=recipe["effect"]).reindex(index=recipe["atoms"], columns=nets)
        significant = df.pivot(index="atom", columns="network_id", values="significant").reindex_like(pivot).fillna(False).to_numpy(bool)
        fig, axes = grid(1, width=10, height=8 if len(pivot) > 2 else 3.2, columns=1)
        ax = axes[0]
        masked = np.where(significant, pivot.to_numpy(), np.nan)
        im = matrix(ax, masked, "", diverging=True, cmap=EFFECT_CMAP, colorbar=False, aspect="auto")
        effect_name = "Hedges’ g" if recipe["effect"] == "hedges_g" else "Cohen’s d"
        if im is not None:
            cbar = fig.colorbar(im, ax=ax, shrink=.9)
            cbar.set_label(f"{effect_name} (SD units)")
            cbar.outline.set_visible(False)
        ax.set_xticks(range(len(nets)), [SHORT_NETWORKS[n] for n in nets], rotation=0, fontsize=9)
        ax.set_yticks(range(len(pivot)), pivot.index, rotation=0)
        ax.set(xlabel="Networks", ylabel="Atoms")
        ax.xaxis.labelpad = 30
        bars(ax, nets, offset=-.16 if len(pivot) <= 2 else -.065)
        ax.set_xticks(np.arange(-.5, len(nets), 1), minor=True)
        ax.set_yticks(np.arange(-.5, len(pivot), 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=.5)
        ax.tick_params(which="minor", length=0)
        for r, c in zip(*np.where(significant & np.isfinite(pivot.to_numpy()))):
            ax.text(c, r, "+" if pivot.iloc[r, c] > 0 else "−" if pivot.iloc[r, c] < 0 else "0",
                    ha="center", va="center", fontsize=12, weight="bold")
        ax.set_title(f"Significant intranetwork changes ({effect_name})\n"
                     f"{recipe['group_b']} − {recipe['group_a']} · {threshold_name}")
    elif kind == "significance_strip":
        df = bundle.table("effects")
        atom = df.loc[df.atom == recipe["atom"]]
        anchor = atom.loc[(atom.group_a == "CN") & (atom.group_b == "AD")].sort_values("difference", ascending=False, kind="stable")
        sub = atom.loc[(atom.group_a == recipe["group_a"]) & (atom.group_b == recipe["group_b"])].set_index("parcel_index").reindex(anchor.parcel_index)
        fig, axes = grid(1, width=12, height=1.5, columns=1)
        ax = axes[0]
        sig = sub.significant.fillna(False).to_numpy(bool)
        rgba = np.array([matplotlib.colors.to_rgba(NETWORK_COLORS[int(n)] if valid else "white")
                         for n, valid in zip(sub.network_id, sig)])
        ax.imshow(rgba[None, :, :], aspect="auto", interpolation="nearest")
        ax.set_xticks([]); ax.set_yticks([]); clean_axes(ax)
        ax.set_title(f"{recipe['atom']} · {recipe['group_b']} − {recipe['group_a']} · {threshold_name}")
        ax.set_xlabel("Parcels ordered by AD − CN; coloured = significant, white = not significant")
        network_legend(fig, sub.network_id)
    elif kind == "network_counts":
        df = bundle.table("network_counts")
        df = df.loc[df.atom == recipe["atom"]]
        panels = list(df.groupby(["group_a", "group_b"], sort=False))
        fig, axes = grid(len(panels))
        for ax, ((a, b), sub) in zip(axes, panels):
            nets = sorted(sub.network_id.unique())
            for offset, direction in ((-.18, "decrease"), (.18, "increase")):
                counts = sub.loc[sub.direction == direction].set_index("network_id").reindex(nets)
                ax.bar(np.arange(len(nets)) + offset, counts.n_significant, .36, label=direction)
            ax.set_xticks(range(len(nets)), [str(n) for n in nets]); ax.set_title(f"{b} − {a}")
            ax.set(xlabel="Network ID", ylabel="Significant parcels"); ax.legend(fontsize=8)
    elif kind == "profiles":
        df = bundle.table("profiles")
        measure = recipe["measure"]
        if recipe.get("order") == "difference":
            pivot = df.pivot(index="parcel_index", columns="group", values=measure.split("_")[0])
            order = (pivot[bundle.details["groups"][-1]] - pivot[recipe["anchor"]]).sort_values(ascending=False, kind="stable").index
            ordering = f"{bundle.details['groups'][-1]} − {recipe['anchor']}"
        else:
            order = df.loc[df.group == recipe["anchor"]].sort_values(measure, kind="stable").parcel_index
            ordering = recipe["anchor"]
        networks = sorted(df.network_id.unique())
        fig = plt.figure(figsize=(11, 8), layout="constrained")
        outer = fig.add_gridspec(1, 2, width_ratios=[1, 1.2], wspace=.18)
        left = outer[0].subgridspec(2, 1, height_ratios=[18, 1])
        main = fig.add_subplot(left[0])
        strip = fig.add_subplot(left[1])
        right = outer[1].subgridspec(math.ceil(len(networks)/2), 2, hspace=.15)
        panels = [fig.add_subplot(right[i//2, i%2]) for i in range(len(networks))]
        ordered_nets = df.loc[df.group == recipe["anchor"]].set_index("parcel_index").reindex(order).network_id
        for gi, group in enumerate(bundle.details["groups"]):
            sub = df.loc[df.group == group].set_index("parcel_index").reindex(order)
            main.scatter(range(len(sub)), sub[measure], label=group, color=color(group, gi),
                         s=20, alpha=.65, edgecolors="none")
            for ax, network in zip(panels, networks):
                panel = sub.loc[sub.network_id == network]
                ax.scatter(range(len(panel)), panel[measure], s=18, color=color(group, gi),
                           alpha=.65, edgecolors="none")
                ax.set_title(SHORT_NETWORKS[network], color=NETWORK_COLORS[network], fontsize=10)
                ax.set_xticks([]); ax.set_yticks([]); clean_axes(ax)
        main.set(title=f"{ordering} ordering (n={len(order)})", xlabel=f"Parcel rank ({ordering} order)",
                 ylabel=measure_label(measure))
        clean_axes(main)
        main.legend(frameon=False, fontsize=9)
        rgba = np.array([matplotlib.colors.to_rgba(NETWORK_COLORS[int(n)]) for n in ordered_nets])
        strip.imshow(rgba[None, :, :], aspect="auto", interpolation="nearest")
        strip.axis("off")
        network_legend(fig, networks)
    elif kind == "correlation_intervals":
        df = bundle.table("correlations")
        fig, axes = grid(1, width=6, height=4, columns=1)
        ax = axes[0]
        for i, row in df.iterrows():
            if np.isfinite(row.pearson_r):
                ax.scatter(i, row.pearson_r, color="black")
                if np.isfinite(row.ci_lower) and np.isfinite(row.ci_upper):
                    ax.vlines(i, row.ci_lower, row.ci_upper, color="black")
        ax.set_xticks(range(len(df)), [f"{r.group_a} / {r.group_b}" for _, r in df.iterrows()], rotation=20)
        ax.set(ylabel="Hierarchy Pearson r (95% participant bootstrap CI)", ylim=(-1.05, 1.05))
    elif kind == "taxa_panels":
        df, tests = bundle.table("observations"), bundle.table("omnibus")
        # JSON is key-sorted on disk; restore the reference's scientific display order.
        preferred = ("causation", "copy_erasure", "transfer", "storage", "upward_causation",
                     "downward_causation", "erasure", "copy")
        if len(bundle.details["taxonomy"]) == 6:
            preferred = ("upward_causation", "downward_causation", "erasure", "copy", "transfer", "storage")
        taxa = [t for t in preferred if t in bundle.details["taxonomy"]]
        taxa += [t for t in bundle.details["taxonomy"] if t not in taxa]
        nets, groups = sorted(df.network_id.unique()), bundle.details["groups"]
        fig, axes = grid(len(taxa), columns=math.ceil(len(taxa)/2), width=4.5, height=3.5)
        width = .65/max(1, len(groups))
        for ax, taxon in zip(axes, taxa):
            for gi, group in enumerate(groups):
                for ni, net in enumerate(nets):
                    vals = df.loc[(df.taxon == taxon) & (df.network_id == net) & (df.group == group), "value"]
                    vals = vals[np.isfinite(vals)].to_numpy()
                    if not len(vals):
                        continue
                    parts = ax.boxplot([vals], positions=[ni + (gi-(len(groups)-1)/2)*width],
                                       widths=width*.9, notch=True, showfliers=False, patch_artist=True,
                                       manage_ticks=False, medianprops={"color": ".25", "linewidth": 1},
                                       whiskerprops={"color": ".4"}, capprops={"color": ".4"})
                    parts["boxes"][0].set(facecolor=color(group, gi), edgecolor=".4", linewidth=.7)
            for ni, net in enumerate(nets):
                test = tests.loc[(tests.taxon == taxon) & (tests.network_id == net)].iloc[0]
                # Annotation uses saved corrected tests, never a recomputed ANOVA.
                if bool(test.significant):
                    ax.text(ni, 1.01, stars(test.q) or f"{test.q:.3g}", transform=ax.get_xaxis_transform(),
                            ha="center", va="bottom", fontsize=10)
            ax.set_xticks(range(len(nets)), [SHORT_NETWORKS[n] for n in nets], rotation=45, ha="right", fontsize=8)
            ax.set_xlim(-.6, len(nets)-.4)
            ax.set_title(taxon.replace("_", " ").title(), pad=22)
            ax.set_ylabel("Within-network mean")
            ax.spines["left"].set_visible(False)
        label = "q" if correction == "fdr_bh" else "p"
        method = "Welch" if bundle.manifest["configuration"]["statistics"]["welch"] else "One-way"
        fig.suptitle(f"Within-network means by taxon and group\n"
                     f"{method} ANOVA · {threshold_name}: * {label}<.05, ** {label}<.01, *** {label}<.001",
                     fontsize=12)
        fig.legend(handles=[Patch(facecolor=color(g, i), label=g) for i, g in enumerate(groups)],
                   loc="outside right center", frameon=False, title="Group")
    elif kind == "reference_scatter":
        df = bundle.table(recipe["table"])
        if "group" in recipe:
            df = df.loc[df.group == recipe["group"]]
        fig, axes = grid(1, width=6, height=5, columns=1)
        ax = axes[0]
        if "group" in recipe and "network_id" in df:
            for net, sub in df.groupby("network_id"):
                ax.scatter(sub.rtr, sub.sts, label=NETWORKS[net], c=NETWORK_COLORS[net], s=20, alpha=.75)
        else:
            for gi, (group, sub) in enumerate(df.groupby("group", sort=False)):
                ax.scatter(sub.rtr, sub.sts, label=group, color=color(group, gi), s=24, alpha=.4)
        ax.scatter([0], [0], color="gray", marker="x", s=28, zorder=5)
        ax.set(xlabel="Redundancy: CN z-score", ylabel="Synergy: CN z-score")
        ax.legend(fontsize=8, frameon=False)
        clean_axes(ax)
    elif kind == "pca":
        group = recipe["group"]
        df = bundle.table("regions").query("group == @group")
        pcs = bundle.table("pca_scores").query("group == @group")
        pcs = pcs.set_index("parcel_index").reindex(df.parcel_index).reset_index()
        model = next(m for m in bundle.details["pca"] if m["group"] == group)
        fig = plt.figure(figsize=(9.5, 5.2), layout="constrained")
        layout = fig.add_gridspec(1, 2, width_ratios=[3, 1.2], wspace=.2)
        axes = [fig.add_subplot(layout[0]), fig.add_subplot(layout[1])]
        if model["status"] != "ok":
            for ax in axes: no_data(ax, model["status"])
        else:
            centre, pc1 = np.array(model["centre"]), np.array(model["pc1"])
            valid = np.isfinite(pcs.PC1)
            colors = [NETWORK_COLORS[n] for n in df.network_id]
            axes[0].scatter(df.rtr, df.sts, c=colors, s=18)
            extent = pcs.loc[valid, "PC1"].agg(["min", "max"]).to_numpy()
            line = centre + extent[:, None]*pc1
            axes[0].plot(line[:, 0], line[:, 1], "--", color="gray")
            axes[0].set(xlabel="CN-reference redundancy", ylabel="CN-reference synergy",
                        title=f"{group} · PC1 explains {model['explained_variance']:.1%}")
            projection = centre + pcs.PC1.to_numpy()[:, None]*pc1
            axes[1].scatter(projection[:, 0], projection[:, 1], c=colors, s=18)
            axes[1].set_title("Flattened along PC1")
            for ax in axes:
                ax.set_aspect("equal", adjustable="box")
                clean_axes(ax)
            axes[1].set_xticks([]); axes[1].set_yticks([])
            network_legend(fig, df.network_id)
    elif kind == "network_intervals":
        df = bundle.table("networks")
        sub = df.loc[(df.atom == recipe["atom"]) & (df.group_a == recipe["group_a"]) & (df.group_b == recipe["group_b"])].sort_values("network_id")
        fig, axes = grid(1, width=8, height=4.5, columns=1)
        ax = axes[0]
        for i, (_, row) in enumerate(sub.iterrows()):
            ax.bar(i, row.observed, color=NETWORK_COLORS[row.network_id], alpha=.8)
            if np.isfinite(row.ci_lower) and np.isfinite(row.ci_upper):
                ax.vlines(i, row.ci_lower, row.ci_upper, color="black")
            if row.significant and np.isfinite(row.observed):
                ax.annotate("*", (i, max(row.observed, row.ci_upper) if np.isfinite(row.ci_upper) else row.observed), ha="center", va="bottom")
        ax.axhline(0, color=".5", linewidth=.8)
        ax.set_xticks(range(len(sub)), [SHORT_NETWORKS[n] for n in sub.network_id], rotation=30, ha="right")
        ax.set(ylabel=f"{recipe['atom']} regional-sum change", title=f"{recipe['group_b']} − {recipe['group_a']} · exploratory regional shuffle")
    elif kind == "group_distribution":
        df = bundle.table(recipe["table"])
        fig, axes = grid(1, width=6, height=5, columns=1)
        distributions(axes[0], df, recipe["value"])
        axes[0].legend(handles=[Line2D([], [], marker="D", linestyle="none", color="black",
                                      markerfacecolor="white", label="Group mean")], frameon=False)
        test = bundle.table("omnibus").query("primary == True").iloc[0]
        axes[0].set(ylabel="Raw global TDMI (input units)", title=f"{test.method} ANOVA · p={test.p:.3g}")
    elif kind == "lengths":
        audit, observations = bundle.table("length_audit"), bundle.table("observations")
        fig, axes = grid(2, columns=2, width=6, height=4.5)
        for gi, group in enumerate(GROUPS):
            sub = audit.loc[audit.group == group]
            counts = sub.groupby("n_time").size()
            axes[0].scatter(counts.index, counts.values, label=group, color=color(group, gi), s=25)
            pts = observations.loc[observations.group == group]
            axes[1].scatter(pts.n_time, pts.value, label=group, color=color(group, gi), s=14, alpha=.5)
        c = bundle.manifest["configuration"]["cohort"]
        axes[0].axvspan(c["time_min"], c["time_max"], color="gray", alpha=.15)
        axes[0].set(xlabel="Retained time points", ylabel="Scans", title="Inventory length audit (protocol unverified)")
        axes[1].set(xlabel="Retained time points", ylabel="Raw global TDMI", title="Selected eligible scans")
        for ax in axes: ax.legend()
    elif kind == "prediction_distribution":
        df = bundle.table(recipe["table"])
        splits = [s for s in ("train", "test", "projection_only") if s in df.split.unique()]
        fig, axes = grid(len(splits), columns=3, width=4.5, height=3.5)
        for ax, split in zip(axes, splits):
            panel = df.loc[df.split == split]
            groups = [g for g in GROUPS if g in panel.group.unique()]
            for gi, group in enumerate(groups):
                vals = panel.loc[panel.group == group, recipe["value"]].to_numpy(float)
                vals = vals[np.isfinite(vals)]
                if len(vals) >= 2 and np.ptp(vals) > 0:
                    orientation = ({"orientation": "horizontal"} if "orientation" in signature(ax.violinplot).parameters
                                   else {"vert": False})  # Support Matplotlib 3.9 as well.
                    parts = ax.violinplot(vals, positions=[gi], showextrema=False, **orientation)
                    for body in parts["bodies"]:
                        body.set(facecolor=color(group, gi), edgecolor=".3", alpha=.8)
                    q1, median, q3 = np.percentile(vals, [25, 50, 75])
                    ax.hlines(gi, q1, q3, color=".25", linewidth=4)
                    ax.scatter(median, gi, color="white", s=12, zorder=4)
                elif len(vals):
                    ax.scatter(vals, np.full(len(vals), gi), color=color(group, gi), s=20)
            ax.set_yticks(range(len(groups)), [f"{g} (n={np.isfinite(panel.loc[panel.group == g, recipe['value']]).sum()})" for g in groups])
            ax.invert_yaxis()
            ax.set_xlabel("Hierarchy correlation" if recipe["value"] == "correlation" else "LDA projection")
            clean_axes(ax)
            ax.set_title({"train": "Training (in-sample)", "test": "Held-out test", "projection_only": "MCI projection only"}[split])
    elif kind == "confusion":
        values = np.asarray(bundle.array("confusion"))
        fig, axes = grid(1, width=5, height=4.5, columns=1)
        matrix(axes[0], values, "Held-out CN/AD predictions")
        axes[0].set_xticks([0, 1], ["CN", "AD"]); axes[0].set_yticks([0, 1], ["CN", "AD"])
        axes[0].set(xlabel="Predicted", ylabel="Observed")
        for i in range(2):
            for j in range(2): axes[0].text(j, i, str(values[i, j]), ha="center", va="center", color="white")
    elif kind == "longitudinal":
        df = bundle.table(recipe["table"])
        sub = df.loc[(df.axis == recipe["axis"]) & (df.measure == recipe["measure"])]
        networks = sorted(bundle.table("visits").network_id.unique())
        fig, axes = grid(len(networks), columns=2, width=5, height=3.2)
        for ax, net in zip(axes, networks):
            panel = sub.loc[sub.network_id == net]
            if panel.empty or not np.isfinite(panel["mean"]).any():
                no_data(ax, "No valid dates/ages/values for this panel")
            else:
                for gi, group in enumerate(GROUPS):
                    data = panel.loc[panel.group == group].sort_values("x")
                    ax.plot(data.x, data["mean"], ".-", label=group, color=color(group, gi), linewidth=1)
                ax.legend(fontsize=8)
            ax.set_title(NETWORKS[net]); ax.set_ylabel(recipe["measure"])
            ax.set_xlabel(recipe["axis"].replace("_", " "), fontsize=8)
    else:
        raise ValueError(f"Unknown saved figure kind: {kind}")
    if bundle.manifest.get("exploratory"):
        fig.text(.995, .002, "Exploratory analysis", ha="right", va="bottom", fontsize=7, color=".4")
    return fig


def _bundles(path):
    path = Path(path)
    if (path / "manifest.json").exists():
        return [path]
    if (path / "run.json").exists():
        return sorted(p for p in path.iterdir() if p.is_dir() and (p / "manifest.json").is_file())
    raise FileNotFoundError(f"Not a result bundle or run directory: {path}")


def plot(figure, results):
    output, found = [], False
    for path in _bundles(results):
        if not verify(path):
            raise ValueError(f"Result bundle is incomplete or corrupt: {path}")
        bundle = Bundle(path)
        if bundle.manifest["status"] == "insufficient_data":
            if figure == "all" or figure in bundle.manifest["declared_figures"]:
                found = True
                output.append({"analysis": bundle.manifest["analysis_id"], "status": "insufficient_data", "reason": bundle.manifest["reason"]})
            continue
        recipes = json.loads((path / "figures.json").read_text())
        # Derive a presentation-only recipe for old bundles, without resealing them.
        if bundle.manifest["analysis_id"] == "coefficients" and not any(
                r["id"] == "coefficients.combined" for r in recipes):
            recipes.append({"id": "coefficients.combined", "kind": "coefficients", "combined": True})
        for recipe in recipes:
            if figure != "all" and recipe["id"] != figure:
                continue
            found = True
            fig = render(bundle, recipe)
            config = bundle.manifest["configuration"]
            cohort, transform = config["cohort"], config["transform"]
            normalisation = bundle.details.get("tdmi_normalisation", transform["tdmi"])
            cohort_name = "all-visits" if bundle.manifest["analysis_id"].startswith("longitudinal") else cohort["visit"]
            timestamp = stamp()
            stem = (f"{recipe['id']}_{cohort_name}_p{cohort['parcels']}_{cohort['estimator']}_"
                    f"tp{cohort['time_min']}-{cohort['time_max']}_tdmi-{normalisation}_{timestamp}_{bundle.manifest['run_id'][:8]}")
            stem = re.sub(r"[^A-Za-z0-9_.-]", "-", stem)
            directory = path / "figures"
            directory.mkdir(exist_ok=True)
            files = []
            try:
                for extension in ("svg", "png", "pdf"):
                    fd, temporary = tempfile.mkstemp(prefix=".render-", suffix="." + extension, dir=directory)
                    os.close(fd)
                    target = directory / f"{stem}.{extension}"
                    try:
                        with plt.rc_context(STYLE):
                            fig.savefig(temporary, format=extension, dpi=300, bbox_inches="tight")
                        os.replace(temporary, target)
                    finally:
                        if os.path.exists(temporary): os.unlink(temporary)
                    files.append(str(target))
            finally:
                plt.close(fig)
            record = {"style_version": STYLE_VERSION, "renderer_sha256": file_hash(Path(__file__)),
                      "figure_id": recipe["id"], "created_at": utc_now(), "recipe": recipe,
                      "source_manifest_sha256": file_hash(path / "manifest.json"),
                      "files": {str(Path(f).name): file_hash(f) for f in files}}
            atomic_json(directory / f"{stem}.json", record)
            output.append({"figure_id": recipe["id"], "status": "complete", "files": files})
    if not found:
        raise ValueError(f"Figure {figure} is not available in {results}; run its analysis first")
    return output
