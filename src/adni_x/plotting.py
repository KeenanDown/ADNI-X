"""Rendering depends exclusively on immutable, saved numerical result bundles."""
from __future__ import annotations

from pathlib import Path
import json
import math
import os
import re
import tempfile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

from .constants import GROUPS, GROUP_COLORS, NETWORK_COLORS, NETWORKS, ATOMS
from .provenance import atomic_json, file_hash, stamp, utc_now, verify

plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none", "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False, "font.size": 10})


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


def matrix(ax, array, title, *, diverging=False, limits=None):
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
    im = ax.imshow(arr, cmap="RdBu_r" if diverging else "magma", vmin=limits[0], vmax=limits[1], interpolation="nearest")
    ax.set_xticks([]); ax.set_yticks([])
    ax.figure.colorbar(im, ax=ax, shrink=.7)


def bars(ax, networks):
    cmap = ListedColormap([NETWORK_COLORS[int(n)] for n in networks])
    n = len(networks)
    bottom = ax.inset_axes([0, -.035, 1, .022])
    bottom.imshow(np.arange(n)[None, :], cmap=cmap, aspect="auto", interpolation="nearest")
    bottom.axis("off")


def network_legend(fig, networks):
    ids = sorted(set(map(int, networks)))
    fig.legend(handles=[Patch(color=NETWORK_COLORS[n], label=NETWORKS[n]) for n in ids],
               loc="outside lower center", ncol=min(4, len(ids)), fontsize=8, frameon=False)


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


def render(bundle, recipe):
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
        fig, axes = grid(4, columns=2, width=5, height=4.3)
        finite = means[np.isfinite(means)]
        limits = (finite.min(), finite.max()) if len(finite) else None
        counts = bundle.table("group_counts").set_index("group")["n"].to_dict()
        for gi, group in enumerate(GROUPS):
            matrix(axes[gi], means[gi][np.ix_(order, order)], f"{group} (n={counts.get(group, 0)})", limits=limits)
            bars(axes[gi], atlas.network_id.to_numpy()[order])
        matrix(axes[3], difference[np.ix_(order, order)], "AD − CN, significant edges", diverging=True)
        bars(axes[3], atlas.network_id.to_numpy()[order])
        network_legend(fig, atlas.network_id)
    elif kind == "coefficients":
        atom = recipe["atom"]
        values = bundle.array("coefficients")[:, :, :, atom]
        names = bundle.details["names"]
        fig, axes = grid(len(values), width=4, height=3.7)
        atlas = bundle.table("atlas")
        order = np.argsort(atlas.network_id, kind="stable")
        for ax, values, name in zip(axes, values, names):
            a = np.array(values[np.ix_(order, order)])
            if bundle.manifest["configuration"]["transform"]["exclude_diagonal"]:
                np.fill_diagonal(a, np.nan)
            matrix(ax, a, name, diverging=True)
            bars(ax, atlas.network_id.to_numpy()[order])
        network_legend(fig, atlas.network_id)
    elif kind == "effects_grid":
        df = bundle.table(recipe["table"])
        df = df.loc[(df.group_a == recipe["group_a"]) & (df.group_b == recipe["group_b"]) & df.atom.isin(recipe["atoms"])]
        nets = sorted(df.network_id.unique())
        pivot = df.pivot(index="atom", columns="network_id", values=recipe["effect"]).reindex(index=recipe["atoms"], columns=nets)
        significant = df.pivot(index="atom", columns="network_id", values="significant").reindex_like(pivot).fillna(False).to_numpy(bool)
        fig, axes = grid(1, width=max(7, len(nets)*1.1), height=max(3, len(pivot)*.42), columns=1)
        ax = axes[0]
        matrix(ax, np.where(significant, pivot.to_numpy(), np.nan), f"{recipe['group_b']} − {recipe['group_a']} · {recipe['effect']}", diverging=True)
        ax.set_xticks(range(len(nets)), [NETWORKS[n] for n in nets], rotation=35, ha="right")
        ax.set_yticks(range(len(pivot)), pivot.index)
        for r, c in zip(*np.where(significant)):
            ax.text(c, r, f"{pivot.iloc[r,c]:.2f}", ha="center", va="center", fontsize=8)
    elif kind == "significance_strip":
        df = bundle.table("effects")
        atom = df.loc[df.atom == recipe["atom"]]
        anchor = atom.loc[(atom.group_a == "CN") & (atom.group_b == "AD")].sort_values("difference", ascending=False, kind="stable")
        sub = atom.loc[(atom.group_a == recipe["group_a"]) & (atom.group_b == recipe["group_b"])].set_index("parcel_index").reindex(anchor.parcel_index)
        fig, axes = grid(1, width=12, height=4, columns=1)
        ax = axes[0]
        x = np.arange(len(sub))
        ax.axhline(0, color=".7", linewidth=1)
        ax.plot(x, sub.difference, color=".75", linewidth=.8)
        for sig, marker in ((False, "."), (True, "o")):
            mask = sub.significant.fillna(False).to_numpy(bool) == sig
            ax.scatter(x[mask], sub.difference.to_numpy()[mask], c=[NETWORK_COLORS[int(n)] for n in sub.network_id.to_numpy()[mask]],
                       marker=marker, s=22 if sig else 8, alpha=1 if sig else .35)
        ax.set(xlabel="Parcels ordered by AD − CN", ylabel=f"{recipe['atom']} regional-sum difference",
               title=f"{recipe['group_b']} − {recipe['group_a']} · larger points pass {threshold_name}")
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
        fig, axes = grid(1 + len(networks), columns=3, width=4.5, height=3.5)
        for gi, group in enumerate(bundle.details["groups"]):
            sub = df.loc[df.group == group].set_index("parcel_index").reindex(order)
            axes[0].plot(range(len(sub)), sub[measure], label=group, color=color(group, gi), linewidth=1)
            for ax, network in zip(axes[1:], networks):
                panel = sub.loc[sub.network_id == network]
                ax.scatter(range(len(panel)), panel[measure], s=16, color=color(group, gi), label=group, alpha=.7)
                ax.set_title(NETWORKS[network], color=NETWORK_COLORS[network])
        axes[0].set_title("All parcels")
        for ax in axes:
            ax.set(xlabel=f"Parcel rank ({ordering})", ylabel=measure.replace("_", " "))
        axes[0].legend(fontsize=8)
        network_legend(fig, df.network_id)
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
        taxa = list(bundle.details["taxonomy"])
        nets = sorted(df.network_id.unique())
        fig, axes = plt.subplots(len(taxa), len(nets), figsize=(max(8, 2.7*len(nets)), 2.9*len(taxa)), squeeze=False, layout="constrained")
        for ti, taxon in enumerate(taxa):
            for ni, net in enumerate(nets):
                ax = axes[ti, ni]
                sub = df.loc[(df.taxon == taxon) & (df.network_id == net)]
                distributions(ax, sub, "value", bundle.details["groups"])
                test = tests.loc[(tests.taxon == taxon) & (tests.network_id == net)].iloc[0]
                label = "q" if correction == "fdr_bh" else "p"
                ax.set_title(f"{taxon.replace('_', ' ')}\n{NETWORKS[net]} · {label}={test.q:.3g}", fontsize=9)
                ax.tick_params(labelsize=8); ax.set_ylabel("Intranetwork value")
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
                ax.scatter(sub.rtr, sub.sts, label=group, color=color(group, gi), s=20, alpha=.65)
        ax.set(xlabel="CN-reference redundancy (standard deviations)", ylabel="CN-reference synergy (standard deviations)")
        ax.legend(fontsize=8)
    elif kind == "pca":
        group = recipe["group"]
        df = bundle.table("regions").query("group == @group")
        pcs = bundle.table("pca_scores").query("group == @group")
        model = next(m for m in bundle.details["pca"] if m["group"] == group)
        fig, axes = grid(2, width=5, height=4.5, columns=2)
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
            axes[1].set(title="Projection onto PC1", xlabel="Redundancy coordinate", ylabel="Synergy coordinate")
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
        ax.set_xticks(range(len(sub)), [NETWORKS[n] for n in sub.network_id], rotation=30, ha="right")
        ax.set(ylabel=f"{recipe['atom']} regional-sum change", title=f"{recipe['group_b']} − {recipe['group_a']} · exploratory regional shuffle")
    elif kind == "group_distribution":
        df = bundle.table(recipe["table"])
        fig, axes = grid(1, width=6, height=5, columns=1)
        distributions(axes[0], df, recipe["value"])
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
        fig, axes = grid(len(splits), columns=3, width=4.5, height=4)
        for ax, split in zip(axes, splits):
            distributions(ax, df.loc[df.split == split], recipe["value"])
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
        fig, axes = grid(len(networks), width=4.5, height=3.5)
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
    fig.suptitle(recipe["id"].replace(".", " · "), fontsize=12)
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
                for extension in ("svg", "png"):
                    fd, temporary = tempfile.mkstemp(prefix=".render-", suffix="." + extension, dir=directory)
                    os.close(fd)
                    target = directory / f"{stem}.{extension}"
                    try:
                        fig.savefig(temporary, format=extension, dpi=300, bbox_inches="tight")
                        os.replace(temporary, target)
                    finally:
                        if os.path.exists(temporary): os.unlink(temporary)
                    files.append(str(target))
            finally:
                plt.close(fig)
            record = {"figure_id": recipe["id"], "created_at": utc_now(), "recipe": recipe,
                      "source_manifest_sha256": file_hash(path / "manifest.json"),
                      "files": {str(Path(f).name): file_hash(f) for f in files}}
            atomic_json(directory / f"{stem}.json", record)
            output.append({"figure_id": recipe["id"], "status": "complete", "files": files})
    if not found:
        raise ValueError(f"Figure {figure} is not available in {results}; run its analysis first")
    return output
