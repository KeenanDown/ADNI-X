import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pytest

from adni_x.catalogue import FIGURES
from adni_x.plotting import Bundle, render, plot, STYLE_VERSION
from adni_x.provenance import file_hash, verify
from adni_x.runner import run


@pytest.mark.parametrize("figure", list(FIGURES))
def test_every_registered_renderer_draws_from_saved_results(bundles, figure):
    path = bundles / FIGURES[figure]
    recipes = json.loads((path / "figures.json").read_text())
    recipe = next(r for r in recipes if r["id"] == figure)
    fig = render(Bundle(path), recipe)
    try:
        fig.canvas.draw()
        assert fig.axes
        assert fig.get_size_inches().prod() > 0
    finally:
        plt.close(fig)


def test_effect_grid_preserves_saved_values_and_significance(bundles):
    bundle = Bundle(bundles / "network-effects")
    recipe = json.loads((bundle.path / "figures.json").read_text())[0]
    df = bundle.table(recipe["table"])
    df = df.loc[(df.group_a == recipe["group_a"]) & (df.group_b == recipe["group_b"])
                & df.atom.isin(recipe["atoms"])]
    values = df.pivot(index="atom", columns="network_id", values=recipe["effect"]).reindex(recipe["atoms"])
    flags = df.pivot(index="atom", columns="network_id", values="significant").reindex_like(values)
    before = dict(plt.rcParams)
    fig = render(bundle, recipe)
    try:
        expected = np.where(flags, values, np.nan)
        if np.isfinite(expected).any():
            image = fig.axes[0].images[0]
            np.testing.assert_allclose(image.get_array().filled(np.nan), expected, equal_nan=True)
            assert image.norm.vmin == -image.norm.vmax
        else:
            assert not fig.axes[0].images
        assert dict(plt.rcParams) == before
    finally:
        plt.close(fig)


def test_restyle_old_bundle_exports_editable_files_without_changing_results(bundles):
    path = bundles / "coefficients"
    recipes = json.loads((path / "figures.json").read_text())
    assert "coefficients.combined" not in [r["id"] for r in recipes]
    before = {p.name: file_hash(p) for p in path.iterdir() if p.is_file()}
    result = plot("coefficients.combined", path)[0]
    files = {Path(f).suffix: Path(f) for f in result["files"]}
    assert set(files) == {".svg", ".png", ".pdf"}
    assert "<text" in files[".svg"].read_text()
    assert files[".pdf"].read_bytes().startswith(b"%PDF")
    record = json.loads(files[".svg"].with_suffix(".json").read_text())
    assert record["style_version"] == STYLE_VERSION
    assert len(record["renderer_sha256"]) == 64
    assert before == {p.name: file_hash(p) for p in path.iterdir() if p.is_file()}
    assert verify(path)


def test_notebook_displays_newest_render_and_refreshes_explicit_saved_run(tmp_path, monkeypatch):
    import ast
    import adni_x

    notebook = json.loads((Path(__file__).parents[1] / "notebooks/ADNI_X.ipynb").read_text())
    source = next("".join(c["source"]) for c in notebook["cells"]
                  if "def show_figures(" in "".join(c["source"]))
    tree = ast.parse(source)
    tree.body = [node for node in tree.body if not (
        isinstance(node, ast.ImportFrom) and node.module == "IPython.display")]
    displayed, refreshed = [], []
    namespace = {"display": displayed.append, "Image": lambda filename: filename}
    exec(compile(tree, "notebook-helper", "exec"), namespace)
    namespace["RESULTS_RUN"] = tmp_path
    figures = tmp_path / "coefficients" / "figures"
    figures.mkdir(parents=True)
    for name, figure_id, timestamp in [
        ("old", "coefficients.combined", "2026-01-01"),
        ("new", "coefficients.combined", "2026-02-01"),
        ("single", "coefficients.rtr", "2026-02-01"),
    ]:
        (figures / f"{name}.png").touch()
        (figures / f"{name}.json").write_text(json.dumps({
            "figure_id": figure_id, "created_at": timestamp, "files": {f"{name}.png": "hash"}}))
    namespace["show_figures"]("coefficients")
    assert displayed == [str(figures / "new.png")]
    monkeypatch.setattr(adni_x, "plot", lambda figure, path: refreshed.append((figure, path)))
    namespace["refresh_figures"]("coefficients")
    assert refreshed == [("all", tmp_path / "coefficients")]
