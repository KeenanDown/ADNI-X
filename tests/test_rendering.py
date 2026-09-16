import json
from pathlib import Path
import matplotlib.pyplot as plt
import pytest

from adni_x.catalogue import FIGURES
from adni_x.plotting import Bundle, render
from adni_x.runner import run


@pytest.fixture(scope="session")
def bundles(prepared):
    report = run("all", prepared.config, render=False, include_sensitivity=True)
    assert all(item["status"] == "complete" for item in report["targets"].values())
    return Path(report["path"])


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
