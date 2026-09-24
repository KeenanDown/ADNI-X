"""Independent, provenance-tracked analyses of completed PhiID tensors."""

__version__ = "0.1.0"


def prepare(config):
    from .preparation import prepare as impl
    return impl(config)


def validate(config):
    from .inputs import inventory
    return inventory(config)


def run(analysis, config, *, render=True, include_sensitivity=False):
    from .runner import run as impl
    return impl(analysis, config, render=render, include_sensitivity=include_sensitivity)


def plot(figure, results):
    from .plotting import plot as impl
    return impl(figure, results)

