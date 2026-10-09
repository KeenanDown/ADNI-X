"""Independent, provenance-tracked analyses of completed PhiID tensors."""

__version__ = "0.2.0"


def prepare(config):
    from .preparation import prepare as impl
    return impl(config)


def validate(config):
    from .inputs import inventory
    return inventory(config)


def run(analysis, config):
    from .runner import run as impl
    return impl(analysis, config)


def plot(figure, results):
    from .plotting import plot as impl
    return impl(figure, results)

