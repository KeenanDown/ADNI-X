from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np
import pandas as pd


class InsufficientData(ValueError):
    """An expected, explicitly reported data limitation, not a programming error."""


@dataclass
class Result:
    tables: dict[str, pd.DataFrame] = field(default_factory=dict)
    arrays: dict[str, np.ndarray] = field(default_factory=dict)
    figures: list[dict] = field(default_factory=list)
    metadata: dict = field(default_factory=dict)
    matlab: dict[str, dict] = field(default_factory=dict)

    def figure(self, name, kind, **parameters):
        self.figures.append({"id": name, "kind": kind, **parameters})

