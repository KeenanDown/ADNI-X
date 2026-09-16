from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
import json
import yaml


@dataclass(frozen=True)
class Inputs:
    tensors: str = ""
    series: str = ""
    clinical: str = ""
    atlas: str = ""
    freesurfer: str | None = None
    resolutions: str | None = None
    tensor_variable: str = "PhiIDTensor"
    series_variable: str | None = None
    # Explicit axis order after decoding MATLAB's storage representation.
    tensor_axes: tuple[int, int, int] = (0, 1, 2)
    tensor_glob: str | None = None
    series_glob: str | None = None


@dataclass(frozen=True)
class Cohort:
    parcels: int = 232
    estimator: str = "mmi"
    time_min: int = 150
    time_max: int = 250
    visit: str = "latest"


@dataclass(frozen=True)
class Transform:
    tdmi: str = "none"
    epsilon: float = 1e-12
    confounds: tuple[str, ...] = ("sex", "education", "visit_age", "ICV")
    missing: str = "complete"
    adjustment: str = "mean"
    age: str = "visit"
    exclude_diagonal: bool = True
    regional_axis: int = 1
    network_edges: str = "both"
    block_features: int = 4096


@dataclass(frozen=True)
class Statistics:
    seed: int = 42
    resamples: int = 10000
    permutations: int = 10000
    alpha: float = .05
    welch: bool = True
    correction: str = "fdr_bh"
    test_fraction: float = .25


@dataclass(frozen=True)
class Config:
    inputs: Inputs = field(default_factory=Inputs)
    cohort: Cohort = field(default_factory=Cohort)
    transform: Transform = field(default_factory=Transform)
    statistics: Statistics = field(default_factory=Statistics)
    cache: str = "cache"
    results: str = "results"

    def __post_init__(self):
        c, t, s = self.cohort, self.transform, self.statistics
        if c.parcels < 2 or c.estimator not in {"mmi", "ccs"}:
            raise ValueError("parcels must be >=2; estimator must be mmi or ccs")
        if not 0 < c.time_min <= c.time_max:
            raise ValueError("Require 0 < time_min <= time_max")
        if c.visit not in {"latest", "baseline"} and not (
            c.visit.startswith("M") and c.visit[1:].isdigit()
        ):
            raise ValueError("visit must be latest, baseline, or M###")
        if sorted(self.inputs.tensor_axes) != [0, 1, 2]:
            raise ValueError("tensor_axes must be a permutation of [0,1,2]")
        if t.tdmi not in {"none", "tdmi"} or t.epsilon <= 0:
            raise ValueError("tdmi must be none or tdmi; epsilon must be positive")
        if len(set(t.confounds)) != len(t.confounds) or not set(t.confounds) <= {
            "sex", "education", "visit_age", "ICV"
        }:
            raise ValueError("Unknown or duplicate confounds")
        if t.missing not in {"complete", "median"}:
            raise ValueError("missing must be complete or median")
        if t.adjustment not in {"mean", "intercept", "residuals"}:
            raise ValueError("adjustment must be mean, intercept, or residuals")
        if t.age not in {"visit", "baseline"} or t.regional_axis not in {0, 1}:
            raise ValueError("Invalid age or regional_axis")
        if t.network_edges not in {"both", "upper"} or t.block_features < 1:
            raise ValueError("Invalid network_edges or block_features")
        if s.correction not in {"fdr_bh", "none"}:
            raise ValueError("correction must be fdr_bh or none")
        if s.resamples < 2 or s.permutations < 2 or not 0 < s.alpha < 1:
            raise ValueError("Require >=2 resamples/permutations and 0<alpha<1")
        if s.seed < 0 or not 0 < s.test_fraction < 1:
            raise ValueError("Invalid seed or test_fraction")

    def to_dict(self):
        return json.loads(json.dumps(asdict(self)))

    @classmethod
    def load(cls, path: str | Path):
        path = Path(path).resolve()
        raw = yaml.safe_load(path.read_text()) or {}
        allowed = {"inputs", "cohort", "transform", "statistics", "cache", "results"}
        if set(raw) - allowed:
            raise ValueError(f"Unknown configuration keys: {set(raw) - allowed}")
        inp = dict(raw.get("inputs", {}))
        for key in ("tensors", "series", "clinical", "atlas", "freesurfer", "resolutions"):
            if inp.get(key):
                inp[key] = str((path.parent / Path(inp[key]).expanduser()).resolve())
        if "tensor_axes" in inp:
            inp["tensor_axes"] = tuple(inp["tensor_axes"])
        trans = dict(raw.get("transform", {}))
        if "confounds" in trans:
            trans["confounds"] = tuple(trans["confounds"])
        return cls(
            inputs=Inputs(**inp), cohort=Cohort(**raw.get("cohort", {})),
            transform=Transform(**trans), statistics=Statistics(**raw.get("statistics", {})),
            cache=str((path.parent / raw.get("cache", "cache")).resolve()),
            results=str((path.parent / raw.get("results", "results")).resolve()),
        )


def ensure_config(config):
    return config if isinstance(config, Config) else Config.load(config)

