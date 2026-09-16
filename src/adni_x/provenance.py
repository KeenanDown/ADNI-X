from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import importlib.metadata
import json
import os
import platform
import shutil
import subprocess
import tempfile
import numpy as np


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def stamp():
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")


def _json_default(value):
    if isinstance(value, (np.integer, np.floating, np.bool_)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(type(value).__name__)


def _clean(value):
    if isinstance(value, dict):
        return {str(k): _clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_clean(v) for v in value]
    if isinstance(value, (np.integer, np.floating, np.bool_)):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def json_text(data):
    return json.dumps(_clean(data), sort_keys=True, indent=2, default=_json_default, allow_nan=False)


def digest(data):
    return hashlib.sha256(json_text(data).encode()).hexdigest()


def file_hash(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def code_hash():
    root = Path(__file__).parent
    return digest({str(p.relative_to(root)): file_hash(p) for p in sorted(root.rglob("*.py"))})


def environment():
    packages = {}
    for name in ("adni-x", "numpy", "scipy", "pandas", "matplotlib", "seaborn", "scikit-learn", "h5py", "PyYAML"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = "source checkout"
    try:
        revision = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=Path(__file__).parent,
            stderr=subprocess.DEVNULL, text=True,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        revision = None
    return {"python": platform.python_version(), "platform": platform.platform(),
            "packages": packages, "git_revision": revision, "implementation_sha256": code_hash()}


def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".json-", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as f:
            f.write(json_text(data) + "\n")
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


@contextmanager
def transaction(destination):
    """Publish a complete directory; an interrupted temporary directory is never a cache hit."""
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix=f".{destination.name}-", dir=destination.parent))
    try:
        yield tmp
        if destination.exists():
            raise FileExistsError(f"Concurrent or existing output: {destination}")
        os.replace(tmp, destination)
    finally:
        if tmp.exists():
            shutil.rmtree(tmp)


def seal(directory, metadata):
    directory = Path(directory)
    files = {str(p.relative_to(directory)): file_hash(p)
             for p in sorted(directory.rglob("*")) if p.is_file() and p.name != "manifest.json"}
    atomic_json(directory / "manifest.json", {**metadata, "files": files})


def verify(directory):
    directory = Path(directory)
    path = directory / "manifest.json"
    if not path.is_file():
        return False
    meta = json.loads(path.read_text())
    return bool(meta.get("files")) and all(
        (directory / name).is_file() and file_hash(directory / name) == h
        for name, h in meta["files"].items()
    )


def rng_for(seed, name):
    words = np.frombuffer(hashlib.sha256(name.encode()).digest()[:16], dtype="<u4")
    return np.random.default_rng(np.random.SeedSequence([seed, *map(int, words)]))
