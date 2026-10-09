#!/usr/bin/env python3
"""Check release contents and install the wheel into a temporary target offline."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import zipfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--sdist", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    wheel, sdist = args.wheel.resolve(), args.sdist.resolve()
    with zipfile.ZipFile(wheel) as archive:
        assert "adni_x/inference.py" in archive.namelist()
        metadata = next(n for n in archive.namelist() if n.endswith(".dist-info/METADATA"))
        assert "Version: 0.2.0" in archive.read(metadata).decode()
    with tarfile.open(sdist) as archive:
        members = archive.getnames()
        for required in ("docs/hpc-rerun-0.2.0.md", "notebooks/ADNI_X.ipynb", "tests/test_rerun.py",
                         "scripts/qmul_slurm.sh", "scripts/run_hpc.sh", "requirements.lock"):
            assert any(n.endswith("/"+required) for n in members), required
        assert not any(n.endswith(".mat") or "/.venv/" in n or "/reports/" in n for n in members)
        for script in ("scripts/qmul_slurm.sh", "scripts/run_hpc.sh"):
            name = next(n for n in members if n.endswith("/"+script))
            assert b"\r\n" not in archive.extractfile(name).read()
    with tempfile.TemporaryDirectory(prefix="adnix-wheel-") as target:
        subprocess.run([sys.executable, "-m", "pip", "install", "--no-index", "--no-deps", "--no-compile",
                        "--target", target, str(wheel)], check=True, stdout=subprocess.PIPE, text=True)
        check = '''
import contextlib, io, json, sys
from pathlib import Path
sys.path.insert(0, sys.argv[1])
import adni_x
from adni_x.catalogue import catalogue
from adni_x.cli import main
from adni_x.inference import hc3_coefficients
import numpy as np
assert adni_x.__version__ == "0.2.0"
assert Path(adni_x.__file__).is_relative_to(Path(sys.argv[1]))
stream = io.StringIO()
with contextlib.redirect_stdout(stream):
    try: main(["--version"])
    except SystemExit as e: assert e.code == 0
assert stream.getvalue().strip() == "ADNI-X 0.2.0"
targets = {t["id"]: t for t in catalogue()}
assert len(targets) == 22
assert not targets["taxa-principal-moca"]["opt_in"]
assert targets["taxa-six-mmse"]["opt_in"]
X = np.column_stack([np.ones(8), np.arange(8)])
Y = (.4+3*np.arange(8)+np.array([.1, -.1]*4))[:, None]
beta = np.linalg.lstsq(X, Y, rcond=None)[0]
assert (hc3_coefficients(X, Y, beta)["status"] == "ok").all()
print(json.dumps({"version": adni_x.__version__, "cli_version": stream.getvalue().strip(),
                  "isolated_import": True, "catalogue_targets": len(targets), "hc3_smoke": "passed"}))
'''
        result = subprocess.run([sys.executable, "-I", "-c", check, target], check=True,
                                stdout=subprocess.PIPE, text=True)
        report = json.loads(result.stdout)
    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (wheel, sdist)}
    report["sha256"] = hashes
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+"\n")
    (wheel.parent/"SHA256SUMS-0.2.0.txt").write_text("".join(f"{digest}  {name}\n" for name, digest in hashes.items()))
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
