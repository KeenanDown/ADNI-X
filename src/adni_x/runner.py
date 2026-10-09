from __future__ import annotations

from pathlib import Path
import json
import logging
import traceback
import fcntl
import numpy as np
from scipy.io import savemat

from .catalogue import ANALYSES
from .config import ensure_config
from .preparation import prepare
from .provenance import atomic_json, digest, environment, seal, transaction, utc_now, verify
from .results import Result, InsufficientData

log = logging.getLogger(__name__)


def publish_status(root, report, target):
    """Independent HPC targets may finish concurrently; never lose another target's status."""
    root = Path(root)
    path = root / "run.json"
    with (root / ".run.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            current = json.loads(path.read_text()) if path.exists() else {"targets": {}}
            current.update({k: v for k, v in report.items() if k != "targets"})
            current["targets"][target] = report["targets"][target]
            atomic_json(path, current)
            return current
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)


def run(analysis, config, *, render=True, include_sensitivity=False):
    config = ensure_config(config)
    if analysis != "all" and analysis not in ANALYSES:
        raise ValueError(f"Unknown analysis {analysis}; use adni-x list")
    prepared = prepare(config)
    prep = prepared.manifest
    run_id = digest({"preparation": prep["cache_id"], "statistics": config.to_dict()["statistics"]})
    root = Path(config.results) / run_id[:16]
    root.mkdir(parents=True, exist_ok=True)
    chosen = [a for a in ANALYSES.values() if analysis == a.id or (analysis == "all" and (not a.opt_in or include_sensitivity))]
    report_path = root / "run.json"
    report = json.loads(report_path.read_text()) if report_path.exists() else {"run_id": run_id, "targets": {}}
    report.update(configuration=config.to_dict(), preparation=str(prepared.path), updated_at=utc_now())
    if analysis == "all" and not include_sensitivity:
        report["opt_in_not_requested"] = [a.id for a in ANALYSES.values() if a.opt_in]
    for spec in chosen:
        cohort = prepared.for_group(spec.grouping)
        destination = root / spec.id
        log.info("Computing %s", spec.id)
        try:
            if destination.exists():
                if not verify(destination):
                    raise ValueError(f"Corrupt or incomplete result bundle: {destination}; move it aside and retry")
                manifest = json.loads((destination / "manifest.json").read_text())
                status = manifest["status"]
            else:
                reason, status = None, "complete"
                try:
                    result = spec.function(prepared)
                except InsufficientData as e:
                    result, reason, status = Result(), str(e), "insufficient_data"
                if status == "complete":
                    actual = {f["id"] for f in result.figures}
                    if actual != set(spec.figures):
                        raise ValueError(f"Figure catalogue mismatch for {spec.id}: {actual ^ set(spec.figures)}")
                with transaction(destination) as out:
                    for name, table in result.tables.items():
                        table.to_csv(out / f"{name}.csv", index=False)
                    for name, array in result.arrays.items():
                        np.save(out / f"{name}.npy", array, allow_pickle=False)
                    for name, content in result.matlab.items():
                        savemat(out / f"{name}.mat", content, do_compression=True)
                    atomic_json(out / "figures.json", result.figures)
                    atomic_json(out / "details.json", result.metadata)
                    # Each bundle is self-contained for provenance and plotting.
                    selected = cohort.table("selected")
                    selected.to_csv(out / "cohort_selected.csv", index=False)
                    cohort.table("subjects").to_csv(out / "cohort_adjusted.csv", index=False)
                    if spec.dependency.startswith("all eligible"):
                        cohort.table("visits").to_csv(out / "cohort_visits.csv", index=False)
                    atomic_json(out / "preparation_audit.json", json.loads((prepared.path / "audit.json").read_text()))
                    seal(out, {"kind": "analysis", "analysis_id": spec.id, "grouping": spec.grouping,
                               "groups": cohort.groups, "reference_group": cohort.groups[0], "run_id": run_id, "status": status,
                               "reason": reason, "created_at": utc_now(), "configuration": config.to_dict(),
                               "preparation_id": prep["cache_id"], "sources": prep["sources"],
                               "environment": environment(), "exploratory": spec.exploratory,
                               "notebook_cells_zero_based": spec.cells, "declared_figures": list(spec.figures),
                               "tables": list(result.tables), "arrays": list(result.arrays)})
            report["targets"][spec.id] = {"status": status, "path": str(destination)}
            if render and status == "complete":
                from .plotting import plot
                try:
                    figures = plot("all", destination)
                    report["targets"][spec.id]["figures"] = figures
                except Exception as e:
                    # Completed numerical output survives a rendering failure.
                    report["targets"][spec.id]["plot_status"] = "failed"
                    report["targets"][spec.id]["plot_error"] = str(e)
                    log.exception("Rendering failed for %s", spec.id)
        except Exception as e:
            report["targets"][spec.id] = {"status": "failed", "error": str(e)}
            atomic_json(root / "failures" / f"{spec.id}.json", {"time": utc_now(), "error": str(e), "traceback": traceback.format_exc()})
            log.exception("Analysis failed: %s", spec.id)
        report = publish_status(root, report, spec.id)
    return {"path": str(root), **report}
