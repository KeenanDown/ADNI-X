from __future__ import annotations

import argparse
from pathlib import Path
import logging
import sys

from .catalogue import catalogue
from .config import Config
from .provenance import atomic_json, json_text, stamp


def main(argv=None):
    parser = argparse.ArgumentParser(prog="adni-x", description="Reproducible analyses from completed S8b PhiID tensors")
    parser.add_argument("--version", action="version", version="ADNI-X 0.1.0")
    subs = parser.add_subparsers(dest="command", required=True)
    subs.add_parser("list", help="List analyses, independent figures and preparation dependencies")
    for name in ("validate", "prepare", "run"):
        sub = subs.add_parser(name)
        sub.add_argument("--config", required=True, type=Path)
        if name == "validate":
            sub.add_argument("--report", type=Path, help="Directory for inventory, atlas and audit report")
        if name == "run":
            sub.add_argument("analysis", help="Analysis identifier or all")
            sub.add_argument("--no-plots", action="store_true")
            sub.add_argument("--include-sensitivity", action="store_true", help="Include opt-in TDMI-ratio targets in all")
    sub = subs.add_parser("plot")
    sub.add_argument("figure", help="Figure identifier or all")
    sub.add_argument("--results", required=True, type=Path, help="Analysis bundle or run directory; raw data not required")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    try:
        if args.command == "list":
            print(json_text(catalogue()))
        elif args.command == "plot":
            from .plotting import plot
            print(json_text(plot(args.figure, args.results)))
        else:
            config = Config.load(args.config)
            if args.command == "validate":
                from .inputs import inventory
                inv = inventory(config)
                root = args.report or Path(config.results) / "validation" / stamp()
                root.mkdir(parents=True, exist_ok=True)
                inv.scans.to_csv(root / "inventory.csv", index=False)
                inv.atlas.to_csv(root / "atlas.csv", index=False)
                atomic_json(root / "audit.json", inv.summary())
                atomic_json(root / "sources.json", inv.sources)
                print(json_text({"report": str(root), **inv.summary()}))
                if not inv.scans.selected.any():
                    return 2
            elif args.command == "prepare":
                from .preparation import prepare
                dataset = prepare(config)
                print(json_text({"prepared": str(dataset.path), "dimensions": dataset.manifest["dimensions"]}))
            else:
                from .runner import run
                report = run(args.analysis, config, render=not args.no_plots, include_sensitivity=args.include_sensitivity)
                print(json_text(report))
                if any(s["status"] == "failed" or s.get("plot_status") == "failed" for s in report["targets"].values()):
                    return 1
    except (ValueError, FileNotFoundError, KeyError) as e:
        print(f"adni-x: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

