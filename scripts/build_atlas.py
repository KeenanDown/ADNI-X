"""Construct a full mapping without inventing subcortical parcel identities."""
import argparse
from pathlib import Path
import pandas as pd
from adni_x.inputs import read_atlas

TOKENS = {"Vis": 1, "SomMot": 2, "DorsAttn": 3, "SalVentAttn": 4,
          "Limbic": 5, "Cont": 6, "Default": 7}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cortex", type=Path, required=True, help="Schaefer centroid CSV with ROI Label and ROI Name")
    parser.add_argument("--subcortex", type=Path, help="Explicit parcel_index,parcel_id,name CSV; indices follow cortex")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cortex = pd.read_csv(args.cortex)
    cortex = cortex.sort_values("ROI Label")
    labels = cortex["ROI Label"].to_numpy()
    if list(labels) != list(range(1, len(cortex) + 1)):
        raise ValueError("Cortical ROI Label must be the contiguous one-based tensor order")
    networks = []
    for name in cortex["ROI Name"]:
        matches = [value for token, value in TOKENS.items() if token in str(name).split("_")]
        if len(matches) != 1:
            raise ValueError(f"Cannot unambiguously map cortical network: {name}")
        networks.append(matches[0])
    atlas = pd.DataFrame({"parcel_index": labels - 1, "parcel_id": cortex["ROI Name"].to_numpy(),
                          "name": cortex["ROI Name"].to_numpy(), "network_id": networks})
    if args.subcortex:
        sub = pd.read_csv(args.subcortex)
        if not {"parcel_index", "parcel_id", "name"} <= set(sub):
            raise ValueError("Subcortex CSV requires explicit parcel_index,parcel_id,name")
        sub["network_id"] = 8
        atlas = pd.concat([atlas, sub], ignore_index=True)
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite atlas: {args.output}")
    # Validate before publishing; the function checks order, length and identifiers.
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        candidate = Path(tmp) / "atlas.csv"
        atlas.to_csv(candidate, index=False)
        read_atlas(candidate, len(atlas))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    atlas.to_csv(args.output, index=False)
    print(f"Wrote {len(atlas)} parcels to {args.output}. Verify order against S8b before analysis.")


if __name__ == "__main__":
    main()
