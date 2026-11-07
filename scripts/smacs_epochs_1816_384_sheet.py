"""Two-epoch cutouts of SMACS ``f200w_1816`` / ``f200w_384`` (F150W, F444W; 2736 vs 6882).

S3 byte-range reads of the level-3 ``_i2d.fits``; one row per target, band and epoch; north-up
contact sheet with brightness shared per band. Visual check only (no stored measurement).

    python scripts/smacs_epochs_1816_384_sheet.py <scratch dir> \\
        --out results/smacs_jwst3/epochs_1816_384   # reads summary.json there
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from astropy.table import Table, vstack

from jwst_anomaly import cutouts, viz

sys.path.insert(0, str(Path(__file__).resolve().parent))
from smacs_epochs_1816_384 import BANDS, E1, E2  # noqa: E402

# S3 directories of the two level-3 products whose CAT names are E1 / E2.
EPOCHS = {"2022": ("jw02736/L3/t/o001/", E1), "2026": ("jw06882/L3/t/o057/", E2)}
SIZE_ARCSEC = 4.0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("work", type=Path, help="scratch directory for the cutout FITS files")
    ap.add_argument("--out", type=Path, required=True, help="dir holding summary.json")
    args = ap.parse_args(argv)
    # Target positions as measured by smacs_epochs_1816_384.py (summary.json).
    measured = json.loads((args.out / "summary.json").read_text())["targets"]
    targets = Table(
        {
            "source_uid": list(measured),
            "ra": [r["ra"] for r in measured.values()],
            "dec": [r["dec"] for r in measured.values()],
        }
    )
    rows = []
    for epoch, (s3dir, name) in EPOCHS.items():
        for band in BANDS:
            i2d = name.format(b=band).replace("_cat.ecsv", "_i2d.fits")
            t = cutouts.make_cutouts(
                f"s3://stpubdata/jwst/public/{s3dir}{i2d}",
                targets,
                size_arcsec=SIZE_ARCSEC,
                out_dir=args.work,
                storage_options={"anon": True},
            )
            t["epoch"] = epoch
            t.meta = {}
            rows.append(t)
    allc = vstack(rows)
    # One sheet cell per (uid, band): prefix the epoch so the two epochs get separate cells.
    pairs = zip(allc["epoch"], allc["source_uid"], strict=True)
    allc["source_uid"] = [f"{e} {u}" for e, u in pairs]
    viz.contact_sheet(allc, args.out / "contact_sheet.png", north_up=True, scale="band")
    allc["source_uid", "epoch", "band", "quality_flag", "frac_nan", "wht_rel"].write(
        args.out / "cutouts.ecsv", overwrite=True
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
