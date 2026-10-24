"""MARS MACS1149 map vs the Schuldt et al. 2024 new multiple images (held-out stage).

Frozen plan: ``results/lens_consistency/mars_m1149_heldout/PLAN.md``. Rules 1-5 of the
``mars_hff_heldout`` stage, unchanged, with the Schuldt+2024 Table 2 list in place of a CANUCS
list. No flag, flux or extra-image reading.

    python scripts/mars_m1149_heldout.py --zip DIR

``DIR`` receives the MARS zip (md5 checked); it is deleted after the run unless ``--keep``.
"""

from __future__ import annotations

import argparse
import json
import zipfile
from pathlib import Path

import mars_hff_heldout as mo
import mars_hff_inputs as mh
import requests
from m0417_canucsref_posthoc import _clean

OUT = Path("results/lens_consistency/mars_m1149_heldout")
FIELD = "MACS J1149.5-2223"  # MARS catalogue field name
IMAGES = OUT / "schuldt24_new_images.ecsv"
K_TRAIN = 0.742  # MACS1149 training gate K rms, D-034 addendum ao


def run(zdir: Path, out: Path = OUT) -> dict:
    meta = mh.fetch_zip(requests.Session(), zdir)
    out.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zdir / mh.ZIP) as zf:
        cat = mh.parse_catalogue(zf.read(mh.CATALOGUE).decode())
        block = mo.field_block(zf, FIELD, cat, out, images=(IMAGES, K_TRAIN))
    summary = {"plan": str(out / "PLAN.md"), "zip": meta, "images": str(IMAGES), "field": block}
    (out / "summary.json").write_text(json.dumps(_clean(summary), indent=1, allow_nan=False) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--zip", type=Path, required=True, help="directory for the zip")
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--keep", action="store_true", help="keep the downloaded zip")
    a = ap.parse_args(argv)
    try:
        s = run(a.zip, a.out)
    finally:
        if not a.keep:
            (a.zip / mh.ZIP).unlink(missing_ok=True)
    print(json.dumps(_clean(s["field"]), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
