"""`arcs` on MACS J0308.9+2645: an independent test of the S1063 anti excess.

Plan (frozen before the run): ``results/arcs/macs0308-relics/PLAN.md``. Same chain as
``arcs_spt0615.py`` (RELICS BPZ -> zout-like table, 95 % interval as the 68 % one, an ASSUMPTION;
offset removed above 0.1"), with the RELICS GLAFIC v2 (primary) and Zitrin LTM-Gauss v1 (second)
maps on the program 5293 o005 F200W catalogue, and the frozen outcome rule through
``arcs_rxj2129.decide``. Shapes are ``observed``, photo-z ``derived`` (RELICS), predicted axes
``model_prediction``, classes, counts and P values ``derived``.

    python scripts/arcs_macs0308.py --catalog <jw05293-o005_t001 F200W _cat.ecsv> \\
        --bpz <hlsp_relics_hst_acs-wfc3ir_macs0308p26_multi_v1_cat.txt>
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
from arcs_rxj2129 import decide  # noqa: E402
from arcs_spt0615 import bpz_zout, offset_arcsec, read_bpz, shift  # noqa: E402

PRIMARY = "macs0308-relics-glafic"
SECOND = "macs0308-relics-zitrin"
OUT = Path(__file__).resolve().parents[1] / "results" / "arcs" / "macs0308-relics"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--catalog", type=Path, required=True)
    ap.add_argument("--bpz", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--zout", type=Path, help="where to write the converted photo-z table")
    args = ap.parse_args(argv)
    zout = bpz_zout(read_bpz(args.bpz))
    off = offset_arcsec(zout, Table.read(args.catalog))
    off["shift_applied"] = shift(zout, off)
    zpath = args.zout or args.bpz.with_name("macs0308_bpz_zout.fits")
    zout.write(zpath, overwrite=True)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "photoz_offset.json").write_text(
        json.dumps({"bpz": args.bpz.name, "n_bpz_galaxies": len(zout), **off}, indent=1) + "\n"
    )
    print(json.dumps(off))
    decide(PRIMARY, SECOND, args.out, args.catalog, zpath, "results/arcs/macs0308-relics/PLAN.md")
    return 0


if __name__ == "__main__":
    np.seterr(invalid="ignore")
    raise SystemExit(main())
