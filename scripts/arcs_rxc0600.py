"""`arcs` on RXC J0600-2007 (MACS J0600.1-2008): an independent test of the S1063 anti excess.

Plan (frozen before the run): ``results/arcs/rxc0600-relics/PLAN.md``. Runs
``lens_consistency.py validate`` and ``arcs`` (D-024 defaults) with the RELICS GLAFIC v3 (primary)
and Lenstool v1 (second) maps on the program 6796 o010 F200W catalogue and applies the frozen
outcome rule (``arcs_macs0647.outcome``) through ``arcs_rxj2129.decide``. Shapes are ``observed``,
predicted axes ``model_prediction``, classes, counts and P values ``derived``; thresholds are
ASSUMPTIONs.

    python scripts/arcs_rxc0600.py --catalog <jw06796-o010_t001 F200W _cat.ecsv> \\
        --photoz <rxcj0600-grizli-v7.2-fix.eazypy.zout.fits>
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from arcs_rxj2129 import decide  # noqa: E402

PRIMARY = "rxc0600-relics-glafic"
SECOND = "rxc0600-relics-lenstool"
OUT = Path(__file__).resolve().parents[1] / "results" / "arcs" / "rxc0600-relics"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--catalog", type=Path, required=True)
    ap.add_argument("--photoz", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    decide(
        PRIMARY, SECOND, args.out, args.catalog, args.photoz, "results/arcs/rxc0600-relics/PLAN.md"
    )
    return 0


if __name__ == "__main__":
    np.seterr(invalid="ignore")
    raise SystemExit(main())
