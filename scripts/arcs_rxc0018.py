"""`arcs` on RXC J0018.5+1626 (MACS J0018.5+1626): an independent test of the S1063 anti excess.

Plan (frozen before the run): ``results/arcs/rxc0018-relics/PLAN.md``. Runs
``lens_consistency.py validate`` and ``arcs`` (D-024 defaults) with the RELICS Lenstool v1 maps
(the only public model) on the program 6882 o001 F200W catalogue and applies the frozen outcome
rule (``arcs_macs0647.outcome``) through ``arcs_rxj2129.decide``, with the one model passed as
primary and second. Shapes are ``observed``, predicted axes ``model_prediction``, classes, counts
and P values ``derived``; thresholds are ASSUMPTIONs.

    python scripts/arcs_rxc0018.py --catalog <jw06882-o001_t001 F200W _cat.ecsv> \\
        --photoz <rj0018p1626-grizli-v7.5-fix.eazypy.zout.fits>
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from arcs_rxj2129 import decide  # noqa: E402

PRIMARY = "rxc0018-relics-lenstool"
OUT = Path(__file__).resolve().parents[1] / "results" / "arcs" / "rxc0018-relics"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--catalog", type=Path, required=True)
    ap.add_argument("--photoz", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    decide(
        PRIMARY, PRIMARY, args.out, args.catalog, args.photoz, "results/arcs/rxc0018-relics/PLAN.md"
    )
    return 0


if __name__ == "__main__":
    np.seterr(invalid="ignore")
    raise SystemExit(main())
