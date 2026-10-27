"""Post hoc checks of the MACS J0138 two-epoch grid (results/transients/macs0138_grid; not planned).

The same three checks as ``a1689_grid_posthoc`` (hit-group coverage and bright-star labels, the G1
split of #578, recovery by expected calibrated significance), on the F150W / F444W pair
2345 o005 → 6549 o003. Thresholds are the post hoc ASSUMPTIONs of #576 / #577 / #579. All numbers
are ``derived``.

    python scripts/m0138_grid_posthoc.py --data data/transients/macs0138 \
        --run results/transients/macs0138_grid
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from astropy.io import fits
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import a1689_grid_posthoc as ap  # noqa: E402
from m2129_grid_posthoc import GHOST_PA, GHOST_PA_REF, classify, first_rows  # noqa: E402

E1 = "jw02345-o005_t003_nircam_clear-"
E2 = "jw06549-o003_t001_nircam_clear-"
BANDS = ("F150W", "F444W")
STAR_BAND = "F444W"


def configure() -> None:
    """Point the reused Abell 1689 functions at this pair (they read module-level names)."""
    ap.E1, ap.E2, ap.BANDS, ap.STAR_BAND = E1, E2, BANDS, STAR_BAND


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--data", type=Path, required=True)
    p.add_argument("--run", type=Path, required=True)
    a = p.parse_args(argv)
    configure()
    hits, grid = Table.read(a.run / "hits.ecsv"), Table.read(a.run / "grid.ecsv")
    rows = first_rows(hits)
    cov, star = ap.coverage(a.data, rows), ap.bright_stars(grid, rows)
    pa = {
        tag: round(float(fits.getheader(a.data / f"{s}f444w_i2d.fits", "SCI")["PA_APER"]), 2)
        for tag, s in (("e1", E1), ("e2", E2))
    }
    res = {
        "provenance": "derived",
        "post_hoc": True,
        "coverage": cov,
        "bright_star": star,
        "pa_aper": pa,
        "ghost_pa_rotated_deg": round((GHOST_PA + pa["e1"] - GHOST_PA_REF) % 360.0, 2),
        "classification": classify(cov, star, pa["e1"]),
        "injection": ap.injection_sigma(a.data, a.run),
        "mosaic_bytes": {q.name: q.stat().st_size for q in sorted(a.data.glob("*_i2d.fits"))},
    }
    (a.run / "posthoc.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: res[k] for k in ("classification", "injection")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
