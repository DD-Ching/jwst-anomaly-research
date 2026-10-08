"""Bright Gaia DR3 stars around a field, for the radial screen's spike veto (D-043).

    python scripts/gaia_stars.py RA DEC --radius-arcmin 4 --gmax 17 --out stars.ecsv

One VizieR cone query (I/355/gaiadr3). Writes ``ra, dec, mag`` (G) as ECSV, the format
``exotic_screens.py radial --spike-stars`` reads. Provenance ``observed``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("ra", type=float)
    ap.add_argument("dec", type=float)
    ap.add_argument("--radius-arcmin", type=float, default=4.0)
    ap.add_argument("--gmax", type=float, default=17.0, help="faintest G kept (ASSUMPTION)")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    from astroquery.vizier import Vizier

    v = Vizier(
        columns=["RA_ICRS", "DE_ICRS", "Gmag"],
        column_filters={"Gmag": f"<{args.gmax}"},
        row_limit=-1,
    )
    res = v.query_region(
        SkyCoord(args.ra * u.deg, args.dec * u.deg),
        radius=args.radius_arcmin * u.arcmin,
        catalog="I/355/gaiadr3",
    )
    g = res[0] if len(res) else Table({"RA_ICRS": [], "DE_ICRS": [], "Gmag": []})
    out = Table(
        {
            "ra": np.asarray(g["RA_ICRS"], float),
            "dec": np.asarray(g["DE_ICRS"], float),
            "mag": np.asarray(g["Gmag"], float),
        }
    )
    out.meta.update(
        provenance="observed",
        source=(
            f"VizieR I/355/gaiadr3 cone {args.ra} {args.dec} r={args.radius_arcmin}' G<{args.gmax}"
        ),
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.write(args.out, format="ascii.ecsv", overwrite=True)
    print(f"{len(out)} stars -> {args.out}")


if __name__ == "__main__":
    main()
