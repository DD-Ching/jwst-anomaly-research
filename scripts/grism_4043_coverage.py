"""Which SMACS 0723 sources fall inside the program 4043 NIRCam F444W grism footprints.

Inputs are recorded in ``results/grism_4043/RESULT.md``:
- an ObsCore VOTable from the MAST CAOM TAP service (``obs_id LIKE 'jw04043%'``), whose
  level-3 grism rows carry ``s_region`` polygons (one row per footprint piece);
- the 2736 NIRCam F200W level-3 catalog (jwst 3.0.0) for the positions of the named sources;
- ``results/smacs_jwst3/candidates.csv`` for the candidate records.

Being inside the footprint is necessary, not sufficient, for a spectrum (ASSUMPTION; the trace
can leave the detector and the 4043 extraction uses its own direct-image catalog).
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import numpy as np

GRISM_OBS = ("jw04043-o001_t001_nircam_f444w-grismr", "jw04043-o001_t001_nircam_f444w-grismc")


def parse_polygons(s_region: str) -> list[np.ndarray]:
    """Split an STC-S ``s_region`` string into (N, 2) arrays of (ra, dec) vertices."""
    polys = []
    for chunk in re.split(r"POLYGON", s_region):
        vals = [float(v) for v in chunk.split() if re.fullmatch(r"[-+0-9.eE]+", v)]
        if len(vals) >= 6:
            polys.append(np.asarray(vals).reshape(-1, 2))
    return polys


def _tangent(ra, dec, ra0, dec0):
    """Gnomonic projection in degrees; adequate for arcminute fields."""
    ra, dec, ra0, dec0 = map(np.radians, (ra, dec, ra0, dec0))
    cosc = np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0)
    x = np.cos(dec) * np.sin(ra - ra0) / cosc
    y = (np.cos(dec0) * np.sin(dec) - np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)) / cosc
    return np.degrees(x), np.degrees(y)


def inside(ra: float, dec: float, poly: np.ndarray) -> bool:
    """Even-odd point-in-polygon test in the tangent plane at the point."""
    x, y = _tangent(poly[:, 0], poly[:, 1], ra, dec)
    n = len(x)
    hit = False
    for i in range(n):
        j = (i - 1) % n
        if (y[i] > 0) != (y[j] > 0):
            xc = x[i] + (0 - y[i]) * (x[j] - x[i]) / (y[j] - y[i])
            if xc > 0:
                hit = not hit
    return hit


def edge_distance_arcsec(ra: float, dec: float, poly: np.ndarray) -> float:
    """Distance from the point to the nearest polygon edge, in arcsec."""
    x, y = _tangent(poly[:, 0], poly[:, 1], ra, dec)
    best = np.inf
    for i in range(len(x)):
        j = (i + 1) % len(x)
        dx, dy = x[j] - x[i], y[j] - y[i]
        seg = dx * dx + dy * dy
        t = 0.0 if seg == 0 else np.clip(-(x[i] * dx + y[i] * dy) / seg, 0, 1)
        best = min(best, np.hypot(x[i] + t * dx, y[i] + t * dy))
    return float(best * 3600)


def coverage(ra: float, dec: float, footprints: dict[str, list[np.ndarray]]) -> dict:
    out = {}
    for obs, polys in footprints.items():
        hits = [p for p in polys if inside(ra, dec, p)]
        out[obs] = {
            "inside": bool(hits),
            "edge_arcsec": round(min(edge_distance_arcsec(ra, dec, p) for p in hits), 1)
            if hits
            else None,
        }
    return out


def load_footprints(votable: Path) -> dict[str, list[np.ndarray]]:
    from astropy.io.votable import parse_single_table

    t = parse_single_table(str(votable)).to_table()
    fp: dict[str, list[np.ndarray]] = {o: [] for o in GRISM_OBS}
    for row in t:
        if row["obs_id"] in fp and int(row["calib_level"]) == 3:
            fp[row["obs_id"]].extend(parse_polygons(str(row["s_region"])))
    return fp


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--obscore", type=Path, required=True)
    ap.add_argument("--f200w-cat", type=Path, required=True)
    ap.add_argument("--candidates", type=Path, default=Path("results/smacs_jwst3/candidates.csv"))
    ap.add_argument("--labels", default="1816,384", help="extra F200W catalog labels to place")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    from astropy.table import Table

    fp = load_footprints(args.obscore)
    cat = Table.read(args.f200w_cat)
    rows = []
    with args.candidates.open() as fh:
        for r in csv.DictReader(fh):
            rows.append(
                {
                    "name": r["record"],
                    "ra": float(r["ra"]),
                    "dec": float(r["dec"]),
                    "source": "candidates.csv",
                }
            )
    for lab in [int(v) for v in args.labels.split(",") if v]:
        m = cat[cat["label"] == lab]
        c = m["sky_centroid"][0]
        rows.append(
            {
                "name": f"jw02736-o001_t001_nircam_f200w_{lab}",
                "ra": float(c.ra.deg),
                "dec": float(c.dec.deg),
                "source": "2736 F200W cat (jwst 3.0.0)",
            }
        )
    # Field-wide fraction of catalog sources inside each footprint (derived).
    sky = cat["sky_centroid"]
    frac = {}
    for obs, polys in fp.items():
        n_in = sum(any(inside(float(s.ra.deg), float(s.dec.deg), p) for p in polys) for s in sky)
        frac[obs] = round(n_in / len(sky), 3)
    for r in rows:
        r["coverage"] = coverage(r["ra"], r["dec"], fp)
    summary = {
        "n_polygons": {k: len(v) for k, v in fp.items()},
        "f200w_catalog_fraction_inside": frac,
        "n_catalog": len(cat),
        "records": rows,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "records"}))
    for r in rows:
        print(r["name"].split("_", 3)[-1], {k.split("_")[-1]: v for k, v in r["coverage"].items()})


if __name__ == "__main__":
    main()
