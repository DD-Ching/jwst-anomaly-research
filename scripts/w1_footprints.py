"""Areas, pairwise overlaps and centroids of a program's NIRCam level-3 footprints (W1 blank-field
pointing choice; results/w1_shear_ceers/PLAN.md).

    python scripts/w1_footprints.py --program 1345 --filter F200W \
        --out results/w1_shear_ceers/footprints.json

One MAST ``Mast.Caom.Filtered`` query; each ``s_region`` (one or more POLYGONs) is rasterised on a
flat tangent-plane grid of ``PIX_ARCMIN``. Footprints ``observed``; areas and overlaps ``derived``
(the grid step and the flat-sky approximation are ASSUMPTIONs; overlaps below 0.05 arcmin² are not
listed).
"""

from __future__ import annotations

import argparse
import json
import re
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
from matplotlib.path import Path as MplPath

MAST = "https://mast.stsci.edu/api/v0/invoke"
PIX_ARCMIN = 0.02  # raster step (1.2"; ASSUMPTION)
HALF_ARCMIN = 25.0  # grid half-width around the footprints' mean centre
MIN_OVERLAP = 0.05  # arcmin²; smaller overlaps are not listed


def query(program: str, band: str) -> list[dict]:
    req = {
        "service": "Mast.Caom.Filtered",
        "format": "json",
        "params": {
            "columns": "obs_id,filters,t_exptime,s_region",
            "filters": [
                {"paramName": "proposal_id", "values": [program]},
                {"paramName": "obs_collection", "values": ["JWST"]},
                {"paramName": "instrument_name", "values": ["NIRCAM/IMAGE"]},
                {"paramName": "calib_level", "values": [3]},
                {"paramName": "filters", "values": [band]},
            ],
        },
    }
    data = urllib.parse.urlencode({"request": json.dumps(req)}).encode()
    with urllib.request.urlopen(MAST, data=data, timeout=120) as resp:  # noqa: S310
        return [r for r in json.load(resp)["data"] if r["filters"] == band]


def polygons(s_region: str) -> list[np.ndarray]:
    """(ra, dec) vertex arrays of every POLYGON in an STC-S string."""
    out = []
    for chunk in re.split(r"POLYGON", s_region)[1:]:
        v = [float(x) for x in re.findall(r"[-+]?[0-9]*\.?[0-9]+", chunk)]
        out.append(np.asarray(v).reshape(-1, 2))
    return out


def footprints(rows: list[dict]) -> dict:
    """Area per pointing, pairwise overlaps and footprint centroids (degrees)."""
    polys = {r["obs_id"].split("_nircam")[0]: polygons(r["s_region"]) for r in rows}
    allv = np.vstack([p for ps in polys.values() for p in ps])
    ra0, de0 = allv.mean(0)
    cosd = np.cos(np.radians(de0))
    n = int(2 * HALF_ARCMIN / PIX_ARCMIN)
    ax = -HALF_ARCMIN + PIX_ARCMIN * np.arange(n)
    grid = np.stack(np.meshgrid(ax, ax, indexing="ij"), -1).reshape(-1, 2)
    masks = {}
    for name, ps in polys.items():
        m = np.zeros(len(grid), bool)
        for p in ps:
            xy = np.c_[(p[:, 0] - ra0) * cosd * 60, (p[:, 1] - de0) * 60]
            box = np.all((grid > xy.min(0) - PIX_ARCMIN) & (grid < xy.max(0) + PIX_ARCMIN), 1)
            m[box] |= MplPath(xy).contains_points(grid[box])
        masks[name] = m
    names = sorted(masks)
    px2 = PIX_ARCMIN**2
    overlap = {}
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            o = float((masks[a] & masks[b]).sum() * px2)
            if o > MIN_OVERLAP:
                overlap[f"{a}|{b}"] = round(o, 2)
    centroid = {}
    for name in names:
        x, y = grid[masks[name]].mean(0)
        centroid[name] = [round(ra0 + x / 60 / cosd, 5), round(de0 + y / 60, 5)]
    return {
        "area_arcmin2": {k: round(float(masks[k].sum() * px2), 2) for k in names},
        "overlap_arcmin2": overlap,
        "centroid": centroid,
        "exptime_s": {r["obs_id"].split("_nircam")[0]: r["t_exptime"] for r in rows},
        "assumptions": {"pix_arcmin": PIX_ARCMIN, "min_overlap_arcmin2": MIN_OVERLAP},
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--program", required=True)
    ap.add_argument("--filter", default="F200W")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    res = footprints(query(args.program, args.filter))
    args.out.write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res["overlap_arcmin2"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
