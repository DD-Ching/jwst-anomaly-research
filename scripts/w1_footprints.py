"""Areas, pairwise overlaps and centroids of a program's NIRCam level-3 footprints (W1 blank-field
pointing choice; results/w1_shear_ceers/PLAN.md).

    python scripts/w1_footprints.py --program 1345 --filter F200W \
        --out results/w1_shear_ceers/footprints.json

One MAST ``Mast.Caom.Filtered`` query; each ``s_region`` (one or more POLYGONs) is rasterised on a
flat tangent-plane grid of ``PIX_ARCMIN``. Footprints ``observed``; areas and overlaps ``derived``
(the grid step and the flat-sky approximation are ASSUMPTIONs; overlaps below 0.05 arcmin² are not
listed). ``dropped``: the pointings ``drop_overlapping`` leaves out.
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
TILE_HALF = 140.0  # arcsec; half-diagonal 198" stays inside the screen's 200" radius (ASSUMPTION)
TILE_MIN_FRAC = 0.3  # a tile is kept when its pointing covers >= this share (ASSUMPTION)
MAX_OVERLAP = 0.6  # arcmin²; pairwise rule of results/w1_shear_goodsn/PLAN.md (ASSUMPTION)


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
        # frame words (ICRS, J2000, ...) are skipped: only plain numbers are vertices
        v = [float(t) for t in chunk.split() if re.fullmatch(r"[-+]?[0-9]*\.?[0-9]+", t)]
        out.append(np.asarray(v).reshape(-1, 2))
    return out


def first_ra(s_region: str) -> float:
    """RA (deg) of the first vertex of an STC-S region."""
    return float(polygons(s_region)[0][0, 0])


def drop_overlapping(overlap: dict[str, float], max_overlap: float | None = None) -> list[str]:
    """Pointings to leave out so that no kept pair overlaps by more than ``max_overlap``: drop the
    pointing with the largest summed excess overlap, repeat (ties: name order)."""
    max_overlap = MAX_OVERLAP if max_overlap is None else max_overlap
    pairs = {tuple(k.split("|")): v for k, v in overlap.items() if v > max_overlap}
    dropped = []
    while pairs:
        load: dict[str, float] = {}
        for (a, b), v in pairs.items():
            load[a] = load.get(a, 0.0) + v
            load[b] = load.get(b, 0.0) + v
        worst = max(sorted(load), key=load.get)
        dropped.append(worst)
        pairs = {k: v for k, v in pairs.items() if worst not in k}
    return sorted(dropped)


def _masks(rows: list[dict]):
    """Rasterised footprints: (grid offsets in arcmin, {name: bool mask}, ra0, de0, cos dec0)."""
    polys = {r["obs_id"].split("_nircam")[0]: polygons(r["s_region"]) for r in rows}
    allv = np.vstack([p for ps in polys.values() for p in ps])
    ra0 = allv[0, 0]
    ra0 += np.mean((allv[:, 0] - ra0 + 180.0) % 360.0 - 180.0)  # RA mean safe across 0/360
    de0 = allv[:, 1].mean()
    cosd = np.cos(np.radians(de0))
    n = int(2 * HALF_ARCMIN / PIX_ARCMIN)
    ax = -HALF_ARCMIN + PIX_ARCMIN * np.arange(n)
    grid = np.stack(np.meshgrid(ax, ax, indexing="ij"), -1).reshape(-1, 2)
    masks = {}
    for name, ps in polys.items():
        m = np.zeros(len(grid), bool)
        for p in ps:
            dra = (p[:, 0] - ra0 + 180.0) % 360.0 - 180.0
            xy = np.c_[dra * cosd * 60, (p[:, 1] - de0) * 60]
            if np.abs(xy).max() > HALF_ARCMIN - PIX_ARCMIN:
                raise ValueError(f"{name}: footprint beyond the +-{HALF_ARCMIN} arcmin grid")
            box = np.all((grid > xy.min(0) - PIX_ARCMIN) & (grid < xy.max(0) + PIX_ARCMIN), 1)
            m[box] |= MplPath(xy).contains_points(grid[box])
        masks[name] = m
    return grid, masks, ra0, de0, cosd


def footprints(rows: list[dict]) -> dict:
    """Area per pointing, pairwise overlaps and footprint centroids (degrees)."""
    grid, masks, ra0, de0, cosd = _masks(rows)
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
        centroid[name] = [round((ra0 + x / 60 / cosd) % 360.0, 5), round(de0 + y / 60, 5)]
    return {
        "area_arcmin2": {k: round(float(masks[k].sum() * px2), 2) for k in names},
        "overlap_arcmin2": overlap,
        "centroid": centroid,
        "exptime_s": {r["obs_id"].split("_nircam")[0]: r["t_exptime"] for r in rows},
        "dropped": drop_overlapping(overlap),
        "assumptions": {
            "pix_arcmin": PIX_ARCMIN,
            "min_overlap_arcmin2": MIN_OVERLAP,
            "max_overlap_arcmin2": MAX_OVERLAP,
        },
    }


def tile_grid(rows: list[dict], keep: list[str], half_arcsec: float, min_frac: float) -> list[dict]:
    """Square tiles of half-width ``half_arcsec`` on one RA/Dec-aligned grid over the union of the
    ``keep`` footprints (results/w1_shear_primer/PLAN.md). Each tile goes to the kept pointing that
    covers most of it (ties: name order) and is listed when that share is >= ``min_frac``. Tiles do
    not overlap, so no sky is screened twice; sky of a tile that only another pointing covers is
    lost (conservative)."""
    grid, masks, ra0, de0, cosd = _masks(rows)
    masks = {k: masks[k] for k in sorted(keep)}
    union = np.logical_or.reduce(list(masks.values()))
    side = 2 * half_arcsec / 60.0  # arcmin
    pts = grid[union]
    x0, y0 = pts.min(0)
    nx, ny = (np.ceil((pts.max(0) - pts.min(0)) / side).astype(int)) + 1
    tiles = []
    for ix in range(nx):
        for iy in range(ny):
            cx, cy = x0 + (ix + 0.5) * side, y0 + (iy + 0.5) * side
            dx, dy = grid[:, 0] - cx, grid[:, 1] - cy
            h = side / 2
            inside = (dx >= -h) & (dx < h) & (dy >= -h) & (dy < h)  # half-open, as in_box
            n = inside.sum()
            share = {k: float((m & inside).sum() / n) for k, m in masks.items()}
            best = max(sorted(share), key=share.get)  # first maximum: name order
            if share[best] >= min_frac:
                tiles.append(
                    {
                        "tile": f"x{ix}y{iy}",
                        "pointing": best,
                        "centre": [
                            round((ra0 + cx / 60 / cosd) % 360.0, 5),
                            round(de0 + cy / 60, 5),
                        ],
                        "share": round(share[best], 3),
                    }
                )
    return tiles


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--program", required=True)
    ap.add_argument("--filter", default="F200W")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument(
        "--ra-range",
        nargs=2,
        type=float,
        metavar=("MIN", "MAX"),
        help="keep pointings whose first vertex RA is in [MIN, MAX) deg (multi-field programs)",
    )
    ap.add_argument("--tile-keep", nargs="+", help="pointings to tile (large mosaics)")
    ap.add_argument("--tile-half", type=float, default=TILE_HALF, help="tile half-width, arcsec")
    ap.add_argument("--tile-min-frac", type=float, default=TILE_MIN_FRAC)
    args = ap.parse_args(argv)
    rows = query(args.program, args.filter)
    if args.ra_range:
        rows = [r for r in rows if args.ra_range[0] <= first_ra(r["s_region"]) < args.ra_range[1]]
    res = footprints(rows)
    if args.tile_keep:
        res["tiles"] = tile_grid(rows, args.tile_keep, args.tile_half, args.tile_min_frac)
        res["assumptions"] |= {
            "tile_half_arcsec": args.tile_half,
            "tile_min_frac": args.tile_min_frac,
        }
    args.out.write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res["overlap_arcmin2"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
