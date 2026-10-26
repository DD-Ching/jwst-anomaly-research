"""Footprint check of the unused Abell 2744 NIRCam epoch pairs (metadata only).

Plan: results/transients/a2744_epoch_footprint/PLAN.md. Reads the pairs of
``results/transients/epoch_inventory/inventory.json``, makes one MAST query for the level-3
NIRCam footprints (``s_region``) of their epochs, and measures each pair's overlap and core
overlap per shared filter on a tangent-plane point grid. Every area is ``derived``; the grid step,
core radius and class thresholds are ASSUMPTIONs of the plan.

    python scripts/a2744_epoch_footprint.py --out results/transients/a2744_epoch_footprint
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from matplotlib.path import Path as MplPath

sys.path.insert(0, str(Path(__file__).resolve().parent))
from epoch_inventory import epoch_key, filter_name  # noqa: E402

CENTRE = (3.58750, -30.39750)  # plan rule 1 fallback (the inventory stored no centre)
GRID_STEP_ARCSEC = 1.0  # ASSUMPTION (rule 2)
HALF_WIDTH_ARCMIN = 5.0  # ASSUMPTION (rule 2)
CORE_RADIUS_ARCMIN = 1.0  # ASSUMPTION (rule 3)
CORE_FRAC = 0.5  # ASSUMPTION (rule 4): half the core disc
MIN_OVERLAP_ARCMIN2 = 1.0  # ASSUMPTION (rule 4)
FRAME_WORDS = {"ICRS", "FK5", "J2000", "GALACTIC", "ECLIPTIC", "UNKNOWNFRAME", "TOPOCENTER"}


def parse_s_region(s: str) -> list[np.ndarray]:
    """Polygons (N x 2 arrays of RA, Dec in deg) of an STC-S ``POLYGON`` string."""
    polys = []
    for part in re.split(r"POLYGON", s, flags=re.IGNORECASE)[1:]:
        nums = [float(t) for t in part.split() if t.upper() not in FRAME_WORDS]
        if len(nums) < 6 or len(nums) % 2:
            raise ValueError(f"bad polygon in s_region: {s[:80]}")
        polys.append(np.array(nums).reshape(-1, 2))
    if not polys:
        raise ValueError(f"no POLYGON in s_region: {s[:80]}")
    return polys


def gnomonic(ra: np.ndarray, dec: np.ndarray, ra0: float, dec0: float) -> np.ndarray:
    """Tangent-plane (xi, eta) in arcmin of (ra, dec) about (ra0, dec0), all in degrees."""
    ra, dec, ra0, dec0 = (np.radians(np.asarray(v, float)) for v in (ra, dec, ra0, dec0))
    cosc = np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0)
    xi = np.cos(dec) * np.sin(ra - ra0) / cosc
    eta = (np.cos(dec0) * np.sin(dec) - np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)) / cosc
    return np.degrees(np.column_stack([xi, eta])) * 60.0


def grid_points(half_width_arcmin: float, step_arcsec: float) -> tuple[np.ndarray, float]:
    """Cell centres of a square tangent-plane grid and the area of one cell (arcmin^2)."""
    step = step_arcsec / 60.0
    n = int(round(2 * half_width_arcmin / step))
    c = -half_width_arcmin + step * (np.arange(n) + 0.5)
    xx, yy = np.meshgrid(c, c)
    return np.column_stack([xx.ravel(), yy.ravel()]), step**2


def inside(points: np.ndarray, polys_xy: list[np.ndarray]) -> np.ndarray:
    """Mask of points inside any of the polygons (rule 2)."""
    m = np.zeros(len(points), bool)
    for p in polys_xy:
        m |= MplPath(p).contains_points(points)
    return m


def classify(per_filter: dict[str, dict]) -> str:
    """Rule 4 class of a pair from its per-filter {'a_ov', 'a_core'} areas (arcmin^2)."""
    core_min = CORE_FRAC * np.pi * CORE_RADIUS_ARCMIN**2
    if any(v["a_core"] >= core_min for v in per_filter.values()):
        return "core"
    if any(v["a_ov"] >= MIN_OVERLAP_ARCMIN2 for v in per_filter.values()):
        return "overlap"
    return "none"


def outcome(pairs: list[dict]) -> tuple[str, dict | None]:
    """Plan outcome G / O / N and, for G, the pair named as the next grid input."""
    core_min = CORE_FRAC * np.pi * CORE_RADIUS_ARCMIN**2
    g = [p for p in pairs if sum(v["a_core"] >= core_min for v in p["per_filter"].values()) >= 2]
    if g:
        best = max(
            g, key=lambda p: (sum(v["a_core"] for v in p["per_filter"].values()), p["dt_days"])
        )
        return "G", best
    if any(p["class"] in {"core", "overlap"} for p in pairs):
        return "O", None
    return "N", None


def query_footprints(epochs: set[str]) -> dict[tuple[str, str], list[str]]:
    """s_region strings per (epoch, filter) from one MAST query (level-3 NIRCam images)."""
    from astroquery.mast import Observations

    proposals = sorted({e.split("-")[0].lstrip("0") for e in epochs})
    tab = Observations.query_criteria(
        objectname="Abell 2744",
        radius="3 arcmin",
        obs_collection="JWST",
        instrument_name="NIRCAM/IMAGE",
        calib_level=3,
        dataproduct_type="image",
        proposal_id=proposals,
    )
    out: dict[tuple[str, str], list[str]] = defaultdict(list)
    for r in tab:
        key = epoch_key(str(r["obs_id"]))
        f = filter_name(str(r["filters"]))
        if key is None or not f or "-".join(key) not in epochs:
            continue
        out[("-".join(key), f)].append(str(r["s_region"]))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument(
        "--inventory",
        type=Path,
        default=Path("results/transients/epoch_inventory/inventory.json"),
    )
    args = ap.parse_args()

    pairs_in = json.loads(args.inventory.read_text())["fields"]["abell2744"]["new_pairs"]
    epochs = {p["early"] for p in pairs_in} | {p["late"] for p in pairs_in}
    regions = query_footprints(epochs)

    pts, cell = grid_points(HALF_WIDTH_ARCMIN, GRID_STEP_ARCSEC)
    core = np.hypot(pts[:, 0], pts[:, 1]) <= CORE_RADIUS_ARCMIN
    masks: dict[tuple[str, str], np.ndarray] = {}
    for key, srs in regions.items():
        polys = [gnomonic(p[:, 0], p[:, 1], *CENTRE) for s in srs for p in parse_s_region(s)]
        masks[key] = inside(pts, polys)

    pairs = []
    for p in pairs_in:
        per, missing = {}, []
        for f in p["shared_filters"]:
            a, b = masks.get((p["early"], f)), masks.get((p["late"], f))
            if a is None or b is None:
                missing.append(f)  # rule 5; recorded so a "none" cannot hide a missing footprint
                continue
            ov = a & b
            per[f] = {
                "a_ov": round(float(ov.sum() * cell), 3),
                "a_core": round(float((ov & core).sum() * cell), 3),
            }
        pairs.append({**p, "per_filter": per, "missing_filters": missing, "class": classify(per)})

    code, best = outcome(pairs)
    epoch_area = {
        f"{e}:{f}": {
            "area": round(float(m.sum() * cell), 3),
            "core": round(float((m & core).sum() * cell), 3),
        }
        for (e, f), m in sorted(masks.items())
    }
    counts = {c: sum(p["class"] == c for p in pairs) for c in ("core", "overlap", "none")}
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "footprints.json").write_text(
        json.dumps(
            {
                "centre": CENTRE,
                "core_disc_arcmin2": round(np.pi * CORE_RADIUS_ARCMIN**2, 3),
                "n_pairs": len(pairs),
                "class_counts": counts,
                "n_pairs_with_missing_filters": sum(bool(p["missing_filters"]) for p in pairs),
                "outcome": code,
                "next_grid_pair": None
                if best is None
                else {k: best[k] for k in ("early", "late", "dt_days", "per_filter")},
                "epoch_filter_area": epoch_area,
                "pairs": pairs,
            },
            indent=1,
        )
        + "\n"
    )
    print(f"outcome {code}; classes {counts}")
    if best:
        print(f"next grid pair {best['early']} -> {best['late']} ({best['dt_days']} d)")


if __name__ == "__main__":
    main()
