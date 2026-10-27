"""Two-epoch NIRCam inventory with a core-footprint test on 15 more lensing clusters (metadata).

Plan: results/transients/epoch_inventory_ext/PLAN.md. Same query, epoch definition and new-pair
rule (a)-(d) as ``epoch_inventory.py``; each new pair is then classed ``core`` / ``overlap`` /
``none`` from the overlap of both epochs' ``s_region`` footprints with a 1' disc about the
MAST-resolved field centre (``a2744_epoch_footprint.py`` geometry). Areas are ``derived``.

    python scripts/epoch_inventory_ext.py --out results/transients/epoch_inventory_ext
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from a2744_epoch_footprint import gnomonic, grid_points, inside, parse_s_region  # noqa: E402
from epoch_inventory import (  # noqa: E402
    EPOCH_RE,
    build_epochs,
    epoch_key,
    filter_name,
    new_pairs,
    used_epochs,
)

FIELDS_EXT = {
    "macs0138": "MACS J0138.0-2155",
    "rxj0152": "RX J0152.7-1357",
    "macs0358": "MACS J0358.8-2955",
    "rxcj2211": "RXC J2211.7-0350",
    "spt0546": "SPT-CL J0546-5345",
    "macs2129": "MACS J2129.4-0741",
    "macs1206": "MACS J1206.2-0847",
    "abell1689": "Abell 1689",
    "macs0329": "MACS J0329.6-0211",
    "abell383": "Abell 383",
    "abell209": "Abell 209",
    "macs0308": "MACS J0308.9+2645",
    "sdss1004": "SDSS J1004+4112",
    "sdss1029": "SDSS J1029+2623",
    "sdss2222": "SDSS J2222+2745",
}
RADIUS = "3 arcmin"
CORE_RADIUS_ARCMIN = 1.0  # ASSUMPTION (plan)
CORE_MIN_ARCMIN2 = (
    1.0  # ASSUMPTION (plan): "core" when the pair overlap inside the disc > 1 arcmin^2
)
GRID_STEP_ARCSEC = 2.0  # ASSUMPTION: cell 1/900 arcmin^2, far below the 1 arcmin^2 threshold
HALF_WIDTH_ARCMIN = 6.0  # covers a 3' cone plus a full NIRCam mosaic edge
INVENTORY_DIRS = ("epoch_inventory", "epoch_inventory_ext")  # list every epoch, not used ones


def classify(per_filter: dict[str, dict]) -> str:
    """Plan class of a pair from its per-filter {'a_ov', 'a_core'} areas (arcmin^2).

    ``core``: a_core > CORE_MIN_ARCMIN2 in at least one shared filter; else ``overlap``: a_ov > 0
    in at least one shared filter; else ``none``.
    """
    if any(v["a_core"] > CORE_MIN_ARCMIN2 for v in per_filter.values()):
        return "core"
    if any(v["a_ov"] > 0 for v in per_filter.values()):
        return "overlap"
    return "none"


def used_excluding_inventories(root: Path) -> set[tuple[str, str]]:
    """Used epochs under ``root`` minus everything named only in the inventory directories."""
    used = set()
    for sub in root.iterdir():
        if sub.name in INVENTORY_DIRS:
            continue
        if sub.is_dir():
            used |= used_epochs(sub)
        elif sub.suffix in {".md", ".json", ".ecsv", ".csv", ".txt"}:
            used |= set(EPOCH_RE.findall(sub.read_text(errors="ignore")))
    return used


def query_field(name: str) -> dict:
    """Resolved centre and level-3 NIRCam image rows (with s_region) of one field."""
    from astroquery.mast import Mast, Observations

    try:
        c = Mast.resolve_object(name)
    except Exception as exc:  # recorded as unresolved, per the plan
        return {"unresolved": f"{type(exc).__name__}: {exc}"[:200]}
    tab = Observations.query_criteria(
        coordinates=c,
        radius=RADIUS,
        obs_collection="JWST",
        instrument_name="NIRCAM/IMAGE",
        calib_level=3,
        dataproduct_type="image",
    )
    rows = []
    for r in tab:
        if np.ma.is_masked(r["t_min"]) or np.ma.is_masked(r["t_obs_release"]):
            continue
        rows.append(
            {
                "obs_id": str(r["obs_id"]),
                "filters": str(r["filters"]),
                "t_min": float(r["t_min"]),
                "t_obs_release": float(r["t_obs_release"]),
                "s_region": str(r["s_region"]),
            }
        )
    return {"centre": (float(c.ra.deg), float(c.dec.deg)), "rows": rows}


def pair_areas(
    pairs: list[dict], rows: list[dict], centre: tuple[float, float]
) -> tuple[list[dict], dict]:
    """Per-filter pair overlap and core-overlap areas, and each pair's class."""
    regions: dict[tuple[str, str], list[str]] = defaultdict(list)
    for r in rows:
        k, f = epoch_key(r["obs_id"]), filter_name(r["filters"])
        if k and f:
            regions[("-".join(k), f)].append(r["s_region"])
    pts, cell = grid_points(HALF_WIDTH_ARCMIN, GRID_STEP_ARCSEC)
    core = np.hypot(pts[:, 0], pts[:, 1]) <= CORE_RADIUS_ARCMIN
    masks = {}
    for key, srs in regions.items():
        polys = [gnomonic(p[:, 0], p[:, 1], *centre) for s in srs for p in parse_s_region(s)]
        masks[key] = inside(pts, polys)
    out = []
    for p in pairs:
        per, missing = {}, []
        for f in p["shared_filters"]:
            a, b = masks.get((p["early"], f)), masks.get((p["late"], f))
            if a is None or b is None:
                missing.append(f)
                continue
            ov = a & b
            per[f] = {
                "a_ov": round(float(ov.sum() * cell), 3),
                "a_core": round(float((ov & core).sum() * cell), 3),
            }
        out.append({**p, "per_filter": per, "missing_filters": missing, "class": classify(per)})
    epoch_core = {
        f"{e}:{f}": round(float((m & core).sum() * cell), 3) for (e, f), m in sorted(masks.items())
    }
    return out, epoch_core


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--used-root", type=Path, default=Path("results/transients"))
    args = ap.parse_args()

    from astropy.time import Time

    now = float(Time.now().mjd)
    used = used_excluding_inventories(args.used_root)
    with ThreadPoolExecutor(max_workers=4) as pool:
        res = dict(zip(FIELDS_EXT, pool.map(query_field, FIELDS_EXT.values()), strict=True))

    summary, upcoming = {}, {}
    for field, r in res.items():
        if "unresolved" in r:
            summary[field] = {"name": FIELDS_EXT[field], "unresolved": r["unresolved"]}
            continue
        epochs = build_epochs(r["rows"])
        pairs, epoch_core = pair_areas(new_pairs(epochs, used, now), r["rows"], r["centre"])
        up = sorted(
            ("-".join(k), round(e["release"], 1)) for k, e in epochs.items() if e["release"] > now
        )
        if up:
            upcoming[field] = [{"epoch": k, "release_mjd": m} for k, m in up]
        summary[field] = {
            "name": FIELDS_EXT[field],
            "centre_deg": [round(v, 5) for v in r["centre"]],
            "n_rows": len(r["rows"]),
            "epochs": {
                "-".join(k): {
                    "t_min_mjd": round(e["t_min"], 2),
                    "public": e["release"] <= now,
                    "used": k in used,
                    "filters": sorted(e["filters"]),
                }
                for k, e in sorted(epochs.items(), key=lambda kv: kv[1]["t_min"])
            },
            "epoch_filter_core_arcmin2": epoch_core,
            "class_counts": {
                c: sum(p["class"] == c for p in pairs) for c in ("core", "overlap", "none")
            },
            "new_pairs": pairs,
        }
    args.out.mkdir(parents=True, exist_ok=True)
    meta = {
        "queried_mjd": round(now, 4),
        "radius": RADIUS,
        "core_radius_arcmin": CORE_RADIUS_ARCMIN,
        "core_min_arcmin2": CORE_MIN_ARCMIN2,
        "n_used_epochs": len(used),
    }
    (args.out / "inventory.json").write_text(
        json.dumps({"meta": meta, "fields": summary}, indent=1) + "\n"
    )
    (args.out / "upcoming.json").write_text(json.dumps(upcoming, indent=1) + "\n")
    for f, s in summary.items():
        if "unresolved" in s:
            print(f"{f}: unresolved")
        else:
            print(f"{f}: {len(s['epochs'])} epochs, classes {s['class_counts']}")


if __name__ == "__main__":
    main()
