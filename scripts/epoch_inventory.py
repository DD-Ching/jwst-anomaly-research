"""Inventory public JWST NIRCam imaging epochs on lensing-cluster fields (metadata only).

Plan: results/transients/epoch_inventory/PLAN.md. One batched MAST query per field; an epoch is
one (proposal, observation). A pair is "new" when both epochs are public, share a filter, are
>= 30 d apart and at least one epoch is absent from the epochs used by the D-027 transient runs.

    python scripts/epoch_inventory.py --out results/transients/epoch_inventory
"""

from __future__ import annotations

import argparse
import itertools
import json
import re
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

# MAST-resolvable names of the fields the transient line and the lens-consistency line have used.
FIELDS = {
    "abell2744": "Abell 2744",
    "abell370": "Abell 370",
    "abells1063": "Abell S1063",
    "bullet": "1E 0657-56",
    "elgordo": "ACT-CL J0102-4915",
    "macs0416": "MACS J0416.1-2403",
    "macs0417": "MACS J0417.5-1154",
    "macs0647": "MACS J0647.7+7015",
    "macs0717": "MACS J0717.5+3745",
    "macs1149": "MACS J1149.5+2223",
    "macs1423": "MACS J1423.8+2404",
    "macs1931": "MACS J1931.8-2635",
    "ms0451": "MS 0451.6-0305",
    "plckg165": "PLCK G165.7+67.0",
    "rxcj0018": "RXC J0018.5+1626",
    "rxj2129": "RX J2129.7+0005",
    "smacs0723": "SMACS J0723.3-7327",
    "spt0615": "SPT-CL J0615-5746",
    "sunrise": "WHL J013719.8-082841",
}

EPOCH_RE = re.compile(r"jw(\d{5})-?o?(\d{3})")
MIN_DT_DAYS = 30.0


def epoch_key(obs_id: str) -> tuple[str, str] | None:
    """(proposal, observation) of a level-3 obs_id (``jw02282-o010_t001_nircam_clear-f200w``)."""
    m = EPOCH_RE.match(obs_id)
    return (m.group(1), m.group(2)) if m else None


def filter_name(filters: str) -> str:
    """First science filter of a MAST ``filters`` value (``CLEAR;F200W`` -> ``F200W``)."""
    parts = [p.strip().upper() for p in re.split(r"[;,]", filters) if p.strip()]
    parts = [p for p in parts if p not in {"CLEAR", "GR150R", "GR150C"}]
    return parts[0] if parts else ""


def build_epochs(rows: list[dict]) -> dict[tuple[str, str], dict]:
    """Group rows into epochs: filters, earliest t_min (MJD) and latest t_obs_release."""
    epochs: dict[tuple[str, str], dict] = {}
    for r in rows:
        key = epoch_key(r["obs_id"])
        if key is None:
            continue
        e = epochs.setdefault(
            key, {"filters": set(), "t_min": r["t_min"], "release": r["t_obs_release"]}
        )
        if f := filter_name(r["filters"]):
            e["filters"].add(f)
        e["t_min"] = min(e["t_min"], r["t_min"])
        e["release"] = max(e["release"], r["t_obs_release"])
    return epochs


def new_pairs(
    epochs: dict[tuple[str, str], dict], used: set[tuple[str, str]], now_mjd: float
) -> list[dict]:
    """Pairs passing the frozen rule (a)-(d) of the plan, latest epoch first."""
    public = {k: e for k, e in epochs.items() if e["release"] <= now_mjd}
    out = []
    for (ka, a), (kb, b) in itertools.combinations(
        sorted(public.items(), key=lambda kv: kv[1]["t_min"]), 2
    ):
        shared = a["filters"] & b["filters"]
        dt = b["t_min"] - a["t_min"]
        if not shared or dt < MIN_DT_DAYS or (ka in used and kb in used):
            continue
        out.append(
            {
                "early": "-".join(ka),
                "late": "-".join(kb),
                "t_early_mjd": round(a["t_min"], 2),
                "t_late_mjd": round(b["t_min"], 2),
                "dt_days": round(dt, 1),
                "shared_filters": sorted(shared),
                "early_used": ka in used,
                "late_used": kb in used,
            }
        )
    return sorted(out, key=lambda p: (-p["t_late_mjd"], -p["dt_days"]))


def used_epochs(root: Path) -> set[tuple[str, str]]:
    """Every (proposal, observation) named in the D-027 transient results under ``root``."""
    used = set()
    for path in root.rglob("*"):
        if path.is_file() and path.suffix in {".md", ".json", ".ecsv", ".csv", ".txt"}:
            used.update(EPOCH_RE.findall(path.read_text(errors="ignore")))
    return used


def query_field(name: str, radius: str) -> list[dict]:
    from astroquery.mast import Observations

    tab = Observations.query_criteria(
        objectname=name,
        radius=radius,
        obs_collection="JWST",
        instrument_name="NIRCAM/IMAGE",
        calib_level=3,
        dataproduct_type="image",
    )
    cols = ("obs_id", "filters", "t_min", "t_obs_release", "dataRights")
    return [
        {c: (str(r[c]) if c in ("obs_id", "filters", "dataRights") else float(r[c])) for c in cols}
        for r in tab
    ]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--radius", default="3 arcmin")
    ap.add_argument("--used-root", type=Path, default=Path("results/transients"))
    args = ap.parse_args()

    from astropy.time import Time

    now = float(Time.now().mjd)
    used = used_epochs(args.used_root)
    with ThreadPoolExecutor(max_workers=8) as pool:
        rows = dict(
            zip(
                FIELDS,
                pool.map(lambda n: query_field(n, args.radius), FIELDS.values()),
                strict=True,
            )
        )

    summary: dict[str, dict] = {}
    upcoming: dict[str, list] = defaultdict(list)
    for field, rs in rows.items():
        epochs = build_epochs(rs)
        for k, e in epochs.items():
            if e["release"] > now:
                upcoming[field].append(
                    {"epoch": "-".join(k), "release_mjd": round(e["release"], 1)}
                )
        pairs = new_pairs(epochs, used, now)
        summary[field] = {
            "n_rows": len(rs),
            "epochs": {
                "-".join(k): {
                    "t_min_mjd": round(e["t_min"], 2),
                    "release_mjd": round(e["release"], 2),
                    "public": e["release"] <= now,
                    "used": k in used,
                    "filters": sorted(e["filters"]),
                }
                for k, e in sorted(epochs.items(), key=lambda kv: kv[1]["t_min"])
            },
            "new_pairs": pairs,
        }
    args.out.mkdir(parents=True, exist_ok=True)
    meta = {"queried_mjd": round(now, 4), "radius": args.radius, "n_used_epochs": len(used)}
    (args.out / "inventory.json").write_text(
        json.dumps({"meta": meta, "fields": summary}, indent=1) + "\n"
    )
    (args.out / "upcoming.json").write_text(json.dumps(upcoming, indent=1) + "\n")
    for field, s in summary.items():
        unused = [k for k, e in s["epochs"].items() if e["public"] and not e["used"]]
        n_ep, n_pairs = len(s["epochs"]), len(s["new_pairs"])
        print(f"{field:11s} epochs {n_ep:3d} unused-public {len(unused):3d} new-pairs {n_pairs:4d}")


if __name__ == "__main__":
    main()
