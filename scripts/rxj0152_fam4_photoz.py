"""RX J0152 family 4: RELICS BPZ photo-z of the three chance-level matches (addendum bg follow-up).

Plan (frozen before any catalogue row was read):
``results/lens_consistency/rxj0152_fam4_photoz/PLAN.md``. Matches segments 112 / 127 / 131 of
addendum bg, family 4's own images and the 18 SLICE control images of families 1-6 to the RELICS
HST BPZ catalogue (nearest row within 0.6", no offset correction), applies gate P and the
per-segment classes, and writes ``summary.json``. Catalogue values ``observed`` (RELICS), photo-z
``derived`` (RELICS BPZ), family z ``model_prediction`` (SLICE; family 1 spectroscopic), classes
``derived``.

    python scripts/rxj0152_fam4_photoz.py \\
        --bpz <hlsp_relics_hst_acs-wfc3ir_cl0152m13_multi_v1_cat.txt>
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
from arcs_spt0615 import read_bpz  # noqa: E402

ROOT = Path(__file__).resolve().parents[1] / "results" / "lens_consistency"
OUT = ROOT / "rxj0152_fam4_photoz"
SLICE = ROOT / "slice_relics_rxj0152" / "slice_images.ecsv"
BG = ROOT / "rxj0152_fam4_photometry_v2" / "result.json"
SEGMENTS = (112, 127, 131)
Z_FAM4 = 1.69  # SLICE model z of family 4 (photometric)
CONTROL_FAMILIES = ("1", "2", "3", "4", "5", "6")
STAR_STEL = 0.8  # ASSUMPTION (plan): stellarity at or above this is a star
MATCH_ARCSEC = 0.6  # ASSUMPTION (plan)
ODDS_MIN = 0.9  # ASSUMPTION (plan): reliable BPZ peak
GATE_P_MIN_MATCHED = 8
GATE_P_MIN_FRAC = 0.6


def match(ra: float, dec: float, cat: Table) -> dict | None:
    """Nearest non-star catalogue row within MATCH_ARCSEC, or None."""
    gal = cat[np.asarray(cat["stel"], float) < STAR_STEL]
    if len(gal) == 0:
        return None
    sep = (
        SkyCoord(ra, dec, unit="deg").separation(SkyCoord(gal["RA"], gal["Dec"], unit="deg")).arcsec
    )
    k = int(np.argmin(sep))
    if sep[k] > MATCH_ARCSEC:
        return None
    r = gal[k]
    return {
        "relics_id": int(r["id"]),
        "sep_arcsec": round(float(sep[k]), 3),
        "dra_arcsec": float((r["RA"] - ra) * np.cos(np.radians(dec)) * 3600),
        "ddec_arcsec": float((r["Dec"] - dec) * 3600),
        "zb": float(r["zb"]),
        "zbmin": float(r["zbmin"]),
        "zbmax": float(r["zbmax"]),
        "odds": float(r["odds"]),
    }


def contains(m: dict, z: float) -> bool:
    return m["zbmin"] <= z <= m["zbmax"]


def classify(m: dict | None, z: float = Z_FAM4) -> str:
    """Plan step 4: consistent / excluded / undetermined."""
    if m is None:
        return "undetermined"
    if contains(m, z):
        return "consistent"
    return "excluded" if m["odds"] >= ODDS_MIN else "undetermined"


def gate_p(rows: list[dict]) -> dict:
    """Plan step 3: rows carry ``match`` (dict or None) and ``z``."""
    matched = [r for r in rows if r["match"] is not None]
    ok = [r for r in matched if contains(r["match"], r["z"])]
    frac = len(ok) / len(matched) if matched else float("nan")
    return {
        "n_images": len(rows),
        "n_matched": len(matched),
        "n_contain": len(ok),
        "frac": frac,
        "pass": len(matched) >= GATE_P_MIN_MATCHED and frac >= GATE_P_MIN_FRAC,
    }


def outcome(gate_pass: bool, classes: list[str]) -> str:
    if not gate_pass:
        return "U"
    if "consistent" in classes:
        return "C"
    if classes and all(c == "excluded" for c in classes):
        return "E"
    return "U"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bpz", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    cat = read_bpz(args.bpz)
    sl = Table.read(SLICE)
    controls = []
    for r in sl:
        if str(r["system"]) in CONTROL_FAMILIES:
            m = match(float(r["ra"]), float(r["dec"]), cat)
            controls.append({"image": str(r["image_id"]), "z": float(r["z"]), "match": m})
    gp = gate_p(controls)
    off = [(c["match"]["dra_arcsec"], c["match"]["ddec_arcsec"]) for c in controls if c["match"]]
    bg = json.loads(BG.read_text())
    pos = {int(c["label"]): (c["ra"], c["dec"]) for c in bg["region_segments"]}
    segs = []
    for lab in SEGMENTS:
        m = match(*pos[lab], cat)
        segs.append(
            {
                "segment": lab,
                "ra": pos[lab][0],
                "dec": pos[lab][1],
                "match": m,
                "class": classify(m),
            }
        )
    classes = [s["class"] for s in segs]
    summary = {
        "plan": "PLAN.md",
        "catalogue": args.bpz.name,
        "catalogue_sha256": hashlib.sha256(args.bpz.read_bytes()).hexdigest(),
        "n_rows": len(cat),
        "z_fam4": Z_FAM4,
        "gate_P": gp,
        "controls": controls,
        "control_median_offset_arcsec": (
            [float(np.median([o[0] for o in off])), float(np.median([o[1] for o in off]))]
            if off
            else None
        ),
        "family4": [c for c in controls if c["image"].startswith("4.")],
        "segments": segs,
        "outcome": outcome(gp["pass"], classes),
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({"gate_P": gp, "classes": classes, "outcome": summary["outcome"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
