"""MACS J1931 stage 3c: the stage 3 / 3b route with one builder change (background box 64 -> 256
pixels), re-validated on MACS0417 first, and a chance-match control added to frame check F.

Frozen plan: ``results/lens_consistency/macs1931_stage3c/PLAN.md``.

    python scripts/m1931_stage3c.py gate              # rule 3c-1: MACS0417 V2 chain, new box
    python scripts/m1931_stage3c.py run --maps DIR    # rules 3c-2 / 3c-3 (after the gate passes)
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import astropy.units as u
import m0417_canucsref_posthoc as cp
import m1931_deblend as md
import m1931_stage3 as s3
import numpy as np
from astropy.coordinates import SkyCoord

BKG_BOX_3C = 256  # ASSUMPTION (rule 3c-0): 8.0" at 0.0312"/px, longer than the unsegmented arcs
OUT_GATE = Path("results/lens_consistency/macs0417_ownref_3c")
OUT = Path("results/lens_consistency/macs1931_stage3c")
CHANCE_SHIFT_ARCSEC = 5.0  # ASSUMPTION (rule 3c-2)
CHANCE_N_DIRECTIONS = 8  # N, NE, E, ... (position angles 0, 45, ..., 315 deg)
CHANCE_MIN_RATIO = 2.0  # ASSUMPTION (rule 3c-2): matches >= 2 x mean shifted matches


def shifted_counts(ra, dec, cra, cdec, match_arcsec: float = s3.be.F_MATCH_ARCSEC) -> list[int]:
    """Rule 3c-2: images matched within ``match_arcsec`` after moving every image position by
    ``CHANCE_SHIFT_ARCSEC`` along each of ``CHANCE_N_DIRECTIONS`` position angles."""
    g = SkyCoord(np.asarray(ra, float), np.asarray(dec, float), unit="deg")
    c = SkyCoord(np.asarray(cra, float), np.asarray(cdec, float), unit="deg")
    out = []
    for k in range(CHANCE_N_DIRECTIONS):
        pa = 360.0 * k / CHANCE_N_DIRECTIONS
        moved = g.directional_offset_by(pa * u.deg, CHANCE_SHIFT_ARCSEC * u.arcsec)
        sep = moved.match_to_catalog_sky(c)[1].arcsec
        out.append(int(np.sum(sep <= match_arcsec)))
    return out


def chance_passes(n_matched: int, shifted: list[int]) -> bool:
    """Rule 3c-2: ``n_matched >= CHANCE_MIN_RATIO * mean(shifted)``."""
    return bool(n_matched >= CHANCE_MIN_RATIO * float(np.mean(shifted)))


def add_chance(f: dict, images, cat) -> None:
    """Rule 3c-2 in place on the rule 5 dict: F passes only if the stage 3 F passes and the
    chance control passes; the stage 3 outcome is kept as ``passes_stage3_rule``."""
    sh = shifted_counts(images["ra"], images["dec"], cat["ra"], cat["dec"])
    f["passes_stage3_rule"] = f["passes"]
    f["chance_shifted_counts"] = sh
    f["chance_mean"] = float(np.mean(sh))
    f["chance_passes"] = chance_passes(f["n_matched"], sh)
    f["chance_rules"] = {
        "shift_arcsec": CHANCE_SHIFT_ARCSEC,
        "n_directions": CHANCE_N_DIRECTIONS,
        "min_ratio": CHANCE_MIN_RATIO,
    }
    f["passes"] = bool(f["passes_stage3_rule"] and f["chance_passes"])


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("step", choices=["gate", "run", "posthoc"])
    ap.add_argument("--maps", type=Path)
    ap.add_argument("--cache", type=Path, default=Path("data/cache/m1931s3"))
    a = ap.parse_args(argv)
    if a.step == "gate":
        return md.stage_v2(OUT_GATE, BKG_BOX_3C)
    OUT.mkdir(parents=True, exist_ok=True)
    if a.step == "posthoc":
        s = cp._clean(s3.posthoc(a.cache, OUT, True, BKG_BOX_3C))
        (OUT / "posthoc.json").write_text(json.dumps(s, indent=1, allow_nan=False) + "\n")
        print(json.dumps(s, indent=1))
        return 0
    gate = json.loads((OUT_GATE / "summary.json").read_text())
    if not gate["stage_V2_passes"]:
        ap.error("rule 3c-1 gate failed: stage 3c stops (PLAN.md)")
    if a.maps is None:
        ap.error("--maps is required for run")
    s = s3.run(a.maps, a.cache, OUT, True, BKG_BOX_3C, add_chance)
    s["plan"] = str(OUT / "PLAN.md")
    s["rule3c_1_gate"] = {k: gate[k] for k in ("gate_C2_shared", "stage_V2_passes")}
    s = cp._clean(s)
    (OUT / "summary.json").write_text(json.dumps(s, indent=1, allow_nan=False) + "\n")
    keys = ("rule2_gate_K_wls", "rule5_frame_F", "stopped")
    print(json.dumps({k: s[k] for k in keys if k in s}, indent=1))
    if "frozen" in s:
        print(json.dumps({k: s["frozen"][k] for k in ("gate_C2_shared", "flags")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
