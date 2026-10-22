"""MACS0417 family 10: source-redshift scan in the CANUCS v1 maps (is 1115938 the third image?).

``../macs0417_canucsref`` flagged a family-10 predicted image as ``absent`` at the map-fitted
z 2.15; CANUCS DR1 1115938 lies 1.49" away with photo-z 1.81. At each z: the source is the
magnification-weighted least-squares position of the back-traced 10.1 and 10.2
(``a2744_s700.solve``), images by ``lensmodel.find_images``, then a one-to-one assignment of
10.1, 10.2 and 1115938 to distinct predicted images (minimum total separation,
``a2744_s700.assign``).
Frozen plan: ``results/lens_consistency/macs0417_f10_zscan/PLAN.md``.

Everything computed is ``model_prediction``; image and source positions are ``observed``.

    python scripts/m0417_f10_zscan.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import a2744_s700 as s7
import lens_consistency as lc
import m0417_frame as mf
import numpy as np
from astropy.table import Table

from jwst_anomaly import lensmodel

OUT = Path("results/lens_consistency/macs0417_f10_zscan")
FRAME = Path("results/lens_consistency/macs0417_frame")
MODEL = "macs0417-canucs-maps"
FAMILY = "10"
# CANUCS DR1 photometry catalogue (hlsp_canucs_jwst-hst_multi_macs0417-clu_multi_v1_photometry-cat,
# sha256 ca6b48aa...90dd, read 2026-10-22): RA, DEC of SOURCE 1115938 (observed)
THIRD = {"id": "CANUCS 1115938", "ra": 64.38501841936997, "dec": -11.915119106841288}
CAT_SHA256 = "ca6b48aa6404444d2e0b8965a35376ff032917071ea64cbf6306d06ee57d90dd"
# Frozen rules (PLAN.md; ASSUMPTIONs)
Z_GRID = s7.Z_GRID  # 1.000-4.000, step 0.025
MATCH_ARCSEC = 1.0
Z_KNOWN = 2.15  # CANUCS map-fitted z of family 10 (``../macs0417_images``)
PHOTOZ = (1.72, 1.88)  # z160-z840 of 1115938 (``../macs0417_canucsref/posthoc.json``)
GRID_STEP = 0.1  # arcsec, as ``m0417_extra``


def lands(sep) -> bool:
    """10.1, 10.2 and 1115938 each have a distinct predicted image within ``MATCH_ARCSEC``.
    ``sep`` is the assignment's separation per row (empty or NaN when too few images)."""
    s = np.asarray(sep, float)
    return s.size == 3 and bool(np.all(np.isfinite(s))) and bool(np.all(s <= MATCH_ARCSEC))


def known_ok(sep_pair) -> bool:
    """Known-case check: 10.1 and 10.2 each matched to a distinct predicted image within 1.0"."""
    s = np.asarray(sep_pair, float)
    return s.size == 2 and bool(np.all(np.isfinite(s))) and bool(np.all(s <= MATCH_ARCSEC))


def outcome(known: bool, rows: list[dict]) -> dict:
    """Frozen outcome from the known-case check and the per-z ``lands`` flags."""
    if not known:
        return {"outcome": "not read: known-case check failed"}
    good = [r["z"] for r in rows if r["lands"]]
    in_pz = [z for z in good if PHOTOZ[0] <= z <= PHOTOZ[1]]
    finite = [r for r in rows if r["sep_third_arcsec"] is not None]
    best = min(finite, key=lambda r: r["sep_third_arcsec"]) if finite else None
    if in_pz:
        o = "A: z artefact, counter-image identified"
    elif good:
        o = "A': lands only outside the photo-z interval"
    else:
        o = "B: never lands"
    return {
        "outcome": o,
        "z_lands": [min(good), max(good)] if good else None,
        "n_z_lands": len(good),
        "z_lands_in_photoz": [min(in_pz), max(in_pz)] if in_pz else None,
        "min_sep_third": (
            {"z": best["z"], "sep_arcsec": best["sep_third_arcsec"]} if best else None
        ),
    }


def family_images() -> Table:
    """10.1 and 10.2 on the JWST frame (``../macs0417_frame`` offset, as ``m0417_arcref``)."""
    off = json.loads((FRAME / "summary.json").read_text())["measurement"]["offset"]
    t = mf.shifted_loader(off["dra"], off["ddec"])("macs0417")
    t = t[np.asarray(t["system"], str) == FAMILY]
    t.sort("image_id")
    return t


def row(model, grid, pair: Table, three: Table, z: float) -> dict:
    pred, beta = s7.solve(model, pair, grid, z)
    idx, sep = s7.assign(pred, three)
    pidx, psep = s7.assign(pred, pair)
    near = s7.nearest(pred, three)
    mu = [s7._r(float(pred["magnification"][k]), 2) for k in idx]
    return {
        "z": float(z),
        "n_images": len(pred),
        "sep_arcsec": [s7._r(float(s)) for s in sep],
        "sep_third_arcsec": s7._r(float(sep[2])) if idx else None,
        "mu_signed": mu,
        "pair_sep_arcsec": [s7._r(float(s)) for s in psep],
        "nearest_arcsec": [s7._r(float(s)) for s in near],
        "source_spread_arcsec": s7._r(float(np.ptp(beta, axis=0).max())),
        "known_ok": known_ok(psep) if pidx else False,
        "lands": lands(sep) if idx else False,
        "predicted": [[s7._r(float(r["ra"]), 6), s7._r(float(r["dec"]), 6)] for r in pred],
    }


def main(argv: list[str] | None = None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    model, _, _ = lc.load_model(MODEL)
    lc.apply_frame_offset(MODEL, model)
    pair = family_images()
    three = Table(
        {
            "image_id": [*map(str, pair["image_id"]), THIRD["id"]],
            "ra": [*map(float, pair["ra"]), THIRD["ra"]],
            "dec": [*map(float, pair["dec"]), THIRD["dec"]],
        }
    )
    half = lc.grid_half_width(model, three)
    grid = lensmodel.DeflectionGrid.cached(
        model, lc.grid_cache_path(model, half, GRID_STEP), half, GRID_STEP
    )
    known = row(model, grid, pair, three, Z_KNOWN)
    rows = [row(model, grid, pair, three, float(z)) for z in Z_GRID]
    out = {
        "model": MODEL,
        "model_sha256": model.sha256,
        "family": FAMILY,
        "images": [*map(str, pair["image_id"])],
        "images_radec_jwst_frame": [[float(r["ra"]), float(r["dec"])] for r in pair],
        "third": THIRD,
        "canucs_cat_sha256": CAT_SHA256,
        "rules": {
            "match_arcsec": MATCH_ARCSEC,
            "z_known": Z_KNOWN,
            "photoz": list(PHOTOZ),
            "z_grid": [float(Z_GRID[0]), float(Z_GRID[-1]), 0.025],
            "grid_step_arcsec": GRID_STEP,
        },
        "known_case": known,
        **outcome(known["known_ok"], rows),
        "provenance": "model_prediction (image and source positions observed)",
        "rows": rows,
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "summary.json").write_text(json.dumps(out, indent=1, allow_nan=False) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1))
    return out


if __name__ == "__main__":
    main()
