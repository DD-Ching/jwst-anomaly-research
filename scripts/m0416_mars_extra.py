"""MACS0416: are the CANUCS-model extra images also predicted by the MARS maps? (two-model check)

Frozen plan: ``results/lens_consistency/m0416_mars_extra/PLAN.md``. Takes the stored extra-image
predictions of the CANUCS Lenstool model (``macs0416_canucsref/predicted_macs0416.ecsv``), moves
them and the CANUCS images into the MARS frame (``mars_hff_heldout.frame_offset``), solves every
family with the MARS maps (mean source of all its CANUCS images) and classifies each extra
prediction (``classify``). No flux or pixel reading.

    python scripts/m0416_mars_extra.py --zip DIR
"""

from __future__ import annotations

import argparse
import json
import zipfile
from pathlib import Path

import lens_consistency as lc
import mars_hff_heldout as mo
import mars_hff_inputs as mh
import numpy as np
import requests
from astropy.table import Table
from m0417_canucsref_posthoc import _clean

from jwst_anomaly import lensmodel

OUT = Path("results/lens_consistency/m0416_mars_extra")
FIELD = "MACS J0416.1-2403"
PREDICTED = Path("results/lens_consistency/macs0416_canucsref/predicted_macs0416.ecsv")
MATCH_ARCSEC = 3.0  # ASSUMPTION (PLAN.md rule 2; gate H's radius)
HIGH_MU = 10.0  # PLAN.md rule 3 (reported split)
FLAG = "D115b"


def match_observed(px, py, ox, oy, r: float = MATCH_ARCSEC) -> np.ndarray:
    """PLAN.md rule 2: mask of the predicted images matched one-to-one (greedy by distance,
    within ``r``) to observed images."""
    px, py, ox, oy = (np.asarray(v, float) for v in (px, py, ox, oy))
    used = np.zeros(len(px), bool)
    if not len(px) or not len(ox):
        return used
    d = np.hypot(ox[:, None] - px[None, :], oy[:, None] - py[None, :])
    taken_obs = np.zeros(len(ox), bool)
    for flat in np.argsort(d, axis=None):
        i, j = np.unravel_index(flat, d.shape)
        if d[i, j] > r:
            break
        if not taken_obs[i] and not used[j]:
            taken_obs[i] = used[j] = True
    return used


def classify(d_nearest: float, testable: bool, r: float = MATCH_ARCSEC) -> str:
    """PLAN.md rule 2: "both" when the nearest unmatched MARS image is within ``r``."""
    if not testable:
        return "not testable"
    return "both" if np.isfinite(d_nearest) and d_nearest <= r else "CANUCS only"


def family_rows(model, grid, t: Table, pred: Table) -> tuple[list[tuple], dict[str, int]]:
    """Rows per extra prediction (``run`` names the columns) and, per solved family, the number of
    MARS images matched to neither an observed image nor a CANUCS extra prediction."""
    bt = lensmodel.backtrace_images(model, t, {})
    sysv = np.asarray(t["system"], str)
    ox, oy = model.to_frame(np.asarray(t["ra"], float), np.asarray(t["dec"], float))
    ex, ey = model.to_frame(np.asarray(pred["ra"], float), np.asarray(pred["dec"], float))
    rows, extras = [], {}
    for s in dict.fromkeys(np.asarray(pred["system"], str)):
        fam = sysv == s
        k = np.where(np.asarray(pred["system"], str) == s)[0]
        bx, by = np.asarray(bt["beta_x"])[fam], np.asarray(bt["beta_y"])[fam]
        ok = fam.any() and np.isfinite(bx).all() and np.isfinite(by).all()
        if ok:
            p = lensmodel.find_images(
                model, grid, float(bx.mean()), float(by.mean()), float(np.asarray(t["z"])[fam][0])
            )
            px, py = np.asarray(p["x"], float), np.asarray(p["y"], float)
            free = ~match_observed(px, py, ox[fam], oy[fam])
            px, py, mu = px[free], py[free], np.abs(np.asarray(p["magnification"], float))[free]
            hit = np.zeros(len(px), bool)
        for j in k:
            d, m = np.nan, np.nan
            if ok and len(px):
                dd = np.hypot(px - ex[j], py - ey[j])
                n = int(np.argmin(dd))
                d, m = float(dd[n]), float(mu[n])
                hit[n] |= d <= MATCH_ARCSEC
            elif ok:
                d = np.inf
            rows.append(
                (
                    s,
                    float(pred["ra"][j]),
                    float(pred["dec"][j]),
                    float(pred["abs_mu"][j]),
                    str(pred["class"][j]),
                    d,
                    m,
                    classify(d, ok),
                )
            )
        if ok:
            extras[s] = int((~hit).sum())
    return rows, extras


def summarise(r: Table, extras: dict) -> dict:
    cls = np.asarray(r["class_mars"], str)
    hi = np.asarray(r["abs_mu_canucs"], float) >= HIGH_MU
    keys = ("both", "CANUCS only", "not testable")
    flag = r[(np.asarray(r["system"], str) == FLAG) & (np.asarray(r["class_canucs"]) == "absent")]
    return {
        "n_extra": len(r),
        "classes": {k: int(np.sum(cls == k)) for k in keys},
        "classes_abs_mu_ge_10": {k: int(np.sum(cls[hi] == k)) for k in keys},
        "flag_D115b": [
            {
                "class": str(x["class_mars"]),
                "d_nearest_arcsec": float(x["d_nearest_arcsec"]),
                "abs_mu_mars": float(x["abs_mu_mars"]),
            }
            for x in flag
        ],
        "mars_unmatched_extra_per_family": extras,
    }


def run(zdir: Path, out: Path = OUT) -> dict:
    meta = mh.fetch_zip(requests.Session(), zdir)
    folder, tag, z_lens = mh.FIELDS[FIELD]
    base = f"{mh.ROOT}{folder}/deflection_angle_{tag}_"
    with zipfile.ZipFile(zdir / mh.ZIP) as zf:
        cat = mh.parse_catalogue(zf.read(mh.CATALOGUE).decode())
        (ax, hx), (ay, _) = (mh._fits(zf, base + f"{c}_w_header.fits") for c in "xy")
    model = mh.build_model(
        ax, ay, hx, False, z_lens, "see mars_hff_inputs", f"mars-{tag} (A_arcsec)"
    )
    mars = cat[np.asarray(cat["field"]) == FIELD]
    canucs = Table.read(mo.FIELDS[FIELD][0])
    off = mo.frame_offset(canucs, mars)
    t = mo.split_seen(mo.shift(canucs, off), mars)
    pred = mo.shift(Table.read(PREDICTED), off)
    grid = lensmodel.DeflectionGrid.compute(model, lc.grid_half_width(model, t), mh.GRID_STEP)
    rows, extras = family_rows(model, grid, t, pred)
    names = (
        "system",
        "ra",
        "dec",
        "abs_mu_canucs",
        "class_canucs",
        "d_nearest_arcsec",
        "abs_mu_mars",
        "class_mars",
    )
    r = Table(rows=rows, names=names)
    out.mkdir(parents=True, exist_ok=True)
    r.write(out / "extra_classes.ecsv", format="ascii.ecsv", overwrite=True)
    summary = {
        "plan": str(out / "PLAN.md"),
        "zip": meta,
        "frame_offset": off,
        **summarise(r, extras),
    }
    (out / "summary.json").write_text(json.dumps(_clean(summary), indent=1, allow_nan=False) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--zip", type=Path, required=True, help="directory for the MARS zip")
    a = ap.parse_args(argv)
    try:
        s = run(a.zip)
    finally:
        (a.zip / mh.ZIP).unlink(missing_ok=True)
    print(json.dumps(_clean({k: v for k, v in s.items() if k != "zip"}), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
