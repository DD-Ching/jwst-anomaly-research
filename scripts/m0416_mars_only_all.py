"""MACS0416: MARS-only images in every CANUCS family (follow-up of D-034 addendum as)

Frozen plan: ``results/lens_consistency/m0416_mars_only_all/PLAN.md``. Runs the addendum-as stage
(``m0416_mars_only_phot``) over every CANUCS family: the families with a CANUCS extra prediction
reproduce addendum as (gate R), the other families give the new MARS-only targets.

    python scripts/m0416_mars_only_all.py --zip DIR --cat <photometry-cat.fits.gz> \
        --image <f200w bgsub-sci.fits.gz>
"""

from __future__ import annotations

import argparse
import json
import zipfile
from pathlib import Path

import lens_consistency as lc
import m0416_canucsref as mc
import m0416_mars_extra as me
import m0416_mars_only_phot as mp
import m0417_arcref as ar
import m0417_canucsref as cr
import m1423_canucsref as m14
import mars_hff_heldout as mo
import mars_hff_inputs as mh
import numpy as np
import requests
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from m0417_canucsref_posthoc import _clean

from jwst_anomaly import lensmodel
from jwst_anomaly.provenance import file_sha256

OUT = Path("results/lens_consistency/m0416_mars_only_all")
KNOWN = ("D115a", 64.044742, -24.060901)  # addendum as target (its summary.json)
KNOWN_TOL_ARCSEC = 0.1  # PLAN.md rule 3
SHEET_MAX = 12  # PLAN.md stop rules


def split_families(t: Table, pred: Table) -> tuple[list[str], list[str]]:
    """PLAN.md rule 1: systems whose images all back-trace finitely, split into those with a CANUCS
    extra prediction (the addendum-as set) and the others (new), in catalogue order."""
    sysv = np.asarray(t["system"], str)
    bad = ~(
        np.isfinite(np.asarray(t["beta_x"], float)) & np.isfinite(np.asarray(t["beta_y"], float))
    )
    done = set(np.asarray(pred["system"], str))
    ok = [s for s in dict.fromkeys(sysv) if not bad[sysv == s].any()]
    return [s for s in ok if s in done], [s for s in ok if s not in done]


def gate_r(targets: list[dict], known=KNOWN, tol: float = KNOWN_TOL_ARCSEC) -> bool:
    """PLAN.md rule 3: the addendum-as families give exactly one target, in ``known[0]`` within
    ``tol`` of the stored position."""
    if len(targets) != 1:
        return False
    tg = targets[0]
    d = np.hypot((tg["ra"] - known[1]) * np.cos(np.radians(known[2])), tg["dec"] - known[2]) * 3600
    return tg["system"] == known[0] and bool(d <= tol)


def sheet_rows(rows: list[dict], n: int = SHEET_MAX) -> list[dict]:
    """PLAN.md stop rules: targets with a reading, at most ``n`` by largest F_pred / sigma_loc."""
    r = [x for x in rows if x["reading"] != "no_reading"]
    key = [float(x["pred_over_sigma"]) if np.isfinite(x["pred_over_sigma"]) else -np.inf for x in r]
    return [r[i] for i in np.argsort(key, kind="stable")[::-1][:n]]


def run(zdir: Path, cat_path: Path, image: Path, out: Path = OUT) -> dict:
    if file_sha256(cat_path) != mc.CAT_SHA or file_sha256(image) != mc.IMAGE_SHA:
        raise ValueError("--cat / --image are not the pinned CANUCS DR1 MACS0416 files")
    meta = mh.fetch_zip(requests.Session(), zdir)
    folder, tag, z_lens = mh.FIELDS[me.FIELD]
    base = f"{mh.ROOT}{folder}/deflection_angle_{tag}_"
    with zipfile.ZipFile(zdir / mh.ZIP) as zf:
        mcat = mh.parse_catalogue(zf.read(mh.CATALOGUE).decode())
        (ax, hx), (ay, _) = (mh._fits(zf, base + f"{c}_w_header.fits") for c in "xy")
    model = mh.build_model(
        ax, ay, hx, False, z_lens, "see mars_hff_inputs", f"mars-{tag} (A_arcsec)"
    )
    mars = mcat[np.asarray(mcat["field"]) == me.FIELD]
    canucs = Table.read(mo.FIELDS[me.FIELD][0])
    off = mo.frame_offset(canucs, mars)
    t = mo.split_seen(mo.shift(canucs, off), mars)
    pred = mo.shift(Table.read(me.PREDICTED), off)
    grid = lensmodel.DeflectionGrid.compute(model, lc.grid_half_width(model, t), mh.GRID_STEP)
    bt = lensmodel.backtrace_images(model, t, {})
    tb = t.copy()
    tb["beta_x"], tb["beta_y"] = bt["beta_x"], bt["beta_y"]
    old, new = split_families(tb, pred)
    fams_old, tg_old = mp.solve(model, grid, t, pred, off, old)
    r_ok = gate_r(tg_old)
    summary = {
        "plan": str(out / "PLAN.md"),
        "zip": meta,
        "canucs_cat_sha256": mc.CAT_SHA,
        "image_sha256": mc.IMAGE_SHA,
        "frame_offset": off,
        "n_families_old": len(old),
        "n_families_new": len(new),
        "gate_R": {"passes": r_ok, "targets_old": tg_old},
    }
    out.mkdir(parents=True, exist_ok=True)
    if not r_ok:  # PLAN.md rule 3: stop, no reading
        (out / "summary.json").write_text(json.dumps(_clean(summary), indent=1) + "\n")
        return summary
    fams_new, targets = mp.solve(model, grid, t, pred, off, new)
    fams = {**fams_old, **fams_new}

    tab = Table.read(cat_path)
    cat = cr.to_cat(tab)
    depth = cr.depth_mag(cat["flux"], cat["snr"])
    controls, rows = mp.photometry(fams, targets, cat, depth)
    c2 = ar.gate_c2([c["class"] for c in controls], [c["arc_ref"] for c in controls])
    cn = [c for c in controls if c["system"] in fams_new]
    c2_new = ar.gate_c2([c["class"] for c in cn], [c["arc_ref"] for c in cn])
    with fits.open(image) as hl:
        w, d = WCS(hl[0].header), np.asarray(hl[0].data, float)
    src_x, src_y = (
        np.asarray(v, float) for v in w.world_to_pixel(SkyCoord(tab["RA"], tab["DEC"], unit="deg"))
    )
    rng = np.random.default_rng(m14.SEED)
    ctab = Table(rows=controls)
    m14.apply_local(ctab, d, w, src_x, src_y, rng)
    ttab = Table(rows=rows) if rows else Table()
    if rows:
        m14.apply_local(ttab, d, w, src_x, src_y, rng)
        ttab["reading"] = [mp.reading(str(c), c2["passes"]) for c in ttab["class_local"]]
    for tb_, name in ((ctab, "controls.ecsv"), (ttab, "targets.ecsv")):
        tb_.meta["provenance"] = "derived"
        tb_.meta["source"] = f"{out / 'PLAN.md'}; MARS maps (model_prediction) + CANUCS DR1 F200W"
        tb_.write(out / name, format="ascii.ecsv", overwrite=True)
    trows = [dict(zip(ttab.colnames, r, strict=True)) for r in ttab] if rows else []
    shown = sheet_rows(trows)
    mp.cutout_png(d, w, shown, fams, out / "targets.png")
    summary.update(
        {
            "depth_mag": depth,
            "search_arcsec": mp.SEARCH_ARCSEC,
            "n_targets_new": len(rows),
            "gate_C2": c2,
            "gate_C2_new_families_reported": c2_new,
            "control_classes_local": m14._count(ctab["class_local"]),
            "target_classes_local": m14._count(ttab["class_local"]) if rows else {},
            "target_readings": m14._count(ttab["reading"]) if rows else {},
            "sheet": [f"{r['system']} {r['ra']:.6f} {r['dec']:.6f}" for r in shown],
        }
    )
    (out / "summary.json").write_text(json.dumps(_clean(summary), indent=1, allow_nan=False) + "\n")
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--zip", type=Path, required=True, help="directory for the MARS zip")
    ap.add_argument("--cat", type=Path, required=True)
    ap.add_argument("--image", type=Path, required=True)
    a = ap.parse_args(argv)
    try:
        s = run(a.zip, a.cat, a.image)
    finally:
        (a.zip / mh.ZIP).unlink(missing_ok=True)
    print(json.dumps(_clean({k: v for k, v in s.items() if k != "zip"}), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
