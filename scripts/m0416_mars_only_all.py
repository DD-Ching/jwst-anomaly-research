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
NEAR_ARCSEC = 6.0  # post hoc listing radius (not a rule)


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
    return [r[i] for i in np.argsort(-np.asarray(key), kind="stable")[:n]]


def _with_beta(model, t: Table) -> Table:
    """``t`` with its back-traced source columns (input of ``split_families``)."""
    bt = lensmodel.backtrace_images(model, t, {})
    tb = t.copy()
    tb["beta_x"], tb["beta_y"] = bt["beta_x"], bt["beta_y"]
    return tb


def critical_distance(model, x: float, y: float, z: float, half=6.0, step=0.1) -> float:
    """Distance (arcsec) from (x, y) to the nearest grid point where det A is finite and has the
    opposite sign (NaN when det A at (x, y) is not finite, or when none is found and part of the
    window is not finite)."""
    g = np.arange(-half, half + step / 2, step)
    gx, gy = np.meshgrid(x + g, y + g)
    a = model.lens_map(np.array([x, *gx.ravel()]), np.array([y, *gy.ravel()]), z)
    det = a["a11"] * a["a22"] - a["a12"] * a["a21"]
    if not np.isfinite(det[0]):
        return float("nan")
    if det[0] == 0:
        return 0.0
    other = np.isfinite(det[1:]) & (det[1:] * det[0] < 0)
    if not other.any():  # unknown when part of the window is off the map, else none within it
        return float("nan") if not np.isfinite(det[1:]).all() else float("inf")
    return float(np.hypot(gx.ravel() - x, gy.ravel() - y)[other].min())


def posthoc(zdir: Path, out: Path = OUT) -> dict:
    """Post hoc (not a frozen rule), per ``lensing_residual_input`` target: the CANUCS extra
    predictions of its family, the MARS critical-curve distance and the family's MARS images."""
    _, model, t, pred, off, grid = mp.load(zdir)
    psys = np.asarray(pred["system"], str)
    cm, _, _ = lc.load_model(mc.MODEL)  # the pinned CANUCS Lenstool model (D-044)
    lc.apply_frame_offset(mc.MODEL, cm)
    canucs = Table.read(mo.FIELDS[me.FIELD][0])  # observed images, CANUCS frame (targets' frame)
    canucs["system"] = [lensmodel.image_family(str(i)) for i in canucs["image_id"]]
    solved = set(np.asarray(t["image_id"], str))  # images offered to the pairing (``split_seen``)
    old, new = split_families(_with_beta(model, t), pred)
    all_fams = list(dict.fromkeys(np.asarray(canucs["system"], str)))
    cgrid = lensmodel.DeflectionGrid.compute(cm, lc.grid_half_width(cm, canucs), mc.GRID_STEP_K)
    res = []
    for tg in Table.read(out / "targets.ecsv"):
        if tg["reading"] != "lensing_residual_input":
            continue
        s = str(tg["system"])
        z = float(np.asarray(t["z"])[np.asarray(t["system"], str) == s][0])
        fams, tgs = mp.solve(model, grid, t, pred, off, [s])
        dec = float(tg["dec"]) + off["dy_arcsec"] / 3600.0  # ``mp.unshift`` inverse
        ra = float(tg["ra"]) + off["dx_arcsec"] / 3600.0 / np.cos(np.deg2rad(dec))
        x, y = (float(v) for v in model.to_frame(ra, dec))
        fam = np.asarray(t["system"], str) == s
        bt = lensmodel.backtrace_images(model, t[fam], {})
        per_image = []
        for i, (bx, by) in enumerate(zip(bt["beta_x"], bt["beta_y"], strict=True)):
            p = lensmodel.find_images(model, grid, float(bx), float(by), z)
            d = np.hypot(np.asarray(p["x"], float) - x, np.asarray(p["y"], float) - y)
            near = d <= NEAR_ARCSEC
            per_image.append(
                {
                    "source_from": str(t["image_id"][fam][i]),
                    "n_images": len(p),
                    "nearest_to_target_arcsec": float(d.min()) if len(d) else float("inf"),
                    "images_within_6arcsec": [
                        {"d_arcsec": float(dd), "abs_mu": float(abs(mm))}
                        for dd, mm in zip(
                            d[near], np.asarray(p["magnification"])[near], strict=True
                        )
                    ],
                }
            )
        cx, cy = (float(v) for v in cm.to_frame(float(tg["ra"]), float(tg["dec"])))
        cfam = canucs[np.asarray(canucs["system"], str) == s]
        cbt = lensmodel.backtrace_images(cm, cfam, {})
        canucs_per_image = []
        for i, (bx, by) in enumerate(zip(cbt["beta_x"], cbt["beta_y"], strict=True)):
            p = lensmodel.find_images(cm, cgrid, float(bx), float(by), z)
            d = np.hypot(np.asarray(p["x"], float) - cx, np.asarray(p["y"], float) - cy)
            canucs_per_image.append(
                {
                    "source_from": str(cfam["image_id"][i]),
                    "n_images": len(p),
                    "nearest_to_target_arcsec": float(d.min()) if len(d) else float("inf"),
                }
            )
        res.append(
            {
                "system": s,
                "canucs_critical_distance_arcsec": critical_distance(cm, cx, cy, z, half=10.0),
                "canucs_source_per_image": canucs_per_image,
                "source_per_image": per_image,
                "z": z,
                "abs_mu": float(tg["abs_mu"]),
                "n_canucs_extra_predictions": int((psys == s).sum()),
                "mars_critical_distance_arcsec": critical_distance(model, x, y, z),
                "observed_images": [
                    {
                        "id": str(r["image_id"]),
                        "paired_to_mars": str(r["image_id"])
                        in {im["id"] for im in fams.get(s, [])},
                        "d_to_target_arcsec": float(
                            np.hypot(
                                (r["ra"] - tg["ra"]) * np.cos(np.radians(tg["dec"])),
                                r["dec"] - tg["dec"],
                            )
                            * 3600
                        ),
                    }
                    for r in cfam
                    if str(r["image_id"]) in solved
                ],
                "mars_images_paired": [
                    {"id": im["id"], "abs_mu": im["mu"]} for im in fams.get(s, [])
                ],
                "mars_only_images": [
                    {
                        "abs_mu": o["mu"],
                        "d_to_target_arcsec": float(
                            np.hypot(
                                (o["ra"] - tg["ra"]) * np.cos(np.radians(tg["dec"])),
                                o["dec"] - tg["dec"],
                            )
                            * 3600
                        ),
                    }
                    for o in tgs
                ],
            }
        )
    p = {
        "label": "post hoc, not a frozen rule",
        "map_pixel_arcsec": float(model.pixel_arcsec),
        "near_arcsec": NEAR_ARCSEC,
        "canucs_grid_step_arcsec": mc.GRID_STEP_K,
        "families": {
            "total": len(all_fams),
            "solved": len(old) + len(new),
            "dropped": [f for f in all_fams if f not in set(old) | set(new)],
        },
        "targets": res,
    }
    (out / "posthoc.json").write_text(json.dumps(_clean(p), indent=1, allow_nan=False) + "\n")
    return p


def run(zdir: Path, cat_path: Path, image: Path, out: Path = OUT) -> dict:
    if file_sha256(cat_path) != mc.CAT_SHA or file_sha256(image) != mc.IMAGE_SHA:
        raise ValueError("--cat / --image are not the pinned CANUCS DR1 MACS0416 files")
    meta, model, t, pred, off, grid = mp.load(zdir)
    old, new = split_families(_with_beta(model, t), pred)
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
    for name in ("controls.ecsv", "targets.ecsv", "targets.png", "posthoc.json"):
        (out / name).unlink(missing_ok=True)  # no stale outputs from an earlier run
    if not r_ok:  # PLAN.md rule 3: stop, no reading
        txt = json.dumps(_clean(summary), indent=1, allow_nan=False)
        (out / "summary.json").write_text(txt + "\n")
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
    ap.add_argument("--cat", type=Path)
    ap.add_argument("--image", type=Path)
    ap.add_argument("--posthoc", action="store_true", help="only the post hoc table (after a run)")
    a = ap.parse_args(argv)
    if not a.posthoc and (a.cat is None or a.image is None):
        ap.error("--cat and --image are required without --posthoc")
    try:
        if a.posthoc:
            print(json.dumps(posthoc(a.zip), indent=1))
            return 0
        s = run(a.zip, a.cat, a.image)
    finally:
        (a.zip / mh.ZIP).unlink(missing_ok=True)
    print(json.dumps(_clean({k: v for k, v in s.items() if k != "zip"}), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
