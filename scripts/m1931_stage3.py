"""MACS J1931: extra-image test with the VENUS v2.1 maps (WLS source) and the own deblended F200W
reference catalogue (``m1931_deblend.build``, accepted by stage V2 on MACS0417).

Frozen plan: ``results/lens_consistency/macs1931_counter/PLAN.md`` ("Stage 3: frozen rules").
Rule 1 loads the md5-verified maps and the spec-z images; rule 2 repeats gate K-WLS; rule 3
writes the frame tables (``frame/``); rule 4 builds the reference on an S3 section of the GO 6882
F200W mosaic; rule 5 is frame check F; rules 6-9 are the V2-b chain, the footprint rule, gate C2
and the flags.

    python scripts/m1931_stage3.py --maps DIR [--cache data/cache/m1931s3]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import a2744_uncover as au
import bullet_extra as be
import lens_consistency as lc
import m0417_arcref as ar
import m0417_canucsref_posthoc as cp
import m0417_extra as mx
import m0417_images as mi
import m1423_canucsref as m14
import m1931_deblend as md
import macs1931_gatek as mg
import numpy as np
from astropy.table import Table, vstack

from jwst_anomaly import lensmodel

FIELD = "macs1931"
OUT = Path("results/lens_consistency/macs1931_stage3")
FRAME = OUT / "frame"
GATEK_SUMMARY = mg.OUT / "summary.json"
S3_6882 = "s3://stpubdata/jwst/public/jw06882/L3/t/o038/jw06882-o038_t041_nircam_clear-f200w_{0}"
F_MIN_N = 20  # frozen rule 5 (ASSUMPTION; 51 images)
OFF_IMAGE = "off_image"
POSTHOC_RADII = (0.5, 1.0, 1.5)  # post hoc diagnostic only


def footprint_class(cls: str, x: float, y: float, data: np.ndarray) -> str:
    """Rule 7: ``off_image`` when (x, y) lies outside ``data`` or its nearest pixel is NaN,
    else ``cls`` unchanged."""
    if not (np.isfinite(x) and np.isfinite(y)):
        return OFF_IMAGE
    ix, iy = int(np.rint(x)), int(np.rint(y))
    ny, nx = data.shape
    if not (0 <= ix < nx and 0 <= iy < ny) or not np.isfinite(data[iy, ix]):
        return OFF_IMAGE
    return cls


def apply_footprint(tab: Table, data: np.ndarray, w) -> None:
    """Rule 7 in place on ``class_local``; the earlier class is kept in ``class_local_pre_fp``."""
    tab["class_local_pre_fp"] = np.asarray(tab["class_local"], str)
    x, y = (np.asarray(v, float) for v in w.world_to_pixel_values(tab["ra"], tab["dec"]))
    tab["class_local"] = [
        footprint_class(str(c), xi, yi, data)
        for c, xi, yi in zip(tab["class_local_pre_fp"], x, y, strict=True)
    ]


def frame_tables(model, spec: Table, cache: Path) -> tuple[dict, object]:
    """Rules 2-3: gate K-WLS repeated, then predicted / control tables from the WLS source."""
    half = lc.grid_half_width(model, spec)
    grid = lensmodel.DeflectionGrid.compute(model, half, mg.GRID_STEP)
    bt = lensmodel.backtrace_images(model, spec, {})
    k = mi.gate_k(au.imageplane(model, grid, bt))
    stored = json.loads(GATEK_SUMMARY.read_text())["reported"]["gate_K_numbers_wls"]
    res = {"rule2_gate_K_wls": k, "rule2_stored_passes": stored["passes"]}
    if k["passes"] != stored["passes"] or not k["passes"]:
        res["stopped"] = "gate K-WLS outcome differs from the stored one or fails (rule 2)"
        return res, None
    wbt = au.wls_backtrace(model, bt)
    shapes = lc.load_shapes(cache / Path(S3_6882.format("cat.ecsv")).name)
    depth = lc.depth_mag(shapes)
    table, unpredicted = lc.predict_counter_images(model, grid, wbt, shapes, depth=depth)
    control = table[np.asarray(table["image_class"], str) == "observed"]
    control["image_class"] = "control"
    FRAME.mkdir(parents=True, exist_ok=True)
    lc._write(table, FRAME / f"predicted_{FIELD}.ecsv")
    lc._write(control, FRAME / f"control_{FIELD}.ecsv")
    cls = np.asarray(table["image_class"], str)
    frame = {
        "unpredicted": unpredicted,
        "pipeline_depth_mag": depth,
        "image_classes": {c: int((cls == c).sum()) for c in sorted(set(cls))},
        "grid_half_width_arcsec": half,
        "grid_step_arcsec": mg.GRID_STEP,
    }
    (FRAME / "summary.json").write_text(
        json.dumps(cp._clean({FIELD: frame}), indent=1, allow_nan=False) + "\n"
    )
    res["rule3_frame_tables"] = frame
    return res, table


def run(maps: Path, cache: Path) -> dict:
    import s3fs

    shas, model_sha = mg.verified_inputs(maps)  # rule 1
    model = mg.load_models(maps, model_sha, ("P_pixels",))["P_pixels"]
    images = mg.load_images()
    spec = images[np.asarray(images["z_spec"], bool)]
    summary: dict = {"plan": "results/lens_consistency/macs1931_counter/PLAN.md (Stage 3)"}
    summary["inputs_sha256"] = shas
    summary["n_images_used"], summary["n_images_spec"] = len(images), len(spec)
    res, table = frame_tables(model, spec, cache)
    summary.update(res)
    if table is None:
        return summary
    ctrl = Table.read(FRAME / f"control_{FIELD}.ecsv")
    pos = vstack([table["ra", "dec"], ctrl["ra", "dec"]])
    fs = s3fs.S3FileSystem(anon=True, default_block_size=64 * 1024)
    uri = S3_6882.format("i2d.fits")
    sci, err, w, box = md.read_section(fs, uri[5:], pos["ra"], pos["dec"])  # rule 4
    cat, meta, _, data, (sx, sy) = md.build(sci, err, w)
    depth = md.kron_depth(cat)
    f = be.frame_check_f(images["ra"], images["dec"], cat["ra"], cat["dec"])  # rule 5
    f["rules"]["min_n"] = F_MIN_N
    f["passes"] = bool(
        f["n_matched"] >= F_MIN_N
        and abs(f["median_dx_arcsec"]) <= be.F_MAX_MEDIAN_ARCSEC
        and abs(f["median_dy_arcsec"]) <= be.F_MAX_MEDIAN_ARCSEC
    )
    summary["rule4_catalogue"] = {"image": uri, "section_y0_y1_x0_x1": [int(v) for v in box]}
    summary["rule4_catalogue"].update(meta, depth_mag_v2a=depth)
    summary["rule5_frame_F"] = f
    if not f["passes"]:
        summary["stopped"] = "frame check F failed (rule 5)"
        return summary
    mx.load_images = lambda field: spec  # ar.run's image list (m0417_frame.shifted_loader)
    inputs = {"image": uri, "section_y0_y1_x0_x1": [int(v) for v in box], **meta}
    for contain, out in ((False, OUT / "without_contain"), (True, OUT)):  # rule 10 / rules 6-9
        out.mkdir(parents=True, exist_ok=True)
        s = ar.run(
            Path("."),
            cat=cat,
            depth=depth,
            cat_inputs=inputs,
            out=out,
            field=FIELD,
            frame=FRAME,
            offset=(0.0, 0.0),
            contain=contain,
        )
        summary["frozen" if contain else "without_contain"] = chain(s, out, data, w, sx, sy, meta)
    return summary


def chain(s: dict, out: Path, data, w, sx, sy, meta: dict) -> dict:
    """Rules 6-9 after ``m0417_arcref.run``: local noise, footprint, gate C2, flags."""
    rng = np.random.default_rng(m14.SEED)
    ctrl = Table.read(out / f"control_{FIELD}.ecsv")
    pred = Table.read(out / f"predicted_{FIELD}.ecsv")
    pix = meta["pixel_scale_arcsec"]
    for t in (ctrl, pred) if len(pred) else (ctrl,):
        m14.apply_local(t, data, w, sx, sy, rng, pix)
        apply_footprint(t, data, w)
    cl = [str(c) for c in ctrl["class_local"]]
    pairs = [
        (int(a), int(b)) for a, b in zip(ctrl["ref_counterpart"], ctrl["counterpart"], strict=True)
    ]
    c2 = ar.gate_c2_shared(cl, list(ctrl["arc_ref"]), [str(x) for x in ctrl["system"]], pairs)
    pre = [str(c) for c in ctrl["class_local_pre_fp"]]
    flags = [
        f"{r['system']}@{r['ra']:.6f},{r['dec']:.6f}"
        for r in pred
        if c2["passes"] and str(r["class_local"]) == "absent"
    ] if len(pred) else []  # fmt: skip
    ctrl.write(out / f"control_{FIELD}.ecsv", overwrite=True)
    pred.write(out / f"predicted_{FIELD}.ecsv", overwrite=True)
    return {
        "gate_C2_shared": c2,
        "flags": flags,
        "control_classes": m14._count(ctrl["class_local"]),
        "target_classes": m14._count(pred["class_local"]) if len(pred) else {},
        "gate_C2_plain": ar.gate_c2(cl, list(ctrl["arc_ref"])),
        "without_footprint_rule": {
            "gate_C2_shared": ar.gate_c2_shared(
                pre, list(ctrl["arc_ref"]), [str(x) for x in ctrl["system"]], pairs
            ),
            "control_classes": m14._count(pre),
        },
        "without_local_noise_rule": {
            k: s[k] for k in ("gate_C2", "control_classes", "target_classes", "flags")
        },
        "posthoc_absent_controls": s["posthoc_absent_controls"],
        "inputs": s["inputs"],
        "rules": s["rules"],
    }


def match_counts(ra, dec, cra, cdec, radii=POSTHOC_RADII) -> dict:
    """Post hoc (not a rule): images with a catalogue row within each radius (arcsec)."""
    from astropy.coordinates import SkyCoord

    g = SkyCoord(np.asarray(ra, float), np.asarray(dec, float), unit="deg")
    c = SkyCoord(np.asarray(cra, float), np.asarray(cdec, float), unit="deg")
    sep = g.match_to_catalog_sky(c)[1].arcsec
    return {f"{r:g}": int(np.sum(sep <= r)) for r in radii}


def posthoc(cache: Path) -> dict:
    """Why frame check F failed (not a rule): nearest-row counts against the own and the
    pipeline catalogues, and images inside a deblended segment, for all and spec-z images."""
    import s3fs

    images = mg.load_images()
    spec = np.asarray(images["z_spec"], bool)
    table = Table.read(FRAME / f"predicted_{FIELD}.ecsv")
    ctrl = Table.read(FRAME / f"control_{FIELD}.ecsv")
    pos = vstack([table["ra", "dec"], ctrl["ra", "dec"]])
    fs = s3fs.S3FileSystem(anon=True, default_block_size=64 * 1024)
    sci, err, w, _ = md.read_section(fs, S3_6882.format("i2d.fits")[5:], pos["ra"], pos["dec"])
    cat, _, deb, _, _ = md.build(sci, err, w)
    shapes = lc.load_shapes(cache / Path(S3_6882.format("cat.ecsv")).name)
    x, y = (np.rint(v).astype(int) for v in w.world_to_pixel_values(images["ra"], images["dec"]))
    in_seg = deb.data[y, x] > 0
    out = {"note": "post hoc, not a frozen rule (PLAN.md stage 3 rule 5 failed)"}
    for name, m in (("all_used", np.ones(len(images), bool)), ("spec_z", spec)):
        out[name] = {
            "n_images": int(m.sum()),
            "own_within_arcsec": match_counts(
                images["ra"][m], images["dec"][m], cat["ra"], cat["dec"]
            ),
            "pipeline_within_arcsec": match_counts(
                images["ra"][m], images["dec"][m], shapes["ra"], shapes["dec"]
            ),
            "inside_own_segment": int(in_seg[m].sum()),
        }
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--maps", type=Path, required=True)
    ap.add_argument("--cache", type=Path, default=Path("data/cache/m1931s3"))
    ap.add_argument("--posthoc", action="store_true", help="frame check F diagnostic only")
    a = ap.parse_args(argv)
    if a.posthoc:
        s = posthoc(a.cache)
        (OUT / "posthoc.json").write_text(json.dumps(s, indent=1) + "\n")
        print(json.dumps(s, indent=1))
        return 0
    OUT.mkdir(parents=True, exist_ok=True)
    s = cp._clean(run(a.maps, a.cache))
    (OUT / "summary.json").write_text(json.dumps(s, indent=1, allow_nan=False) + "\n")
    keys = ("rule2_gate_K_wls", "rule5_frame_F", "stopped")
    print(json.dumps({k: s[k] for k in keys if k in s}, indent=1))
    if "frozen" in s:
        print(json.dumps({k: s["frozen"][k] for k in ("gate_C2_shared", "flags")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
