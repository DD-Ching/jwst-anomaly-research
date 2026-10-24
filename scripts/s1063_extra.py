"""Abell S1063: extra-image test with the WSLAP+ JWST 2025 maps (WLS source), DJA v7 reference.

Frozen plan: ``results/lens_consistency/s1063_extra/PLAN.md``. Rule 1 (frame check F) compares the
36 Chan+2020 z_spec images with the DJA v7 ``abells1063`` rows; rule 2 re-fetches the two alpha
members (``s1063_wslap.fetch_member``) and repeats gate K-WLS, which must reproduce
``../s1063_wslap/summary.json``; rule 3 builds the frame tables with ``m0417_extra.run_field`` and
the WLS source on GO 3293 ``o001_t001``; rules 4-8 are El Gordo stage 2b's DJA reference
(``elgordo_djaref``), local noise and gate C2 shared (``a2744_uncoverref._reference_run``).

    python scripts/s1063_extra.py --maps DIR [--cache data/cache/s1063x]

The maps are deleted after the run unless ``--keep``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import a2744_uncover as au
import a2744_uncoverref as uref
import bullet_extra as be
import elgordo_djaref as eg
import lens_consistency as lc
import m0417_canucsref as cr
import m0417_extra as mx
import numpy as np
import requests
import s1063_wslap as sw
from astropy.coordinates import SkyCoord
from m0417_canucsref_posthoc import _clean

FIELD = "abells1063"
MODEL = "abells1063-wslap"
OUT = Path("results/lens_consistency/s1063_extra")
FRAME = OUT / "frame"
GATEK = Path("results/lens_consistency/s1063_wslap/summary.json")
IMAGES_SPEC = "chan20_zspec_images.ecsv"  # the 36 gate-K images, written by the run
OBS, TILE = "o001", "t001"  # GO 3293, TARGPROP ACO-S-1063 (PLAN.md footprint)
S3_3293 = "s3://stpubdata/jwst/public/jw03293/L3/t/{0}/jw03293-{0}_{1}_nircam_clear-f200w_{2}"
DJA_PHOT = "https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/abells1063-grizli-v7.0-fix_phot.fits"
GRID_STEP_K = sw.GRID_STEP  # 0.25 arcsec, as gate K of ../s1063_wslap
K_REPEAT_TOL_ARCSEC = 1e-6  # rule 2: same map bytes -> same rms
F_MIN_N = 18  # rule 1 (ASSUMPTION): half of the 36 images; match radius and medians as Bullet


def spec_images_table(images):
    """The gate-K images in ``m1149_sharon.gold_images`` form (``system`` / ``a`` re-derived)."""
    t = images.copy()
    t.remove_columns(["system", "a"])
    return t


def k_repeat(k_wls: dict, stored: dict, tol: float = K_REPEAT_TOL_ARCSEC) -> dict:
    """Rule 2: gate K-WLS passes again and reproduces the stored K-WLS (n, rms, fraction)."""
    same = (
        k_wls["n_images"] == stored["n_images"]
        and abs(k_wls["rms_arcsec"] - stored["rms_arcsec"]) <= tol
        and abs(k_wls["matched_fraction"] - stored["matched_fraction"]) <= 1e-12
    )
    return {
        "K_wls": k_wls,
        "stored": stored,
        "reproduced": bool(same),
        "passes": bool(same and k_wls["passes"]),
    }


def register() -> None:
    """Point ``m0417_extra``'s field table at this field, tile and directory."""
    mx.FIELDS[FIELD] = (IMAGES_SPEC, MODEL, OBS, TILE)
    mx.IMAGES = OUT
    mx.OUT = FRAME
    mx.S3 = S3_3293


def fetch_maps(maps: Path) -> tuple:
    with requests.Session() as session:
        _, directory = sw.zip_directory(session)
    maps.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(len(sw.MAPS)) as pool:
        futs = {
            m: pool.submit(sw.fetch_member, requests.Session(), directory[sw.ZIP_DIR + m], maps / m)
            for m in sw.MAPS
        }
        shas = {m: f.result() for m, f in futs.items()}
    (ax, hx), (ay, _) = (sw._read(maps / m) for m in sw.MAPS)
    model_sha = hashlib.sha256((shas[sw.MAPS[0]] + shas[sw.MAPS[1]]).encode()).hexdigest()
    return sw.build_model(ax, ay, hx, sw.READINGS["A_arcsec"], model_sha, MODEL), shas


def run(maps: Path, cache: Path) -> dict:
    import s3fs

    stored = json.loads(GATEK.read_text())
    images = sw.spec_images()
    OUT.mkdir(parents=True, exist_ok=True)
    spec_images_table(images).write(OUT / IMAGES_SPEC, format="ascii.ecsv", overwrite=True)
    t, sha, nbytes = eg.fetch_table(DJA_PHOT)
    summary: dict = {
        "plan": str(OUT / "PLAN.md"),
        "dja": {"url": DJA_PHOT, "sha256": sha, "bytes": nbytes, "rows": len(t)},
    }
    f = be.frame_check_f(images["ra"], images["dec"], t["ra"], t["dec"], min_n=F_MIN_N)
    summary["rule1_frame_F"] = f
    if not f["passes"]:
        summary["stopped"] = "frame check F failed (PLAN.md rule 1)"
        return summary
    model, shas = fetch_maps(maps)
    summary["maps_sha256"] = shas
    k = sw.k_block(model, images, GRID_STEP_K)
    rep = k_repeat(k["K_wls"], stored["reported"]["K_wls_not_a_gate"])
    summary["rule2_gate_K_wls"] = rep
    if not rep["passes"]:
        summary["stopped"] = "gate K-WLS not reproduced or failed (PLAN.md rule 2)"
        return summary
    lc._write(k["ip_wls"], OUT / f"imageplane_{FIELD}_{MODEL}.ecsv")  # z_used for rule 3
    unit, grid = eg.flux_unit_jy(t, with_tot_corr=True), eg.grid_rotation_deg(t)
    summary["rule4a_unit"], summary["rule4b_grid"] = unit, grid
    if not (unit["passes"] and grid["passes"]):
        summary["stopped"] = "rule 4a / 4b failed (PLAN.md rule 4)"
        return summary
    register()
    fs = s3fs.S3FileSystem(anon=True, default_block_size=64 * 1024)
    cache.mkdir(parents=True, exist_ok=True)
    shapes = cache / Path(S3_3293.format(OBS, TILE, "cat.ecsv")).name
    if not shapes.exists():
        fs.get(S3_3293.format(OBS, TILE, "cat.ecsv")[5:], str(shapes))
    FRAME.mkdir(parents=True, exist_ok=True)
    try:  # rule 3 (the model as gate K-WLS used it: no frame offset)
        frame = mx.run_field(FIELD, cache, bt_hook=au.wls_backtrace, model=model)
    finally:
        shapes.unlink(missing_ok=True)  # cloud disk rule: the pipeline catalogue is not kept
    summary["rule3_frame_tables"] = frame
    (FRAME / "summary.json").write_text(
        json.dumps(_clean({FIELD: frame}), indent=1, allow_nan=False) + "\n"
    )
    w, d, pix = uref._read_sci(fs, S3_3293.format(OBS, TILE, "i2d.fits")[5:])
    rc = eg.to_cat_dja(t, unit["unit_jy"], grid["pixel_scale_arcsec"])
    inputs = {
        "reference": Path(DJA_PHOT).name,
        "reference_sha256": sha,
        "n_rows_kept": int(len(rc["flux"])),
        "flux": f"{eg.BAND}_flux_aper_{eg.APER_IX} x flux_auto / flux_aper_{eg.APER_IX}"
        f" x {eg.BAND}_tot_corr",
        "noise_image": S3_3293.format(OBS, TILE, "i2d.fits"),
        "noise_pixel_scale_arcsec": pix,
    }
    depth = cr.depth_mag(rc["flux"], rc["snr"])
    src_x, src_y = (
        np.asarray(v, float) for v in w.world_to_pixel(SkyCoord(t["ra"], t["dec"], unit="deg"))
    )
    for contain, out in ((False, OUT / "without_contain"), (True, OUT)):  # rule 8 / rules 5-7
        out.mkdir(parents=True, exist_ok=True)
        summary["frozen" if contain else "without_contain"] = uref._reference_run(
            rc, depth, inputs, out, contain, d, w, src_x, src_y, pix, field=FIELD, frame=FRAME
        )
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--maps", type=Path, required=True)
    ap.add_argument("--cache", type=Path, default=Path("data/cache/s1063x"))
    ap.add_argument("--keep", action="store_true", help="keep the fetched maps")
    a = ap.parse_args(argv)
    try:
        s = _clean(run(a.maps, a.cache))
    finally:
        if not a.keep:
            for m in sw.MAPS:
                (a.maps / m).unlink(missing_ok=True)
    (OUT / "summary.json").write_text(json.dumps(s, indent=1, allow_nan=False) + "\n")
    keys = ("rule1_frame_F", "rule4a_unit", "rule4b_grid", "stopped")
    brief = {k: s[k] for k in keys if k in s}
    brief["rule2"] = {k: v for k, v in s.get("rule2_gate_K_wls", {}).items() if k != "stored"}
    brief["C2"] = s.get("frozen", {}).get("gate_C2_shared")
    brief["flags"] = s.get("frozen", {}).get("flags")
    print(json.dumps(brief, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
