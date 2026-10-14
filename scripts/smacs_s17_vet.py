"""SMACS 0723 system 17: is its halo-subtracted ``absent`` counter-image a lens-model residual?

``lens_consistency.py images --halo-filter-arcsec 1.0`` (D-029 revisit) turns the Mahler+22 ICLv2
prediction of a mu ~ -10 image of system 17 near the BCG from a marginal recovery (7.8 sigma on the
BCG's light) into ``absent``. System 17 has no spectroscopic redshift (arcs.dat z = 0: fitted), so
this script asks whether the predicted image survives the ordinary explanations:

- ``zscan``: the ICLv2 images of system 17's mean back-traced source over the model redshifts of
  the published models (Pascale+22 1.72, Caminha+22 1.84-1.92, Mahler+22 2.048 best / 2.12
  median, NIRCam photo-z 2.34), with halo-filtered forced photometry (best S/N within 1" and 2")
  of every image that is not a catalogued one, and the flux each catalogued sibling predicts
  (its halo-filtered peak flux / |mu| x |mu|);
- ``relics``: the same source traced through two independent HST-era RELICS deflection maps
  (Glafic v2, Zitrin LTM-Gauss v2): image count, magnification and distance to the ICLv2 images.

Thresholds are ASSUMPTIONs; image positions and magnifications are model_prediction, fluxes derived.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lens_consistency as lc  # noqa: E402

from jwst_anomaly import lensmodel, paths  # noqa: E402
from jwst_anomaly.photometry import fetch_catalog  # noqa: E402

MODEL = "smacs0723-iclv2"
F200W = "s3://stpubdata/jwst/public/jw02736/L3/t/o001/jw02736-o001_t001_nircam_clear-f200w_i2d.fits"
# model redshifts of system 17 in the published models and the NIRCam photo-z (SOURCES.md)
Z_SCAN = (1.5, 1.72, 1.84, 1.92, 2.048, 2.12, 2.34)
CAT_MATCH_ARCSEC = 1.5  # a predicted image this close to a catalogued one is that image
FILTER_ARCSEC = 1.0  # halo filter (as lens_consistency --halo-filter-arcsec; ASSUMPTION)
RELICS_URL = "https://archive.stsci.edu/hlsps/relics/smacs0723m73/models/"
RELICS = {  # (x, y) arcsec-deflection maps at D_LS/D_S = 1, sha256 as downloaded 2026-10-14
    "glafic_v2": (
        "glafic/v2/",
        "e3d8fc285a33e38ecd1bce9fc50418ce1d3ff9ceb251ba5e41b5ddd38e1ad74c",
        "51ff730475dfd2a541ca0a2f7037540fe2359bfd4701fa7939ae97df8f4d4a75",
    ),
    "zitrin-ltm-gauss_v2": (
        "zitrin-ltm-gauss/v2/",
        "cc85c434e939bcdae1b6fd1dde0b7c99971ae2048141d4cd10b568001bdaadf8",
        "6cdd40c8f070480ffa537935b9aa5109eb8132a72c995a5f000753ea1fab2c06",
    ),
}
RELICS_Z_LENS, RELICS_COSMO = 0.39, (70.0, 0.3)


def _r(v, nd=2):
    return None if v is None or not np.isfinite(v) else round(float(v), nd)


def trace(model, grid, s17, z: float) -> tuple[float, object]:
    """(source-plane spread of the catalogued images, all images of their mean source)."""
    bt = lensmodel.backtrace_images(model, s17, {"17": z})
    bx, by = np.nanmean(bt["beta_x"]), np.nanmean(bt["beta_y"])
    spread = float(np.nanmax(np.hypot(bt["beta_x"] - bx, bt["beta_y"] - by)))
    return spread, lensmodel.find_images(model, grid, bx, by, z), bt


def zscan(s17, stamp) -> list[dict]:
    model, _, _ = lc.load_model(MODEL)
    lc.apply_frame_offset(MODEL, model)
    grid = lensmodel.DeflectionGrid.compute(model, 45.0, 0.25)
    cat = np.array(model.to_frame(np.asarray(s17["ra"], float), np.asarray(s17["dec"], float))).T
    ref = {}
    for r in s17:  # halo-filtered peak flux within 0.4" of each catalogued image
        img, err, xx, yy = stamp(float(r["ra"]), float(r["dec"]))
        ref[str(r["image_id"])] = lc.peak_flux(img, err, xx, yy)[0]
    out = []
    for z in Z_SCAN:
        spread, im, bt = trace(model, grid, s17, z)
        lum = {
            str(i): ref[str(i)] / abs(float(m))
            for i, m in zip(bt["image_id"], bt["magnification"], strict=True)
        }
        rows = []
        for r in im:
            d = float(np.hypot(cat[:, 0] - r["x"], cat[:, 1] - r["y"]).min())
            row = {"xy": [_r(r["x"]), _r(r["y"])], "mu": _r(r["magnification"]), "d_cat": _r(d)}
            if d >= CAT_MATCH_ARCSEC:
                img, err, xx, yy = stamp(float(r["ra"]), float(r["dec"]))
                mu = abs(float(r["magnification"]))
                for rad in (1.0, 2.0):
                    snr, flux, dx, dy = lc.best_within(img, err, xx, yy, rad)
                    row[f"best_{rad:.0f}arcsec"] = {
                        "snr": _r(snr, 1),
                        "flux": _r(flux),
                        "offset": [_r(dx, 1), _r(dy, 1)],
                    }
                row["pred_flux_from"] = {k: _r(v * mu) for k, v in lum.items()}
            rows.append(row)
        out.append({"z": z, "source_spread_arcsec": _r(spread), "images": rows})
    return out


def relics(s17, ref_xy: list) -> list[dict]:
    centre = (float(np.mean(s17["ra"])), float(np.mean(s17["dec"])))
    cat = SkyCoord(s17["ra"], s17["dec"], unit="deg")
    ref = SkyCoord([p[0] for p in ref_xy], [p[1] for p in ref_xy], unit="deg")
    out = []
    for name, (sub, sha_x, sha_y) in RELICS.items():
        base = f"{RELICS_URL}{sub}hlsp_relics_model_model_smacs0723m73_{name}"
        px = fetch_catalog(f"{base}_x-arcsec-deflect.fits", sha_x)
        py = fetch_catalog(f"{base}_y-arcsec-deflect.fits", sha_y)
        h0, om = RELICS_COSMO
        model = lensmodel.MapLensModel.from_fits(
            px, py, RELICS_Z_LENS, FlatLambdaCDM(H0=h0, Om0=om), source=name, centre=centre
        )
        grid = lensmodel.DeflectionGrid.compute(model, 40.0, 0.25)
        for z in Z_SCAN:
            spread, im, _ = trace(model, grid, s17, z)
            ic = SkyCoord(im["ra"], im["dec"], unit="deg")
            rows = [
                {
                    "mu": _r(r["magnification"]),
                    "d_cat": _r(c.separation(cat).arcsec.min()),
                    "d_iclv2_central": [_r(v) for v in c.separation(ref).arcsec],
                }
                for r, c in zip(im, ic, strict=True)
            ]
            out.append({"model": name, "z": z, "source_spread_arcsec": _r(spread), "images": rows})
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=paths.repo_root() / "results" / "smacs_s17")
    ap.add_argument("--no-relics", action="store_true", help="skip the RELICS maps (118 MB)")
    args = ap.parse_args(argv)
    m0, files, par = lc.load_model(MODEL)
    images, _ = lc.image_list(MODEL, files, par)
    lc.apply_frame_offset(MODEL, m0, images)
    s17 = images[np.asarray(images["system"]).astype(str) == "17"]
    raw, close = lc.image_stamper(F200W, lc.MAX_SEARCH_ARCSEC + lc.FORCED_ANNULUS[1] + 0.2)
    try:
        scan = zscan(s17, lc.halo_stamper(raw, FILTER_ARCSEC))
    finally:
        close()
    best = next(s for s in scan if s["z"] == 2.048)
    central = [
        tuple(float(v) for v in m0.to_sky(*row["xy"]))
        for row in best["images"]
        if row["d_cat"] >= CAT_MATCH_ARCSEC
    ]
    summary = {
        "model": MODEL,
        "image": F200W,
        "provenance": "positions and mu: model_prediction; fluxes: derived",
        "assumptions": {
            "halo_filter_arcsec": FILTER_ARCSEC,
            "cat_match_arcsec": CAT_MATCH_ARCSEC,
            "z_scan": list(Z_SCAN),
            "relics_z_lens": RELICS_Z_LENS,
            "relics_cosmology": list(RELICS_COSMO),
        },
        "iclv2_zscan": scan,
        "relics": [] if args.no_relics else relics(s17, central),
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "vet.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({"written": str(args.out / "vet.json")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
