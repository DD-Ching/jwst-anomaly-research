"""SMACS 0723 systems 11, 16 and 26: can any reference flux make their counter-images testable?

D-029 left these systems ``no_reference``: every catalogued image has |mu| > 50 or S/N < 5 in
F200W. This script drops the |mu| and compactness gates and asks the direct question for each
uncatalogued image the model predicts (farther than ``CAT_MATCH_ARCSEC`` from every catalogued
image):

- the flux each catalogued sibling predicts (its peak flux within 0.4" / |mu_sib| x |mu_img|),
  as an S/N against the image's own 0.2" aperture error, in F200W and F444W;
- ``pred_snr_gated``: the same from the siblings that pass the D-034 reference gates of
  ``lens_consistency.forced_check`` other than its |mu| cap (``gates_ok``: S/N > 5 at the peak and
  at the catalogued position, peak <= 2x the position flux, f(0.2")/f(0.4") >= 0.6);
- ``need_mu_factor``: by how much |mu_sib| / |mu_img| would have to be overestimated for the
  brightest sibling to predict 10 sigma (the D-029 ``absent`` floor);
- the best forced S/N within ``SEARCH_ARCSEC`` (what is there).

Models: Mahler+22 ICLv2 at its fitted redshift (arcs.dat z = 0 for all three systems), and the
two RELICS HST-era deflection maps of ``smacs_s17_vet.py`` at the same redshift (ASSUMPTION; the
maps are deleted after use). Thresholds are ASSUMPTIONs; positions and mu are model_prediction,
fluxes derived.
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
import smacs_s17_vet as s17v  # noqa: E402

from jwst_anomaly import lensmodel, paths  # noqa: E402
from jwst_anomaly.photometry import fetch_catalog  # noqa: E402

MODEL = s17v.MODEL
SYSTEMS = ("11", "16", "26")
BANDS = {
    "f200w": s17v.F200W,
    "f444w": s17v.F200W.replace("f200w", "f444w"),
}
CAT_MATCH_ARCSEC = s17v.CAT_MATCH_ARCSEC
SEARCH_ARCSEC = 1.0  # D-029 forced search radius (ASSUMPTION)
ABSENT_SNR = 10.0  # D-029: an absent image must be predicted at >= 10 sigma (ASSUMPTION)
_r = s17v._r


def sibling_flux(cat, stamps) -> dict:
    """Peak flux, its S/N and f(0.2")/f(0.4") of every catalogued image, per band."""
    out = {}
    for band, stamp in stamps.items():
        out[band] = {}
        for r in cat:
            img, err, xx, yy = stamp(float(r["ra"]), float(r["dec"]))
            f, e = lc.peak_flux(img, err, xx, yy)
            f0, e0 = lc.aperture_snr(img, err, xx, yy, 0.0, 0.0)
            px, py = lc.peak_position(img, xx, yy)
            fb, _ = lc.aperture_snr(img, err, xx, yy, px, py, r_ap=0.4)
            snr, snr0 = (f / e if e > 0 else np.nan), (f0 / e0 if e0 > 0 else np.nan)
            conc = f / fb if np.isfinite(fb) and fb > 0 else np.nan
            out[band][str(r["image_id"])] = {
                "flux": f,
                "snr": snr,
                "snr_at_position": snr0,
                "peak_over_position": f / f0 if f0 > 0 else np.nan,
                "conc": conc,
                # the D-034 reference gates of lens_consistency.forced_check, without its |mu| cap
                "gates_ok": bool(
                    snr > 5
                    and snr0 > 5
                    and f <= lc.MAX_PEAK_OVER_POSITION * max(f0, 0)
                    and conc >= lc.MIN_REF_CONCENTRATION
                ),
            }
    return out


def measure(model, grid, cat, z: float, sib: dict, stamps) -> dict:
    """Every image of the system's mean back-traced source under one model, with the sibling
    predictions and forced photometry of the uncatalogued ones."""
    bt = lensmodel.backtrace_images(model, cat, {str(cat["system"][0]): z})
    bx, by = np.nanmean(bt["beta_x"]), np.nanmean(bt["beta_y"])
    spread = float(np.nanmax(np.hypot(bt["beta_x"] - bx, bt["beta_y"] - by)))
    mu_sib = {
        str(i): abs(float(m)) for i, m in zip(bt["image_id"], bt["magnification"], strict=True)
    }
    cat_c = SkyCoord(cat["ra"], cat["dec"], unit="deg")
    rows = []
    for r in lensmodel.find_images(model, grid, bx, by, z):
        c = SkyCoord(float(r["ra"]), float(r["dec"]), unit="deg")
        d = float(c.separation(cat_c).arcsec.min())
        mu = abs(float(r["magnification"]))
        row = {
            "radec": [round(float(r["ra"]), 7), round(float(r["dec"]), 7)],
            "xy": [_r(r["x"]), _r(r["y"])],
            "mu": _r(r["magnification"]),
            "d_cat": _r(d),
        }
        if d >= CAT_MATCH_ARCSEC:
            for band, stamp in stamps.items():
                img, err, xx, yy = stamp(float(r["ra"]), float(r["dec"]))
                _, e0 = lc.aperture_snr(img, err, xx, yy, 0.0, 0.0)
                snr, flux, dx, dy = lc.best_within(img, err, xx, yy, SEARCH_ARCSEC)
                pred = {
                    k: v["flux"] / mu_sib[k] * mu / e0
                    for k, v in sib[band].items()
                    if np.isfinite(v["flux"]) and v["flux"] > 0 and e0 > 0
                }
                top = max(pred.values(), default=np.nan)
                gated = [v for k, v in pred.items() if sib[band][k]["gates_ok"]]
                row[band] = {
                    "pred_snr_from": {k: _r(v, 1) for k, v in pred.items()},
                    "need_mu_factor": _r(ABSENT_SNR / top, 1) if top > 0 else None,
                    "pred_snr_gated": _r(max(gated), 1) if gated else None,
                    "best_snr": _r(snr, 1),
                    "best_flux_over_top_pred": _r(flux / (top * e0), 2) if top > 0 else None,
                    "best_offset": [_r(dx, 1), _r(dy, 1)],
                }
        rows.append(row)
    return {
        "z": z,
        "source_spread_arcsec": _r(spread),
        "mu_catalogued": {k: _r(v, 1) for k, v in mu_sib.items()},
        "images": rows,
    }


def relics_models(centre, keep_maps: bool):
    """Yield ``(name, model)`` for the RELICS maps in the JWST frame (as ``smacs_s17_vet``)."""
    dra, ddec = lc.MODELS[MODEL].get("frame_offset_arcsec", (0.0, 0.0))
    for name, (sub, sha_x, sha_y) in s17v.RELICS.items():
        base = f"{s17v.RELICS_URL}{sub}hlsp_relics_model_model_smacs0723m73_{name}"
        px = fetch_catalog(f"{base}_x-arcsec-deflect.fits", sha_x)
        py = fetch_catalog(f"{base}_y-arcsec-deflect.fits", sha_y)
        try:
            h0, om = s17v.RELICS_COSMO
            model = lensmodel.MapLensModel.from_fits(
                px, py, s17v.RELICS_Z_LENS, FlatLambdaCDM(H0=h0, Om0=om), source=name, centre=centre
            )
        finally:
            if not keep_maps:
                px.unlink(missing_ok=True)
                py.unlink(missing_ok=True)
        model.shift_frame(dra, ddec)
        yield name, model


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=paths.repo_root() / "results" / "smacs_ref")
    ap.add_argument("--no-relics", action="store_true", help="skip the RELICS maps (118 MB)")
    ap.add_argument("--keep-maps", action="store_true", help="keep the RELICS maps after use")
    args = ap.parse_args(argv)
    m0, files, par = lc.load_model(MODEL)
    images, zml = lc.image_list(MODEL, files, par)
    centre = (m0.ra0, m0.dec0)  # the ICLv2 frame before its offset (the RELICS maps get it below)
    lc.apply_frame_offset(MODEL, m0, images)
    sysid = np.asarray(images["system"]).astype(str)
    cats = {s: images[sysid == s] for s in SYSTEMS}
    half = lc.MAX_SEARCH_ARCSEC + lc.FORCED_ANNULUS[1] + 0.2
    opened = {b: lc.image_stamper(u, half) for b, u in BANDS.items()}
    stamps = {b: o[0] for b, o in opened.items()}
    try:
        sib = {s: sibling_flux(cats[s], stamps) for s in SYSTEMS}
        models = [("iclv2", m0, lensmodel.DeflectionGrid.compute(m0, 45.0, 0.25))]
        if not args.no_relics:
            for name, m in relics_models(centre, args.keep_maps):
                models.append((name, m, lensmodel.DeflectionGrid.compute(m, 45.0, 0.25)))
        runs = [
            {"model": name, "system": s, **measure(m, g, cats[s], float(zml[s]), sib[s], stamps)}
            for name, m, g in models
            for s in SYSTEMS
        ]
    finally:
        for _, close in opened.values():
            close()
    summary = {
        "model_family": MODEL,
        "images": BANDS,
        "provenance": "positions and mu: model_prediction; fluxes: derived; z: model-fitted",
        "assumptions": {
            "cat_match_arcsec": CAT_MATCH_ARCSEC,
            "search_arcsec": SEARCH_ARCSEC,
            "absent_snr": ABSENT_SNR,
            "z": {s: round(float(zml[s]), 3) for s in SYSTEMS},
            "relics_z": "the ICLv2 fitted z of each system",
            "relics_frame": "shifted by the ICLv2 frame_offset_arcsec (same HST frame assumed)",
        },
        "siblings": {
            s: {
                b: {
                    k: {q: x if isinstance(x, bool) else _r(x, 2) for q, x in v.items()}
                    for k, v in d.items()
                }
                for b, d in sb.items()
            }
            for s, sb in sib.items()
        },
        "runs": runs,
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "ref_fluxes.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({"written": str(args.out / "ref_fluxes.json")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
