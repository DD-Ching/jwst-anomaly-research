"""MACS1149: sibling reference fluxes and a second-model check for the CATS ``no_reference`` rows.

The CATS v4.1 ``images`` run (docs/fields/macs1149.md step 4) forced-tested 16 predicted images in
F277W; 11 got ``no_reference`` because no catalogued sibling passed the D-034 reference gates.
This script tests them the way ``a2744_xmodel_refs.py`` tested Abell 2744's:

- each catalogued sibling's peak flux (0.4" search, 0.2" aperture, 1" median-filtered stamp) divided
  by its CATS |mu| and times |mu_target| as an S/N against the target's aperture error, per sibling
  and from the summed siblings (``pred_snr_pair``), in F150W, F277W and F444W, and whether each
  sibling passes the D-034 gates (without the |mu| cap);
- the best aperture S/N within the system's forced search radius and
  ``lens_consistency.forced_class``
  (``forced_class`` only when every sibling passes the gates; ``ungated_class`` as a diagnostic);
- Sharon v4cor (HFF Lenstool deflection maps, the D-037 second model) solved for the same system:
  the source is the magnification-weighted least-squares source of the catalogued images
  (``a2744_uncover.family_source``) at the CATS system redshift, and the nearest Sharon image to the
  Sharon images within ``au.MATCH_ARCSEC`` of a catalogued image are that image; the rest are
  extra. D-041 rule: a CATS extra image is "both models" when a Sharon extra lies within
  ``XMODEL_LOOSE_ARCSEC``; else "catalogued image misplaced by CATS" when Sharon's image of a
  catalogued image that CATS does not reproduce lies there; else model-dependent.

Inputs: the ``images_predicted.ecsv`` of a ``lens_consistency.py --model macs1149-cats images
--forced-image <F277W i2d>`` run (``model_prediction``), the CATS ``arcs.txt`` positions
(``observed``) and the two Sharon maps (local paths; sha256 checked against
docs/fields/macs1149.md).
Photometry is ``derived`` from S3 byte-range stamps of CANUCS program 1208. Thresholds are
ASSUMPTIONs.

The tracked input is ``results/lens_consistency/macs1149-cats/images_predicted.ecsv``
(2026-10-14 run, ``--step 0.25``, F200W catalogue, F277W forced image).

Usage: python scripts/m1149_refs.py --sharon-x <fits> --sharon-y <fits> [--pred <ecsv>]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import a2744_uncover as au  # noqa: E402
import a2744_xmodel_refs as xr  # noqa: E402
import lens_consistency as lc  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402
from jwst_anomaly.derived_store import sha256_file  # noqa: E402

S3_I2D = "s3://stpubdata/jwst/public/jw01208/L3/t/o008/jw01208-o008_t004_nircam_clear-{}_i2d.fits"
BANDS = ("f150w", "f277w", "f444w")
SHARON_SHA256 = {
    "x": "44f1d21a78cead051c08295d13d7c9189638971cf1e9445bb2caa641a129324f",
    "y": "54508deb502d9bcafca9733c3b28f3b4d794b64f11167589b56c1fc37d66a15c",
}
SHARON_Z_LENS = 0.542  # map headers carry none (docs/fields/macs1149.md)
XMODEL_LOOSE_ARCSEC = au.XMODEL_LOOSE_ARCSEC  # D-041 addendum loose radius (ASSUMPTION)
DEFAULT_OUT = Path("results/lens_consistency/macs1149-cats")
_r = au._r


def sharon_model(x_path: Path, y_path: Path) -> lensmodel.MapLensModel:
    for k, p in (("x", x_path), ("y", y_path)):
        if sha256_file(p) != SHARON_SHA256[k]:
            raise ValueError(f"{p}: sha256 differs from the pinned Sharon v4cor {k} map")
    return lensmodel.MapLensModel.from_fits(
        x_path, y_path, SHARON_Z_LENS, FlatLambdaCDM(H0=70, Om0=0.3), source="macs1149-sharon-v4cor"
    )


def second_model_images(model, grid, cat: Table, z: float) -> Table:
    """Images ``model`` predicts for the source of the catalogued images ``cat`` at ``z``."""
    ox, oy = model.to_frame(np.asarray(cat["ra"], float), np.asarray(cat["dec"], float))
    m = model.lens_map(ox, oy, np.full(len(cat), z))
    jac = np.stack([np.stack([m["a11"], m["a12"]], -1), np.stack([m["a21"], m["a22"]], -1)], -2)
    beta = np.stack([np.asarray(m["beta_x"], float), np.asarray(m["beta_y"], float)], -1)
    ok = np.all(np.isfinite(beta), axis=1) & np.all(np.isfinite(jac), axis=(1, 2))
    ok &= np.abs(np.linalg.det(np.where(np.isfinite(jac), jac, 0.0))) > 1e-6  # off critical curves
    if ok.sum() < 1:
        return Table(names=("ra", "dec", "magnification"))
    bx, by = au.wls_source(beta[ok], jac[ok])
    return lensmodel.find_images(model, grid, float(bx), float(by), z)


def split_second(sh: Table, cat: Table) -> np.ndarray:
    """Catalogued image id each second-model image reproduces ("" for an extra image): the
    nearest catalogued image within ``au.MATCH_ARCSEC``, as the ``images`` classes."""
    ids = np.full(len(sh), "", dtype=object)
    for k in range(len(sh)):
        d = au._sep_arcsec(cat["ra"], cat["dec"], float(sh["ra"][k]), float(sh["dec"][k]))
        if len(d) and float(np.min(d)) <= au.MATCH_ARCSEC:
            ids[k] = str(cat["image_id"][int(np.argmin(d))])
    return ids


def d041_verdict(sep_extra: float, sep_obs: float, obs_predicted_by_cats: bool) -> str:
    """D-041 cross-model reading of one CATS extra image. ``sep_extra``: arcsec to the second
    model's nearest extra image; ``sep_obs``: to its nearest image of a catalogued one;
    ``obs_predicted_by_cats``: CATS itself reproduces that catalogued image."""
    if np.isfinite(sep_extra) and sep_extra <= XMODEL_LOOSE_ARCSEC:
        return "both models"
    if np.isfinite(sep_obs) and sep_obs <= XMODEL_LOOSE_ARCSEC and not obs_predicted_by_cats:
        return "catalogued image misplaced by CATS"
    return "model-dependent (CATS only)"


def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pred", type=Path, default=DEFAULT_OUT / "images_predicted.ecsv")
    ap.add_argument("--sharon-x", type=Path, required=True)
    ap.add_argument("--sharon-y", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args(argv)
    pred = Table.read(args.pred)
    if "forced_class" not in pred.colnames:
        raise SystemExit(f"{args.pred}: not an images run with --forced-image")
    fc = np.array([au._col(pred, "forced_class", k, "") for k in range(len(pred))]).astype(str)
    tgt = pred[fc == "no_reference"]
    _, files, par = lc.load_model("macs1149-cats")
    arcs, _ = lc.image_list("macs1149-cats", files, par)
    sysid = np.asarray(arcs["system"]).astype(str)
    sharon = sharon_model(args.sharon_x, args.sharon_y)
    systems = list(dict.fromkeys(str(s) for s in tgt["system"]))
    cats_all = arcs[np.isin(sysid, systems)]
    half = lc.grid_half_width(sharon, cats_all)
    grid = lensmodel.DeflectionGrid.compute(sharon, half, 0.25)
    halfst = lc.MAX_SEARCH_ARCSEC + lc.FORCED_ANNULUS[1] + 0.2
    opened = {b: lc.image_stamper(S3_I2D.format(b), halfst) for b in BANDS}
    stamps = {b: lc.halo_stamper(o[0], au.HALO_FILTER_ARCSEC) for b, o in opened.items()}
    rows = []
    try:
        for s in systems:
            cat = arcs[sysid == s]
            z = float(pred["z_sys"][np.asarray(pred["system"]).astype(str) == s][0])
            sibs = {b: xr.sibling_fluxes(cat, st) for b, st in stamps.items()}
            mus = xr.sibling_mu(pred, s)
            has_unpred = any(
                lensmodel.image_family(u) == s for u in (pred.meta.get("unpredicted") or [])
            )
            search = lc.system_search_radius(pred, s, 1.0, has_unpred)
            sh = second_model_images(sharon, grid, cat, z)
            sh_ids = split_second(sh, cat)
            for t in tgt[np.asarray(tgt["system"]).astype(str) == s]:
                ra, dec, mu_t = float(t["ra"]), float(t["dec"]), float(t["magnification"])
                row = {
                    "system": s,
                    "z": z,
                    "radec": [round(ra, 6), round(dec, 6)],
                    "image_class": str(t["image_class"]),
                    "mu_cats": _r(mu_t),
                    "search_arcsec": _r(search),
                    "sibling_mu_cats": {k: _r(v, 1) for k, v in mus.items()},
                    "catalogued": [str(i) for i in cat["image_id"]],
                }
                d = au._sep_arcsec(sh["ra"], sh["dec"], ra, dec) if len(sh) else np.array([])
                ext, obs = np.flatnonzero(sh_ids == ""), np.flatnonzero(sh_ids != "")
                ke = int(ext[np.argmin(d[ext])]) if len(ext) else -1
                ko = int(obs[np.argmin(d[obs])]) if len(obs) else -1
                sep_e = float(d[ke]) if ke >= 0 else np.nan
                sep_o = float(d[ko]) if ko >= 0 else np.nan
                obs_id = str(sh_ids[ko]) if ko >= 0 else ""
                row["sharon"] = {
                    "n_images": len(sh),
                    "n_extra": len(ext),
                    "n_catalogued": len(cat),
                    "nearest_extra_sep_arcsec": _r(sep_e),
                    "nearest_extra_mu": _r(float(sh["magnification"][ke])) if ke >= 0 else None,
                    "nearest_catalogued_sep_arcsec": _r(sep_o),
                    "nearest_catalogued_id": obs_id,
                    "nearest_catalogued_predicted_by_cats": obs_id in mus,
                    "d041": d041_verdict(sep_e, sep_o, obs_id in mus),
                }
                for band, stamp in stamps.items():
                    sib = sibs[band]
                    img, err, xx, yy = stamp(ra, dec)
                    _, e_t = lc.aperture_snr(img, err, xx, yy, 0.0, 0.0)
                    p = xr.predictions(
                        {k: (v["flux"], mus.get(k, np.nan)) for k, v in sib.items()}, mu_t, e_t
                    )
                    snr, flux, dx, dy = lc.best_within(img, err, xx, yy, search)
                    ratio = flux / p["pair_flux"] if p["pair_flux"] > 0 else np.nan
                    fclass, ungated = xr.classify(p, sib, e_t, snr, ratio)
                    row[band] = {
                        "siblings": {
                            k: {
                                "flux": _r(v["flux"], 4),
                                "snr": _r(v["snr"], 1),
                                "conc": _r(v["conc"]),
                                "gates_ok": v["gates_ok"],
                            }
                            for k, v in sib.items()
                        },
                        "pred_snr_from": {k: _r(v, 1) for k, v in p["from"].items()},
                        "pred_snr_pair": _r(p["pair"], 1),
                        "best_snr": _r(snr, 1),
                        "best_offset": [_r(dx, 1), _r(dy, 1)],
                        "flux_over_pair_pred": _r(ratio),
                        "forced_class": fclass,
                        "ungated_class": ungated,
                    }
                rows.append(row)
    finally:
        for _, close in opened.values():
            close()
    summary = {
        "provenance": "positions and mu: model_prediction (CATS v4.1, Sharon v4cor); catalogued "
        "positions: observed; fluxes: derived",
        "images": {b: S3_I2D.format(b) for b in BANDS},
        "assumptions": {
            "halo_filter_arcsec": au.HALO_FILTER_ARCSEC,
            "xmodel_loose_arcsec": XMODEL_LOOSE_ARCSEC,
            "match_arcsec": au.MATCH_ARCSEC,
            "sharon_z_lens": SHARON_Z_LENS,
            "sharon_source": "magnification-weighted least-squares source of the catalogued images",
            "class": "forced_class only when every sibling passes the D-034 gates, else "
            "no_reference (ungated_class is a diagnostic)",
            "mu_cap": "none",
        },
        "rows": rows,
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "no_reference_refs.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


if __name__ == "__main__":
    s = main()
    for r in s["rows"]:
        print(
            r["system"],
            r["radec"],
            r["mu_cats"],
            r["sharon"]["nearest_extra_sep_arcsec"],
            r["sharon"]["nearest_catalogued_sep_arcsec"],
            r["sharon"]["d041"],
            *(
                f"{b}:{r[b]['forced_class']}/{r[b]['ungated_class']} "
                f"pred={r[b]['pred_snr_pair']} best={r[b]['best_snr']}"
                for b in BANDS
            ),
        )
