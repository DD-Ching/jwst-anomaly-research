"""Abell 2744: the UNCOVER v2.0 lens model (Furtak+2023, Price+2025) as an independent second model.

Step 1 of the UNCOVER cross-check (TASKS.md "Supporting" item 1): validate the published v2.0
best-fit deflection maps (0.1"/pix, D_LS/D_S = 1) before any counter-image test uses them.

* **Known cases (own list):** UNCOVER's own multiple-image table (``A2744_multiple-images_v2.0``,
  ids ``S.KI``: source ``S.K``, image ``I``), every image at its catalogued redshift, traced
  through ``lensmodel.MapLensModel`` and re-imaged from one source per family. The release quotes
  a lens-plane rms of 0.60"; ours must come within ``lens_consistency.IMAGE_LIST_RMS_FACTOR`` of
  it (ASSUMPTION, as for the HFF maps, D-035) before the model is used. The source is the
  magnification-weighted least-squares source (``wls_source``): the mean back-traced source
  fails the gate (rms 4.8") because it falls outside the caustic for fold pairs (34.1, 80.1);
  the solver grid step (0.25" vs 0.1") changes nothing.
* **Cross-model (Bergamini+23 list):** the Bergamini+23 ``arcs.dat`` images with their
  ``best.par`` redshifts through the UNCOVER maps: per family the image-plane residuals and the
  number of predicted images against the number catalogued. Families with predicted images that
  are not catalogued (``n_extra``) are where the two models disagree on multiplicity; those are
  the inputs of the D-041 cross-model rule in the next step.

Everything computed is ``model_prediction``; the image positions are ``observed``.

Data: the maps are fetched once from the UNCOVER Google Drive release (sha256-pinned; 216 MB
compressed in three files, ~250 MB on disk while the run lasts; the maps are deleted at the end
unless ``--keep-maps``).

* **4.2c (D-041 second-model test):** the Bergamini+23 ``absent`` flag 4.2c
  (docs/fields/abell2744.md). ``check_42c`` solves family 4.2 in the UNCOVER maps and measures
  1"-high-pass F150W/F277W/F444W peak fluxes (r = 0.2") at the UNCOVER 4.2c and at 4.2a/b, plus
  injections at 4.2c of the most magnified catalogued image's knot, scaled by the mu ratio.

* **Cross-model rule on every forced row (D-041 addendum):** ``--crossmodel`` takes the
  Bergamini+23 ``images_predicted.ecsv`` of ``lens_consistency.py images --forced-image
  --halo-filter-arcsec 1.0``, force-photometers UNCOVER's own predicted images the same way
  (``forced_images``), and lists every extra image of either model with the other model's
  nearest prediction (``crossmodel``).

Usage: python scripts/a2744_uncover.py [--out results/lens_consistency/abell2744-uncover2]
       [--check-42c] [--crossmodel IMAGES_PREDICTED --catalog F200W_CAT --forced-image F277W_I2D]
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
import lens_consistency as lc  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402
from jwst_anomaly.photometry import fetch_catalog  # noqa: E402

_DRIVE = "https://drive.usercontent.google.com/download?export=download&confirm=t&id="
# UNCOVER DR4 lensing maps v2.0 (2024-08-07), Best-model_low-resolution_100mas; accessed 2026-10-14
UNCOVER_V2 = {
    "alpha_x": (
        _DRIVE + "1BHZ63x9nYAd5O0s1K_mO1UndyaekeJjz",
        "7ba107377e1585327d9ad95f8d847941dbcb7864611c442af6b712b455df4803",
    ),
    "alpha_y": (
        _DRIVE + "1Ht5KmqeXifqI8Q0tSRcVzPvEz6iAaH3w",
        "23469d560df560f6c7136cf273695aebe0fedd09853f256599b55c30335d156a",
    ),
    "images": (
        _DRIVE + "1uQQQBA1W_zd-y8SILZfBYKpVGAGoTxMG",
        "0d2649b3991129df0b775965caadf10109714d32cd804dcc76afa17c9c789119",
    ),
}
Z_LENS = 0.308
COSMOLOGY = (70.0, 0.3)  # ASSUMPTION: the Zitrin-code default; maps are at D_LS/D_S = 1 anyway
QUOTED_RMS = 0.60  # README v2.0, lens-plane rms (arcsec)
MATCH_ARCSEC = 1.5  # a predicted image within this of a catalogued one is that image (ASSUMPTION)


def uncover2_model(files: dict[str, Path]) -> lensmodel.MapLensModel:
    """The v2.0 maps (``files["alpha_x"]``, ``files["alpha_y"]``) as a ``MapLensModel``."""
    h0, om = COSMOLOGY
    return lensmodel.MapLensModel.from_fits(
        files["alpha_x"],
        files["alpha_y"],
        Z_LENS,
        FlatLambdaCDM(H0=h0, Om0=om),
        source="abell2744-uncover2",
    )


def wls_source(beta: np.ndarray, jac: np.ndarray) -> np.ndarray:
    """Magnification-weighted least-squares source of one family.

    ``beta`` (n, 2) are the back-traced sources, ``jac`` (n, 2, 2) the lens Jacobians A at the
    images. Minimises sum |A_i^-1 (beta - beta_i)|^2, the image-plane offsets to first order
    (the source-plane chi^2 of Keeton 2001 / Oguri 2010 without position errors)."""
    minv = np.linalg.inv(jac)
    w = np.einsum("nki,nkj->nij", minv, minv)
    return np.linalg.solve(w.sum(0), np.einsum("nij,nj->i", w, beta))


def family_source(model, obs: Table) -> tuple[float, float]:
    """``wls_source`` of a family's back-traced rows (``backtrace_images`` output)."""
    z = float(obs["z_used"][0])
    ox, oy = model.to_frame(np.asarray(obs["ra"], float), np.asarray(obs["dec"], float))
    m = model.lens_map(ox, oy, z)
    jac = np.stack([np.stack([m["a11"], m["a12"]], -1), np.stack([m["a21"], m["a22"]], -1)], -2)
    beta = np.stack([np.asarray(obs["beta_x"], float), np.asarray(obs["beta_y"], float)], -1)
    bx, by = wls_source(beta, jac)
    return float(bx), float(by)


def imageplane(model, grid, bt: Table) -> Table:
    """Per catalogued image: offset to the nearest predicted image of its family's
    ``family_source`` (``lensmodel.imageplane_residuals`` with the WLS source)."""
    systems = np.asarray(bt["system"]).astype(str)
    dtheta = np.full(len(bt), np.nan)
    mu = np.full(len(bt), np.nan)
    n_pred = np.zeros(len(bt), int)
    shared = np.zeros(len(bt), bool)
    for fam in dict.fromkeys(systems):
        sel = np.flatnonzero((systems == fam) & np.isfinite(bt["beta_x"]))
        if len(sel) < 2:
            continue
        obs = bt[sel]
        pred = lensmodel.find_images(
            model, grid, *family_source(model, obs), float(obs["z_used"][0])
        )
        n_pred[sel] = len(pred)
        if not len(pred):
            continue
        ox, oy = model.to_frame(np.asarray(obs["ra"], float), np.asarray(obs["dec"], float))
        d = np.hypot(ox[:, None] - pred["x"][None, :], oy[:, None] - pred["y"][None, :])
        k = np.argmin(d, axis=1)
        dtheta[sel] = d[np.arange(len(sel)), k]
        mu[sel] = np.asarray(pred["magnification"])[k]
        shared[sel] = np.bincount(k, minlength=len(pred))[k] > 1
    return Table(
        {
            "image_id": bt["image_id"],
            "system": systems,
            "z_used": bt["z_used"],
            "dtheta_arcsec": dtheta,
            "magnification": mu,
            "n_predicted": n_pred,
            "shared_match": shared,
        }
    )


def uncover_system(image_id: str) -> str:
    """Source of an UNCOVER image id: "S.KI" (two decimals) -> "S.K"; "S.I" (one decimal,
    single-knot systems such as 83.1-83.3) -> "S". Ids are stored as floats, so "S.K0" cannot
    occur (image numbers start at 1)."""
    whole, frac = image_id.split(".")
    return f"{whole}.{frac[0]}" if len(frac) == 2 else whole


def uncover_images(path: Path) -> Table:
    """UNCOVER's multiple-image table as a ``load_lenstool_images``-like table."""
    t = Table.read(path)
    ids = [repr(float(i)) for i in t["ID"]]
    out = Table(
        {
            "image_id": ids,
            "system": [uncover_system(i) for i in ids],
            "ra": np.asarray(t["RA"], float),
            "dec": np.asarray(t["Dec"], float),
            "a": np.full(len(t), np.nan),
            "z": np.asarray(t["z"], float),
        }
    )
    out.meta.update(provenance="observed", source=str(path))
    return out


def multiplicity(model, grid, bt: Table) -> Table:
    """Per family: catalogued and predicted image counts, and predicted images more than
    ``MATCH_ARCSEC`` from every catalogued image (``n_extra``; their |mu| listed)."""
    rows = []
    for fam in np.unique(np.asarray(bt["system"]).astype(str)):
        sel = bt[(np.asarray(bt["system"]).astype(str) == fam) & np.isfinite(bt["beta_x"])]
        if len(sel) < 2:
            continue
        z_s = float(sel["z_used"][0])
        pred = lensmodel.find_images(model, grid, *family_source(model, sel), z_s)
        ox, oy = model.to_frame(np.asarray(sel["ra"], float), np.asarray(sel["dec"], float))
        if len(pred):
            d = np.hypot(pred["x"][:, None] - ox[None, :], pred["y"][:, None] - oy[None, :])
            extra = np.flatnonzero(d.min(axis=1) > MATCH_ARCSEC)
        else:
            extra = np.array([], int)
        rows.append(
            (
                fam,
                z_s,
                len(sel),
                len(pred),
                len(extra),
                ",".join(f"{abs(float(pred['magnification'][k])):.1f}" for k in extra),
                ",".join(f"{float(pred['ra'][k]):.5f}:{float(pred['dec'][k]):.5f}" for k in extra),
            )
        )
    return Table(
        rows=rows or None,
        names=(
            "system",
            "z_used",
            "n_catalogued",
            "n_predicted",
            "n_extra",
            "extra_mu",
            "extra_radec",
        ),
        dtype=(str, float, int, int, int, str, str),
    )


# Predicted images of two models this close are one image (ASSUMPTION; 1.8x UNCOVER's 0.84" rms);
# the loose tier is ~1.3x UNCOVER's largest own-list residual, 2.31" (ASSUMPTION).
XMODEL_ARCSEC = 1.5
XMODEL_LOOSE_ARCSEC = 3.0
HALO_FILTER_ARCSEC = 1.0  # forced photometry on 1"-median-filtered stamps, as the Bergamini run


def wls_backtrace(model, bt: Table) -> Table:
    """``bt`` with every traced row's source replaced by its family's ``wls_source``, so that
    ``lens_consistency.predict_counter_images`` (which takes the mean source) solves from it."""
    out = bt.copy()
    systems = np.asarray(bt["system"]).astype(str)
    for fam in dict.fromkeys(systems):
        sel = np.flatnonzero((systems == fam) & np.isfinite(bt["beta_x"]))
        if len(sel) < 2:
            continue
        out["beta_x"][sel], out["beta_y"][sel] = family_source(model, bt[sel])
    return out


def forced_images(model, grid, bt: Table, catalog: Path, forced_image: str) -> Table:
    """UNCOVER's predicted images of the image list, classified and force-photometered as
    ``lens_consistency.py images --forced-image --halo-filter-arcsec 1.0`` does for Lenstool
    models (without its ``halo_recovery`` column), but from each family's WLS source."""
    wbt = wls_backtrace(model, bt)
    shapes = lc.load_shapes(catalog)
    table, unpredicted = lc.predict_counter_images(
        model, grid, wbt, shapes, MATCH_ARCSEC, depth=lc.depth_mag(shapes)
    )
    stamp, close = lc.image_stamper(forced_image, lc.MAX_SEARCH_ARCSEC + lc.FORCED_ANNULUS[1] + 0.1)
    try:
        lc.forced_check(
            table, wbt, lc.halo_stamper(stamp, HALO_FILTER_ARCSEC), unpredicted=unpredicted
        )
    finally:
        close()
    table.meta["forced"].update(halo_filter_arcsec=HALO_FILTER_ARCSEC, image=forced_image)
    table.meta["unpredicted"] = unpredicted
    return table


def check_forced_table(t: Table) -> None:
    """Refuse an ``images`` table that was not force-photometered with the same halo filter."""
    forced = t.meta.get("forced") or {}
    if "forced_class" not in t.colnames or forced.get("halo_filter_arcsec") != HALO_FILTER_ARCSEC:
        raise SystemExit(
            f"--crossmodel needs an images run with --forced-image and --halo-filter-arcsec "
            f"{HALO_FILTER_ARCSEC} (got halo filter {forced.get('halo_filter_arcsec')})"
        )


def _sep_arcsec(ra, dec, ra0, dec0) -> np.ndarray:
    """Small-angle separations (arcsec) of (ra, dec) arrays from one position."""
    cos = np.cos(np.deg2rad(dec0))
    return 3600.0 * np.hypot((np.asarray(ra, float) - ra0) * cos, np.asarray(dec, float) - dec0)


def _col(t: Table, name: str, k: int, default):
    if name not in t.colnames or np.ma.is_masked(t[name][k]):
        return default
    return type(default)(t[name][k])


def _nearest(t: Table, idx: np.ndarray, ra: float, dec: float) -> tuple[float, int]:
    """Separation to, and row of, the nearest of rows ``idx`` of ``t`` (NaN, -1 if none)."""
    if not len(idx):
        return np.nan, -1
    d = _sep_arcsec(t["ra"][idx], t["dec"][idx], ra, dec)
    return float(d.min()), int(idx[np.argmin(d)])


def crossmodel(tables: dict[str, Table]) -> Table:
    """D-041 cross-model rule for two models' ``images`` tables of the same image list. One row
    per extra image (``image_class`` != observed) of either model, with its own forced class, and
    of the other model's predictions of the same family:

    * the nearest *extra* image (``other_*``): ``both`` when within ``XMODEL_ARCSEC``,
      ``both_loose`` within ``XMODEL_LOOSE_ARCSEC`` (not one-to-one: two extras may share one);
    * the nearest image matched to a catalogued one (``other_observed_sep_arcsec``): an extra
      within ``XMODEL_LOOSE_ARCSEC`` of it is the other model's catalogued image, misplaced.

    ``other_image_class`` is "none" when the other model predicts no extra image of the family."""
    (na, ta), (nb, tb) = tables.items()
    rows = []
    for name, own, other_name, other in ((na, ta, nb, tb), (nb, tb, na, ta)):
        osys = np.asarray(other["system"]).astype(str)
        oobs = np.asarray(other["image_class"]).astype(str) == "observed"
        for i in np.flatnonzero(np.asarray(own["image_class"]).astype(str) != "observed"):
            r = own[i]
            ra, dec = float(r["ra"]), float(r["dec"])
            fam = osys == str(r["system"])
            sep, k = _nearest(other, np.flatnonzero(fam & ~oobs), ra, dec)
            sep_obs, _ = _nearest(other, np.flatnonzero(fam & oobs), ra, dec)
            rows.append(
                {
                    "model": name,
                    "system": str(r["system"]),
                    "ra": ra,
                    "dec": dec,
                    "magnification": float(r["magnification"]),
                    "image_class": str(r["image_class"]),
                    "forced_class": _col(own, "forced_class", i, ""),
                    "pred_snr": _col(own, "pred_snr", i, np.nan),
                    "flux_ratio": _col(own, "flux_ratio", i, np.nan),
                    "other_model": other_name,
                    "other_sep_arcsec": sep,
                    "other_magnification": _col(other, "magnification", k, np.nan)
                    if k >= 0
                    else np.nan,
                    "other_image_class": _col(other, "image_class", k, "") if k >= 0 else "none",
                    "other_forced_class": _col(other, "forced_class", k, "") if k >= 0 else "",
                    "other_observed_sep_arcsec": sep_obs,
                }
            )
    out = Table(rows=rows) if rows else Table()
    if len(out):
        sep = np.asarray(out["other_sep_arcsec"], float)
        out["both"] = sep <= XMODEL_ARCSEC
        out["both_loose"] = sep <= XMODEL_LOOSE_ARCSEC
    return out


def d041_verdict(r) -> str:
    """D-041 reading of one ``crossmodel`` row."""
    if float(r["other_observed_sep_arcsec"]) <= XMODEL_LOOSE_ARCSEC and not r["both"]:
        return "the other model's catalogued image, misplaced (ordinary)"
    if r["both_loose"]:
        if str(r["other_forced_class"]) == "recovered":
            return "the other model's image is recovered (position offset; ordinary)"
        return "both models" if r["both"] else "both models (loose tier)"
    return "model-dependent (D-041)"


def _r(v: float, nd: int = 2):
    v = float(v)
    return round(v, nd) if np.isfinite(v) else None


def crossmodel_summary(xm: Table, tables: dict[str, Table]) -> dict:
    """Per model: forced classes of all its predicted images, its extra images, how many the
    other model also predicts (both tiers), and every unmet prediction (``absent``/``ambiguous``)
    with its ``d041_verdict``."""
    out = {
        "xmodel_arcsec": XMODEL_ARCSEC,
        "xmodel_loose_arcsec": XMODEL_LOOSE_ARCSEC,
        "match_arcsec": MATCH_ARCSEC,
        "halo_filter_arcsec": HALO_FILTER_ARCSEC,
    }
    for m, t in tables.items():
        fc = np.array([_col(t, "forced_class", k, "") for k in range(len(t))])
        out[m] = {"forced_all": {c: int(np.sum(fc == c)) for c in sorted(set(fc) - {""})}}
    if not len(xm):
        return out
    for m in tables:
        t = xm[xm["model"] == m]
        fc = np.asarray(t["forced_class"]).astype(str)
        out[m].update(
            extra=len(t),
            both=int(np.sum(t["both"])),
            both_loose=int(np.sum(t["both_loose"])),
            extra_forced={
                c or "not_forced": {
                    "rows": int(np.sum(fc == c)),
                    "both": int(np.sum(t["both"][fc == c])),
                    "both_loose": int(np.sum(t["both_loose"][fc == c])),
                }
                for c in sorted(set(fc))
            },
        )
    unmet = xm[np.isin(np.asarray(xm["forced_class"]).astype(str), ["absent", "ambiguous"])]
    out["unmet"] = [
        {
            "model": str(r["model"]),
            "system": str(r["system"]),
            "ra_dec": [round(float(r["ra"]), 5), round(float(r["dec"]), 5)],
            "mu": _r(r["magnification"]),
            "forced_class": str(r["forced_class"]),
            "pred_snr": _r(r["pred_snr"], 1),
            "flux_ratio": _r(r["flux_ratio"]),
            "other_sep_arcsec": _r(r["other_sep_arcsec"]),
            "other_mu": _r(r["other_magnification"]),
            "other_forced_class": str(r["other_forced_class"]),
            "other_observed_sep_arcsec": _r(r["other_observed_sep_arcsec"]),
            "d041": d041_verdict(r),
        }
        for r in unmet
    ]
    out["unmet_in_both_models"] = sum(u["d041"].startswith("both models") for u in out["unmet"])
    return out


S3_I2D = "s3://stpubdata/jwst/public/jw02561/L3/t/o001/jw02561-o001_t003_nircam_clear-{}_i2d.fits"
BERG_42C = (3.57972, -30.40835)  # Bergamini+23 prediction of 4.2c (docs/fields/abell2744.md)
HIGHPASS_ARCSEC = 1.0  # median-filter width removing BCG/ICL light (ASSUMPTION, as halo_phot)
PEAK_RADIUS = {"ref": 0.4, "target": 1.0}  # peak search radius, arcsec (ASSUMPTION)
INJECT_SCALES = (1.0, 0.5, 0.25)


def _highpass(sci, xx, yy):
    return lc._median_subtract(sci, xx, yy, HIGHPASS_ARCSEC)


def _peak_flux(hp, err, xx, yy, radius: float, r_ap: float = 0.2) -> dict:
    inner = np.where((np.hypot(xx, yy) <= radius) & np.isfinite(hp), hp, -np.inf)
    k = int(np.argmax(inner))
    cx, cy = float(xx.flat[k]), float(yy.flat[k])
    ap = np.hypot(xx - cx, yy - cy) <= r_ap
    if not np.isfinite(inner).any():
        return {"dx_west": np.nan, "dy_north": np.nan, "flux": np.nan, "snr": np.nan}
    f, e = float(np.nansum(hp[ap])), float(np.sqrt(np.nansum(err[ap] ** 2)))
    return {"dx_west": cx, "dy_north": cy, "flux": f, "snr": f / e if e > 0 else np.nan}


def check_42c(model, grid, bt: Table) -> dict:
    """UNCOVER's 4.2 solution, high-pass peak photometry at its 4.2c and at 4.2a/b, and
    injections of the most magnified catalogued image's knot (scaled by the mu ratio) at 4.2c."""
    sel = bt[(np.asarray(bt["system"]).astype(str) == "4.2") & np.isfinite(bt["beta_x"])]
    if len(sel) < 2:
        return {"error": "family 4.2 has fewer than two traced images"}
    z = float(sel["z_used"][0])
    # the grid covers the whole cluster (4.2c lies ~44" from 4.2a/b)
    pred = lensmodel.find_images(model, grid, *family_source(model, sel), z)
    if len(pred) <= len(sel):
        return {"error": f"no extra image predicted ({len(pred)} images)"}
    ox, oy = model.to_frame(np.asarray(sel["ra"], float), np.asarray(sel["dec"], float))
    d = np.hypot(pred["x"][:, None] - ox[None, :], pred["y"][:, None] - oy[None, :])
    k_new = int(np.argmax(d.min(axis=1)))  # the predicted image farthest from a and b
    match = np.argmin(d, axis=0)  # predicted image of each catalogued one
    c = pred[k_new]
    bx, by = model.to_frame(*BERG_42C)
    out = {
        "z": z,
        "predicted": [
            {"ra": float(r["ra"]), "dec": float(r["dec"]), "mu": float(r["magnification"])}
            for r in pred
        ],
        "c_ra_dec": [float(c["ra"]), float(c["dec"])],
        "c_mu": float(c["magnification"]),
        "c_offset_from_bergamini_arcsec": float(np.hypot(c["x"] - bx, c["y"] - by)),
        "bands": {},
    }
    pos = {
        str(sel["image_id"][i]): (float(sel["ra"][i]), float(sel["dec"][i]))
        for i in range(len(sel))
    }
    mus = {str(sel["image_id"][i]): float(pred["magnification"][match[i]]) for i in range(len(sel))}
    out["mu_catalogued"] = mus
    ref_id = max(mus, key=lambda i: abs(mus[i]))
    for band in ("f150w", "f277w", "f444w"):
        stamp, close = lc.image_stamper(S3_I2D.format(band), 4.0)
        try:
            res = {}
            for iid, (ra, dec) in pos.items():
                sci, err, xx, yy = stamp(ra, dec)
                res[iid] = _peak_flux(_highpass(sci, xx, yy), err, xx, yy, PEAK_RADIUS["ref"])
            sci, err, xx, yy = stamp(float(c["ra"]), float(c["dec"]))
            res["4.2c_uncover"] = _peak_flux(
                _highpass(sci, xx, yy), err, xx, yy, PEAK_RADIUS["target"]
            )
            # injection: the reference image's positive high-pass knot (1.2" box) scaled to mu(c)
            rsci, _, rx, ry = stamp(*pos[ref_id])
            rhp = _highpass(rsci, rx, ry)
            j = np.unravel_index(
                np.argmax(np.where(np.hypot(rx, ry) < PEAK_RADIUS["ref"], rhp, -np.inf)), rhp.shape
            )
            n = int(round(0.6 / abs(float(rx[0, 1] - rx[0, 0]))))
            knot = np.clip(rhp[j[0] - n : j[0] + n + 1, j[1] - n : j[1] + n + 1], 0, None)
            ci, cj = np.array(sci.shape) // 2
            inj = {}
            for s in INJECT_SCALES:
                im = sci.copy()
                im[ci - n : ci + n + 1, cj - n : cj + n + 1] += (
                    s * knot * abs(out["c_mu"] / mus[ref_id])
                )
                inj[str(s)] = _peak_flux(_highpass(im, xx, yy), err, xx, yy, PEAK_RADIUS["target"])
            res["injections"] = inj
            res["expected_c_flux"] = res[ref_id]["flux"] * abs(out["c_mu"] / mus[ref_id])
            out["bands"][band] = res
        finally:
            close()
    out["reference"] = ref_id
    return out


def run(model, images: Table, zml: dict, label: str, out: Path) -> tuple[dict, object, Table]:
    bt = lensmodel.backtrace_images(model, images, zml)
    traced = images[np.isfinite(bt["beta_x"])]
    half = lc.grid_half_width(model, traced)
    grid = lensmodel.DeflectionGrid.compute(model, half, 0.25)
    ip = imageplane(model, grid, bt)
    d = np.asarray(ip["dtheta_arcsec"], float)
    ok = np.isfinite(d)
    mult = multiplicity(model, grid, bt)
    fam = np.asarray(bt["system"]).astype(str)
    traced = np.isfinite(bt["beta_x"])
    multi = np.array([np.sum(traced & (fam == f)) >= 2 for f in fam])
    for t in (ip, mult):
        t.meta.update(provenance="model_prediction", source=f"abell2744-uncover2 vs {label}")
    ip.write(out / f"imageplane_{label}.ecsv", overwrite=True)
    mult.write(out / f"multiplicity_{label}.ecsv", overwrite=True)
    return (
        {
            "image_list": label,
            "n_images": len(images),
            "n_traced": int(np.isfinite(bt["beta_x"]).sum()),
            "n_reimaged": int(ok.sum()),
            # traced images of families with >= 2 traced images that got no prediction
            "n_unsolved": int(((ip["n_predicted"] == 0) & np.isfinite(bt["beta_x"]) & multi).sum()),
            "n_shared_match": int(np.asarray(ip["shared_match"]).sum()),
            "rms_dtheta_arcsec": float(np.sqrt(np.mean(d[ok] ** 2))) if ok.any() else None,
            "median_dtheta_arcsec": float(np.median(d[ok])) if ok.any() else None,
            "p90_dtheta_arcsec": float(np.percentile(d[ok], 90)) if ok.any() else None,
            "max_dtheta_arcsec": float(d[ok].max()) if ok.any() else None,
            "n_over_2arcsec": int((d[ok] > 2.0).sum()),
            "families": len(mult),
            "families_with_extra": int((mult["n_extra"] > 0).sum()) if len(mult) else 0,
            "families_short": int((mult["n_predicted"] < mult["n_catalogued"]).sum())
            if len(mult)
            else 0,
        },
        grid,
        bt,
    )


def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=Path("results/lens_consistency/abell2744-uncover2"))
    ap.add_argument("--check-42c", action="store_true", help="4.2c photometry (S3 byte ranges)")
    ap.add_argument("--keep-maps", action="store_true", help="keep the cached deflection maps")
    ap.add_argument(
        "--crossmodel",
        type=Path,
        help="Bergamini+23 images_predicted.ecsv (lens_consistency images): D-041 cross-model rule",
    )
    ap.add_argument("--catalog", type=Path, help="F200W pipeline catalogue (with --crossmodel)")
    ap.add_argument("--forced-image", help="F277W _i2d (path or s3://; with --crossmodel)")
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    files = {k: fetch_catalog(url, sha) for k, (url, sha) in UNCOVER_V2.items()}
    model = uncover2_model(files)
    summary = {
        "model": "abell2744-uncover2",
        "model_sha256": model.sha256,
        "quoted_rms_arcsec": QUOTED_RMS,
    }
    own, _, _ = run(model, uncover_images(files["images"]), {}, "uncover2", args.out)
    # every traced multi-image family must re-image; the rms of the rest is not enough
    own["gate_passes"] = bool(
        own["rms_dtheta_arcsec"] is not None
        and own["n_unsolved"] == 0
        and own["rms_dtheta_arcsec"] <= lc.IMAGE_LIST_RMS_FACTOR * QUOTED_RMS
    )
    summary["own_list"] = own
    _, bfiles, par = lc.load_model("abell2744-bergamini23")
    berg = lensmodel.load_lenstool_images(bfiles["arcs.dat"])
    summary["bergamini23_list"], grid, bt = run(
        model, berg, par["z_m_limit"], "bergamini23", args.out
    )
    if args.crossmodel:
        if not (args.catalog and args.forced_image):
            raise SystemExit("--crossmodel needs --catalog and --forced-image")
        berg_forced = Table.read(args.crossmodel)
        check_forced_table(berg_forced)
        unc = forced_images(model, grid, bt, args.catalog, args.forced_image)
        unc.write(args.out / "images_predicted_uncover2.ecsv", overwrite=True)
        tables = {"bergamini23": berg_forced, "uncover2": unc}
        xm = crossmodel(tables)
        xm.meta.update(provenance="model_prediction", source=args.crossmodel.name)
        xm.write(args.out / "crossmodel.ecsv", overwrite=True)
        xsum = crossmodel_summary(xm, tables)
        (args.out / "crossmodel.json").write_text(json.dumps(xsum, indent=1) + "\n")
    if args.check_42c:
        summary["check_42c"] = check_42c(model, grid, bt)
    if not args.keep_maps:  # cloud-disk rule: ~250 MB of maps, deleted after use
        for key in ("alpha_x", "alpha_y"):
            files[key].unlink(missing_ok=True)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    return summary


if __name__ == "__main__":
    main()
