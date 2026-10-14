"""Vetting checks for the El Gordo and Sunrise follow-up sources (docs/fields/*.md).

Sources: El Gordo ``f200w_2828`` and ``f200w_1601``; Sunrise ``f200w_1756``, ``f200w_1907`` and
``f200w_1578``. For each (pipeline F200W catalogue position, observed):

* the DJA matched-aperture photometry (0.5" diameter) of the nearest DJA object, with its S/N;
* spike geometry: catalogue sources brighter than ``SPIKE_MAG`` on the same detector-frame row or
  column of the mosaic (within ``SPIKE_PIX`` pixels), where their diffraction spikes run;
* r = 0.25" aperture AB magnitudes in five bands from the level-3 mosaics (S3 byte ranges);
* the lens model's view: El Gordo with Caminha+23 (signed magnification and the counter-images a
  point source at the traced position would have, with predicted F200W magnitude and the nearest
  catalogue source, as ``a2744_followups``); Sunrise with the four RELICS v1 models (signed
  magnification from finite differences of local deflection patches, ``sunrise_relics``).

``--pm`` adds, for Sunrise ``f200w_1756``, the F115W centroid shift between o010 (2022-07-30) and
o120 (2023-01-10) relative to compact reference sources (proper motion).

Every number is ``derived`` (lens-model values: ``model_prediction``); thresholds are ASSUMPTIONs.

Usage: python scripts/field_followups.py [--out results/field_followups] [--pm]
"""

from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))

import a2744_followups as af  # noqa: E402
import lens_consistency as lc  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402
from jwst_anomaly.photometry import fetch_catalog  # noqa: E402

_MAST = "https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:JWST/product/{}_cat.ecsv"
_DJA = "https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/{}"
FIELDS = {
    "elgordo": {
        "product": "jw01176-o241_t012_nircam_clear-f200w",
        "cat_sha256": "47f60bddc1eb66ed2a276df35120354b049105fe8d9b85b1df12a29dc577b779",
        "i2d": "s3://stpubdata/jwst/public/jw01176/L3/t/o241/jw01176-o241_t012_nircam_clear-{}_i2d.fits",
        "dja": (
            _DJA.format("elgordo-grizli-v7.0-fix_phot.fits"),
            "ac1cf7064edc3f1a907c1b95b7eacd864e19810ff559e7dc39b917081ed66e31",
        ),
        "labels": (2828, 1601),
    },
    "sunrise": {
        "product": "jw02282-o010_t001_nircam_clear-f200w",
        "cat_sha256": "b06930655e770a6addf0c566b2b2399d039bbbb6be49ff9c1757b4183e4889df",
        "i2d": "s3://stpubdata/jwst/public/jw02282/L3/t/o010/jw02282-o010_t001_nircam_clear-{}_i2d.fits",
        "dja": (
            _DJA.format("sunrise-grizli-v7.5-fix_phot.fits"),
            "8bd178e94156d6a9a126e487dec8a850b6df9f6d76ceb9e36931d9960558ebf0",
        ),
        "labels": (1756, 1907, 1578),
    },
}
BANDS = ("f115w", "f150w", "f200w", "f277w", "f444w")
DJA_APER = 1  # ASEC_1 = 0.5" diameter, as configs/*.yaml matched_photometry
DJA_ZP = 23.9  # DJA fluxes are in uJy
DJA_MATCH_ARCSEC = 0.3  # as attach_photoz (ASSUMPTION)
SPIKE_MAG = 21.0  # aper50 AB; a source this bright has visible spikes (ASSUMPTION)
SPIKE_PIX = 6.0  # row/column half-width, mosaic pixels (~0.19"; ASSUMPTION)
SPIKE_MAX_ARCSEC = 30.0  # spike reach considered (ASSUMPTION)
AP_R = 0.25  # stamp aperture radius, arcsec (no aperture correction)
BG_R = (1.0, 2.0)  # stamp background annulus, arcsec
Z_GRID = (1.0, 2.0, 3.0, 6.0)  # source redshifts for the lens-model view (ASSUMPTION)
ELGORDO_MODEL = "elgordo-caminha23"
GRID_STEP = 0.25
GRID_MARGIN = 10.0
FORCED_BANDS = ("f115w", "f200w")  # 1601 is detected in the SW bands only (DJA)
FORCED_SEARCH = 1.5  # counter-image search radius, arcsec (as `lens_consistency images`)
FORCED_UIDS = ("f200w_1601",)  # 2828 is a spike (not a source); see the docs
SUNRISE_HALF = 50.0  # RELICS solver box half-width around the core, arcsec (sunrise_counter)
SUNRISE_COUNTER = ("f200w_1756",)  # inside the strong-lensing zone (|mu| > 8 at z = 2)
PM_EPOCHS = (
    "s3://stpubdata/jwst/public/jw02282/L3/t/o010/jw02282-o010_t001_nircam_clear-f115w_i2d.fits",
    "s3://stpubdata/jwst/public/jw02282/L3/t/o120/jw02282-o120_t001_nircam_clear-f115w_i2d.fits",
)
PM_DT_YR = (59954.54 - 59790.66) / 365.25  # F115W exposure MJDs (o120 - o010)
PM_CEN_R = 0.15  # centroid window radius, arcsec
PM_REF = (22.0, 26.0, 40.0, 20.0)  # reference sources: mag range, max distance ", min nn_dist pix


def dja_photometry(t: Table, ra: float, dec: float, bands=None) -> dict:
    """Nearest DJA object within ``DJA_MATCH_ARCSEC`` of (ra, dec): id, separation and per-band AB
    magnitude and S/N in the ``DJA_APER`` aperture (observed fluxes; magnitudes derived). Bands
    with no valid error are left out; non-positive fluxes have magnitude None. With no match,
    ``id`` is None."""
    cc = SkyCoord(np.asarray(t["ra"], float), np.asarray(t["dec"], float), unit="deg")
    d = SkyCoord(ra, dec, unit="deg").separation(cc).arcsec
    i = int(np.argmin(d))
    if d[i] > DJA_MATCH_ARCSEC:
        return {"id": None, "nearest_sep_arcsec": round(float(d[i]), 2), "bands": {}}
    suffix = f"_flux_aper_{DJA_APER}"
    if bands is None:
        bands = [c[: -len(suffix)] for c in t.colnames if c.endswith(suffix)]
    out = {"id": int(t["id"][i]), "sep_arcsec": round(float(d[i]), 2), "bands": {}}
    for b in bands:
        f = float(t[i][f"{b}_flux_aper_{DJA_APER}"])
        e = float(t[i][f"{b}_fluxerr_aper_{DJA_APER}"])
        if not (np.isfinite(f) and np.isfinite(e) and 0 < e < 1e10):
            continue
        mag = round(DJA_ZP - 2.5 * np.log10(f), 2) if f > 0 else None
        out["bands"][b] = {"mag": mag, "snr": round(f / e, 1)}
    return out


def spike_sources(cat: Table, label: int, pixel_arcsec: float) -> list[dict]:
    """Catalogue sources brighter than ``SPIKE_MAG`` (or with no finite magnitude: saturated
    cores) within ``SPIKE_MAX_ARCSEC`` whose mosaic row or column lies within ``SPIKE_PIX`` pixels
    of source ``label`` (a spike along the pixel axes). Stars outside the catalogue are missed."""
    r = cat[cat["label"] == label][0]
    dx = np.asarray(cat["xcentroid"], float) - float(r["xcentroid"])
    dy = np.asarray(cat["ycentroid"], float) - float(r["ycentroid"])
    sep = np.hypot(dx, dy) * pixel_arcsec
    mag = np.asarray(cat["aper50_abmag"], float)
    out = []
    for axis, d in (("row", dy), ("column", dx)):
        bright = ~np.isfinite(mag) | (mag < SPIKE_MAG)
        k = np.where((np.abs(d) <= SPIKE_PIX) & bright & (sep > 1.0))[0]
        k = k[sep[k] <= SPIKE_MAX_ARCSEC]
        for i in k[np.argsort(mag[k])]:
            out.append(
                {
                    "axis": axis,
                    "label": int(cat["label"][i]),
                    "aper50_abmag": round(float(mag[i]), 2) if np.isfinite(mag[i]) else None,
                    "sep_arcsec": round(float(sep[i]), 2),
                    "offset_pix": round(float(d[i]), 1),
                }
            )
    return out


def aperture_mag(sci, err, xx, yy) -> dict:
    """AB magnitude and S/N in an ``AP_R`` aperture on an i2d stamp (MJy/sr), median background."""
    pix = float(np.hypot(xx[0, 1] - xx[0, 0], yy[0, 1] - yy[0, 0]))
    r = np.hypot(xx, yy)
    bg = np.nanmedian(sci[(r > BG_R[0]) & (r < BG_R[1])])
    m = r <= AP_R
    to_jy = (pix / 206264.806) ** 2 * 1e6
    f = float(np.nansum(sci[m] - bg)) * to_jy
    e = float(np.sqrt(np.nansum(err[m] ** 2))) * to_jy
    if not e > 0:  # no valid pixels (off the mosaic or masked)
        return {"mag": None, "snr": None}
    return {"mag": round(-2.5 * np.log10(f / 3631.0), 2) if f > 0 else None, "snr": round(f / e, 1)}


def centroid(sci, xx, yy, r0: float = PM_CEN_R, iters: int = 3) -> tuple[float, float]:
    """Iterated first-moment centroid (arcsec, model-frame axes) within ``r0`` of the centre."""
    bg = np.nanmedian(sci[np.hypot(xx, yy) > 4 * r0])
    w0 = np.clip(np.nan_to_num(sci - bg), 0, None)
    x0 = y0 = 0.0
    for _ in range(iters):
        w = np.where(np.hypot(xx - x0, yy - y0) < r0, w0, 0.0)
        if w.sum() <= 0:
            return float("nan"), float("nan")
        x0, y0 = float((w * xx).sum() / w.sum()), float((w * yy).sum() / w.sum())
    return x0, y0


def relative_shift(target: np.ndarray, refs: np.ndarray) -> dict:
    """Shift of the target (row 0 of epoch-2 minus epoch-1 centroids, mas) relative to the median
    reference shift, with the references' robust scatter."""
    ok = np.isfinite(refs).all(axis=1)
    ref = refs[ok]
    med = np.median(ref, axis=0)
    mad = 1.4826 * np.median(np.abs(ref - med), axis=0)
    d = target - med
    return {
        "n_ref": int(ok.sum()),
        "frame_shift_mas": np.round(med, 1).tolist(),
        "ref_scatter_mas": np.round(mad, 1).tolist(),
        "shift_mas": np.round(d, 1).tolist(),
        "shift_sigma": np.round(d / mad, 1).tolist(),
    }


def proper_motion(cat: Table, label: int) -> dict:
    """F115W centroid shift of ``label`` between the two Sunrise epochs (S3 byte ranges)."""
    c = cat["sky_centroid"]
    r0 = cat[cat["label"] == label][0]["sky_centroid"]
    sep = c.separation(r0).arcsec
    mag = np.asarray(cat["aper50_abmag"], float)
    lo, hi, rmax, nn = PM_REF
    k = np.where(
        (sep > 1) & (sep < rmax) & (mag > lo) & (mag < hi) & (np.asarray(cat["nn_dist"]) > nn)
    )[0]
    pts = [(r0.ra.deg, r0.dec.deg)] + [(c[i].ra.deg, c[i].dec.deg) for i in k]
    cens = []
    for uri in PM_EPOCHS:
        stamp, close = lc.image_stamper(uri, 0.8)
        try:
            cens.append(
                np.array([centroid(s, xx, yy) for s, _, xx, yy in (stamp(*p) for p in pts)])
            )
        finally:
            close()
    d = (cens[1] - cens[0]) * 1000.0
    out = relative_shift(d[0], d[1:])
    out["dt_yr"] = round(PM_DT_YR, 3)
    # per-axis error: the references' scatter (conservative for a bright target)
    out["pm_mas_per_yr"] = round(float(np.hypot(*out["shift_mas"])) / PM_DT_YR, 1)
    out["pm_err_mas_per_yr"] = round(float(np.hypot(*out["ref_scatter_mas"])) / PM_DT_YR, 1)
    out["pm_2sigma_upper_mas_per_yr"] = round(
        out["pm_mas_per_yr"] + 2 * out["pm_err_mas_per_yr"], 1
    )
    return out


def sunrise_mu(ra: float, dec: float) -> dict:
    """Signed magnification per RELICS model at ``Z_GRID`` (model_prediction)."""
    import requests
    from astropy.wcs import WCS
    from sunrise_relics import (
        BASE,
        HALF,
        KEYS,
        MODELS,
        frame_check,
        local_lens,
        read_header,
        read_patch,
    )
    from sunrise_zfit import dratio

    def one(model):
        s = requests.Session()
        urls = {k: BASE.format(m=model, k=KEYS[k]) for k in ("defx", "defy")}
        heads = {k: read_header(s, u) for k, u in urls.items()}
        wcs = WCS(heads["defx"][0])
        step = frame_check(wcs)
        i, j = (float(v) for v in wcs.world_to_pixel_values(ra, dec))
        i0, j0 = int(round(i)), int(round(j))
        try:
            p = {k: read_patch(s, urls[k], *heads[k], i0, j0) for k in urls}
        except RuntimeError:
            return model, None  # outside the model's map
        c, r = HALF + i - i0, HALF + j - j0
        return model, {
            str(z): round(local_lens(p["defx"], p["defy"], step, c, r, dratio(z))["mu"], 3)
            for z in Z_GRID
        }

    with ThreadPoolExecutor(len(MODELS)) as ex:
        return dict(ex.map(one, MODELS))


def drop_self(cps: list[dict], v: dict) -> list[dict]:
    """A counter-image matched to the source itself is not a detection of the counter-image: its
    match is removed and it is flagged ``near_source`` (forced photometry would see the source)."""
    for cp in cps:
        cp["near_source"] = cp["sep_from_source_arcsec"] <= FORCED_SEARCH + lc.FORCED_R_AP
        if cp["f200w_match"] == v["label"]:
            cp["f200w_match"] = cp["f200w_match_mag"] = None
    return cps


def sunrise_counter_images(rows: dict, uids, cat: Table) -> dict:
    """RELICS v1 counter-images (model_prediction) of ``uids`` at ``Z_GRID``, each model's maps
    over a ``SUNRISE_HALF`` box (``sunrise_counter.load_model``), solved on a ``GRID_STEP`` grid."""
    import sunrise_counter as sc

    cat_coords = SkyCoord(cat["ra"], cat["dec"], unit="deg")
    with ThreadPoolExecutor(len(sc.MODELS)) as ex:
        loaded = ex.map(lambda m: sc.load_model(m, SUNRISE_HALF), sc.MODELS)
        models = dict(zip(sc.MODELS, loaded, strict=True))
    out = {uid: {} for uid in uids}
    for name, model in models.items():
        grid = lensmodel.DeflectionGrid.cached(
            model, lc.grid_cache_path(model, SUNRISE_HALF, GRID_STEP), SUNRISE_HALF, GRID_STEP
        )
        for uid in uids:
            v = rows[uid]
            x, y = model.to_frame(np.array([v["ra"]]), np.array([v["dec"]]))
            c = SkyCoord(v["ra"], v["dec"], unit="deg")
            per_z = []
            for z in Z_GRID:
                mu, cps = af.counter_images(
                    model, grid, float(x[0]), float(y[0]), c, z, v["mag_f200w"], cat, cat_coords
                )
                per_z.append(
                    {"z": z, "mu_source": round(mu, 3), "counter_images": drop_self(cps, v)}
                )
            out[uid][name] = per_z
    return out


def elgordo_lens(rows: dict, cat: Table) -> dict:
    """Caminha+23 signed magnification and counter-images (model_prediction) at ``Z_GRID``."""
    model, _, _ = lc.load_model(ELGORDO_MODEL)
    lc.apply_frame_offset(ELGORDO_MODEL, model)
    cat_coords = SkyCoord(cat["ra"], cat["dec"], unit="deg")
    ras = np.array([v["ra"] for v in rows.values()])
    decs = np.array([v["dec"] for v in rows.values()])
    xs, ys = model.to_frame(ras, decs)
    half = float(np.max(np.hypot(xs, ys))) + GRID_MARGIN
    grid = lensmodel.DeflectionGrid.cached(
        model, lc.grid_cache_path(model, half, GRID_STEP), half, GRID_STEP
    )
    out = {}
    for (uid, v), x, y in zip(rows.items(), xs, ys, strict=True):
        c = SkyCoord(v["ra"], v["dec"], unit="deg")
        per_z = []
        for z in Z_GRID:
            mu, cps = af.counter_images(model, grid, x, y, c, z, v["mag_f200w"], cat, cat_coords)
            per_z.append({"z": z, "mu_source": round(mu, 3), "counter_images": drop_self(cps, v)})
        out[uid] = {"r_model_arcsec": round(float(np.hypot(x, y)), 1), "per_z": per_z}
    return out


def forced_counter_images(i2d: str, v: dict, lens: dict, bands=FORCED_BANDS) -> None:
    """D-029 forced photometry (``sunrise_counter.classify``, with its injection) at every
    predicted counter-image in ``lens["per_z"]``, in place. The predicted flux is the source's own
    aperture flux in the same band times |mu_c / mu_source| (point-source scaling; ASSUMPTION)."""
    import sunrise_counter as sc

    for band in bands:
        stamp, close = lc.image_stamper(
            i2d.format(band), lc.MAX_SEARCH_ARCSEC + lc.FORCED_ANNULUS[1] + 0.1
        )
        try:
            sci, err, xx, yy = stamp(v["ra"], v["dec"])
            f0, e0 = lc.aperture_snr(sci, err, xx, yy, 0.0, 0.0)
            v.setdefault("forced_self", {})[band] = {"flux": f0, "snr": round(f0 / e0, 1)}
            for row in lens["per_z"]:
                for cp in row["counter_images"]:
                    if cp["near_source"]:
                        cp.setdefault("forced", {})[band] = {"class": "near_source"}
                        continue
                    pred = f0 * abs(cp["mu"] / row["mu_source"])
                    res = sc.classify(stamp, cp["ra"], cp["dec"], pred, FORCED_SEARCH)
                    cp.setdefault("forced", {})[band] = {
                        k: (round(x, 3) if isinstance(x, float) else x) for k, x in res.items()
                    }
        finally:
            close()


def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=Path("results/field_followups"))
    ap.add_argument("--pm", action="store_true", help="Sunrise 1756 proper motion (S3 ranges)")
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    out = {
        "provenance": "positions and photometry observed; magnitudes, shifts derived; "
        "magnifications and counter-images model_prediction; thresholds ASSUMPTIONs",
        "sources": {},
    }
    for field, spec in FIELDS.items():
        cat_path = fetch_catalog(_MAST.format(spec["product"]), spec["cat_sha256"])
        cat = Table.read(cat_path)
        shapes = lc.load_shapes(cat_path)
        dja = Table.read(fetch_catalog(*spec["dja"]))
        rows = {}
        for label in spec["labels"]:
            r = cat[cat["label"] == label][0]
            rows[f"f200w_{label}"] = {
                "label": label,
                "ra": round(float(r["sky_centroid"].ra.deg), 6),
                "dec": round(float(r["sky_centroid"].dec.deg), 6),
                "mag_f200w": round(float(r["aper50_abmag"]), 2),
                "ci_50_30": round(float(r["CI_50_30"]), 2),
                "is_extended": bool(r["is_extended"]),
                "iso_snr": round(float(r["isophotal_flux"] / r["isophotal_flux_err"]), 1),
            }

        def band_stamps(band, rows=rows, spec=spec):
            stamp, close = lc.image_stamper(spec["i2d"].format(band), 2.5)
            try:
                return {(uid, band): stamp(v["ra"], v["dec"]) for uid, v in rows.items()}
            finally:
                close()

        stamps = {}
        with ThreadPoolExecutor(len(BANDS)) as ex:  # I/O-bound S3 range reads
            for part in ex.map(band_stamps, BANDS):
                stamps.update(part)
        for uid, v in rows.items():
            sci, err, xx, yy = stamps[(uid, "f200w")]
            pix = float(np.hypot(xx[0, 1] - xx[0, 0], yy[0, 1] - yy[0, 0]))
            v["dja"] = dja_photometry(dja, v["ra"], v["dec"])
            v["spike_sources"] = spike_sources(cat, v["label"], pix)
            v["stamp_ab"] = {b: aperture_mag(*stamps[(uid, b)]) for b in BANDS}
        if field == "elgordo":
            for uid, lens in elgordo_lens(rows, shapes).items():
                rows[uid]["lens_caminha23"] = lens
                if uid in FORCED_UIDS:
                    forced_counter_images(spec["i2d"], rows[uid], lens)
        else:
            for v in rows.values():
                v["lens_relics"] = sunrise_mu(v["ra"], v["dec"])
            for uid, per_model in sunrise_counter_images(rows, SUNRISE_COUNTER, shapes).items():
                rows[uid]["counter_images_relics"] = per_model
                for per_z in per_model.values():
                    forced_counter_images(spec["i2d"], rows[uid], {"per_z": per_z})
            if args.pm:
                rows["f200w_1756"]["proper_motion"] = proper_motion(cat, 1756)
        for uid, v in rows.items():
            out["sources"][f"{field}:{uid}"] = v
    (args.out / "summary.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    main()
