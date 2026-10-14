"""Lens-model checks for the five unvetted Abell 2744 follow-up sources (docs/fields/abell2744.md).

For each source (pipeline catalogue position, observed):

* the nearest published multiple image in Bergamini+23 ``arcs.dat`` and UNCOVER v2.0 (observed);
* the DJA eazy photo-z (z16, z50, z84; derived);
* in the Bergamini+23 model (model_prediction): the signed magnification at the photo-z
  percentiles, and the images a point source at the traced source position would have. Each
  predicted counter-image (> ``SELF_ARCSEC`` from the source) gets a predicted F200W magnitude
  (the source's magnitude scaled by |mu_cp / mu_src|) and the nearest F200W catalogue source.

A source whose model predicts a detectable counter-image with no catalogue source near it would be
a lens-model inconsistency to vet; a source with no counter-image is singly imaged and the lens
model says nothing more about it. Every threshold is an ASSUMPTION.

Usage: python scripts/a2744_followups.py --catalog <f200w_cat.ecsv> --catalog-f277w <f277w_cat.ecsv>
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))

import a2744_uncover as au  # noqa: E402
import lens_consistency as lc  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402
from jwst_anomaly.photometry import fetch_catalog  # noqa: E402

# (band, pipeline label) of jw02561-o001_t003_nircam_clear-<band>_cat.ecsv
SOURCES = {
    "f200w_5904": ("f200w", 5904),
    "f200w_7987": ("f200w", 7987),
    "f200w_4731": ("f200w", 4731),
    "f200w_7298": ("f200w", 7298),
    "f277w_264": ("f277w", 264),
}
# DJA v7.2 eazy zout (docs/fields/abell2744.md "Lens-model and exotic screens")
ZOUT = (
    "https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/abell2744clu-grizli-v7.2-fix.eazypy.zout.fits",
    "436626861c7dbf370419a8787dd79140fad99c2ecbd268257941afe6d493991c",
)
SELF_ARCSEC = 0.3  # a predicted image this close to the source is the source itself (ASSUMPTION)
MATCH_ARCSEC = 1.5  # counter-image search radius, as `lens_consistency images` (ASSUMPTION)
GRID_STEP = 0.1  # solver grid step, arcsec (as `images`)
GRID_MARGIN = 40.0  # grid half-width = farthest source from the model centre + this (ASSUMPTION)
# f277w_264 knot test: F200W peaks of a 1-px-smoothed 2.4" stamp above its 98.5th percentile,
# r = 0.08" apertures, background = stamp median beyond 0.9" (all ASSUMPTIONs)
I2D = "s3://stpubdata/jwst/public/jw02561/L3/t/o001/jw02561-o001_t003_nircam_clear-{}_i2d.fits"
KNOT_BANDS = ("f115w", "f150w", "f200w", "f277w", "f356w", "f410m", "f444w")
KNOT_R_AP = 0.08
KNOT_BG_R = 0.9
KNOT_PCTL = 98.5


def nearest(c: SkyCoord, t: Table, idcol: str) -> dict:
    cc = SkyCoord(np.asarray(t["ra"], float), np.asarray(t["dec"], float), unit="deg")
    sep = c.separation(cc).arcsec
    k = int(np.argmin(sep))
    return {"id": str(t[idcol][k]), "sep_arcsec": round(float(sep[k]), 2)}


def counter_images(model, grid, x, y, c, z, mag, cat, cat_coords):
    """Signed magnification of a source at model-frame ``(x, y)`` (sky ``c``) for redshift ``z``,
    and the other images of its traced source position (model_prediction), each with a predicted
    magnitude and the nearest catalogue source."""
    lm = model.lens_map(np.array([x]), np.array([y]), np.array([z]))
    mu = 1.0 / float(lm["a11"][0] * lm["a22"][0] - lm["a12"][0] * lm["a21"][0])
    pred = lensmodel.find_images(model, grid, float(lm["beta_x"][0]), float(lm["beta_y"][0]), z)
    cps = []
    for p in pred:
        cp = SkyCoord(float(p["ra"]), float(p["dec"]), unit="deg")
        if cp.separation(c).arcsec <= SELF_ARCSEC:
            continue
        d = cp.separation(cat_coords).arcsec
        k = int(np.argmin(d))
        cps.append(
            {
                "ra": round(float(p["ra"]), 6),
                "dec": round(float(p["dec"]), 6),
                "sep_from_source_arcsec": round(float(cp.separation(c).arcsec), 2),
                "mu": round(float(p["magnification"]), 3),
                "pred_mag_f200w": round(
                    mag - 2.5 * np.log10(abs(float(p["magnification"]) / mu)), 2
                ),
                "nearest_f200w": int(cat["label"][k]),
                "nearest_f200w_arcsec": round(float(d[k]), 2),
                "nearest_f200w_mag": round(float(cat["mag"][k]), 2),
            }
        )
    return mu, cps


def knot_photometry(stamps: dict, ref: str = "f200w") -> dict:
    """Aperture fluxes (stamp units x pixels, observed) of the ``ref``-band peaks in every band.

    ``stamps[band] = (sci, err, xx, yy)`` as :func:`lens_consistency.image_stamper` returns."""
    from scipy import ndimage

    sci, _, xx, yy = stamps[ref]
    sm = ndimage.gaussian_filter(np.nan_to_num(sci), 1)
    peaks = np.argwhere(
        (sm == ndimage.maximum_filter(sm, size=5)) & (sm > np.nanpercentile(sm, KNOT_PCTL))
    )
    knots = [(float(xx[j, i]), float(yy[j, i])) for j, i in peaks]
    out = {"knots_arcsec": [[round(a, 3), round(b, 3)] for a, b in knots], "flux": {}}
    for band, (sci, err, xx, yy) in stamps.items():
        bg = float(np.nanmedian(sci[np.hypot(xx, yy) > KNOT_BG_R]))
        row = []
        for kx, ky in knots:
            m = np.hypot(xx - kx, yy - ky) <= KNOT_R_AP
            row.append(
                [
                    round(float(np.nansum(sci[m] - bg)), 4),
                    round(float(np.sqrt(np.nansum(err[m] ** 2))), 4),
                ]
            )
        out["flux"][band] = row
    return out


def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--catalog", type=Path, required=True, help="F200W pipeline _cat.ecsv")
    ap.add_argument("--catalog-f277w", type=Path, required=True, help="F277W pipeline _cat.ecsv")
    ap.add_argument("--out", type=Path, default=Path("results/a2744_followups"))
    ap.add_argument("--knots", action="store_true", help="f277w_264 knot photometry (S3 ranges)")
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    shapes = {"f200w": lc.load_shapes(args.catalog), "f277w": lc.load_shapes(args.catalog_f277w)}
    zout = fetch_catalog(*ZOUT)
    for s in shapes.values():
        lc.attach_photoz(s, zout)
    f200 = shapes["f200w"]
    cat200 = SkyCoord(f200["ra"], f200["dec"], unit="deg")

    model, files, par = lc.load_model("abell2744-bergamini23")
    berg = lensmodel.load_lenstool_images(files["arcs.dat"])
    unc = au.uncover_images(fetch_catalog(*au.UNCOVER_V2["images"]))

    rows = {}
    for uid, (band, label) in SOURCES.items():
        t = shapes[band]
        r = t[t["label"] == label][0]
        rows[uid] = r
    xs, ys = model.to_frame(
        np.array([float(r["ra"]) for r in rows.values()]),
        np.array([float(r["dec"]) for r in rows.values()]),
    )
    half = float(np.max(np.hypot(xs, ys))) + GRID_MARGIN
    grid = lensmodel.DeflectionGrid.cached(
        model, lc.grid_cache_path(model, half, GRID_STEP), half, GRID_STEP
    )

    out = {
        "model": "abell2744-bergamini23",
        "grid_half_width_arcsec": round(half, 1),
        "sources": {},
    }
    for (uid, r), x, y in zip(rows.items(), xs, ys, strict=True):
        c = SkyCoord(float(r["ra"]), float(r["dec"]), unit="deg")
        rec = {
            "ra": round(float(r["ra"]), 6),
            "dec": round(float(r["dec"]), 6),
            "mag_iso": round(float(r["mag"]), 2),
            "is_extended": bool(r["is_extended"]),
            "nearest_bergamini23_image": nearest(c, berg, "image_id"),
            "nearest_uncover2_image": nearest(c, unc, "image_id"),
            "z_phot": [round(float(r[k]), 3) for k in ("z160", "z_phot", "z840")],
            "per_z": [],
        }
        for zkey in ("z160", "z_phot", "z840"):
            z = float(r[zkey])
            if not np.isfinite(z) or z <= model.z_lens + 0.02:
                rec["per_z"].append(
                    {"z": None if not np.isfinite(z) else z, "note": "not behind the cluster"}
                )
                continue
            mu, cps = counter_images(model, grid, x, y, c, z, float(r["mag"]), f200, cat200)
            rec["per_z"].append(
                {"z": round(z, 3), "mu_source": round(mu, 3), "counter_images": cps}
            )
        out["sources"][uid] = rec
    if args.knots:
        r = rows["f277w_264"]
        stamps = {}
        for band in KNOT_BANDS:
            stamp, close = lc.image_stamper(I2D.format(band), 1.2)
            try:
                stamps[band] = stamp(float(r["ra"]), float(r["dec"]))
            finally:
                close()
        out["f277w_264_knots"] = knot_photometry(stamps)
    (args.out / "summary.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    main()
