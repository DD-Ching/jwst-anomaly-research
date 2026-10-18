"""Forced aperture photometry at fixed sky positions in HST drizzled association products.

The HST counterpart of ``hit_lightcurve.py`` (which reads JWST level-3 mosaics). For each
association (e.g. ``jeou05010``) the pipeline ``_drc`` (ACS, WFC3/UVIS; CTE-corrected) or ``_drz``
(WFC3/IR) is read by S3 byte range (a box around the target; nothing stored). Fluxes are e-/s x
PHOTFNU in nJy, from ``hit_lightcurve.aperture_photometry`` (aperture minus annulus median); the
error is the robust std of the same aperture at random box positions away from the target (local,
includes confusion and correlated noise). Compact field sources (local maxima above ``--peak-snr``)
are measured too, for an epoch-to-epoch zero point. Every number is ``derived``; radii are
ASSUMPTIONs (plan: ``results/transients/macs1423_h000_hst/PLAN.md``).

    python scripts/hst_forced.py --ra 215.949099 --dec 24.086948 --assoc jeou05010 ifaa47010 \\
        --out results/transients/macs1423_h000_hst
"""

from __future__ import annotations

import argparse
import json
import socket
import sys
import zlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hit_lightcurve import aperture_photometry  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

APER = 0.2  # arcsec (ASSUMPTION; #323 HST convention)
ANNULUS = (0.6, 1.0)
BOX = 8.0  # half-size of the read box, arcsec
N_RANDOM = 300
EXCLUDE = 1.5  # random apertures at least this far from the target, arcsec
MAX_FIELD = 30


def resolve(assocs: list[str]) -> dict[str, str]:
    """S3 URI of each association's drizzled product (``_drc`` preferred over ``_drz``)."""
    from astroquery.mast import Observations

    from jwst_anomaly import query

    Observations.TIMEOUT = 120
    obs = Observations.query_criteria(obs_id=assocs)
    prods = query.list_products(obs, subgroups=("DRC", "DRZ"), calib_level=3, cloud_uris=True)
    out: dict[str, str] = {}
    for f, u in zip(prods["productFilename"], prods["cloud_uri"], strict=True):
        a, kind = str(f).split("_")[0], str(f).split("_")[1]
        if u and (a not in out or kind == "drc.fits"):
            out[a] = str(u)
    return out


def centroid(img: np.ndarray, x: int, y: int) -> tuple[float, float]:
    """Flux-weighted centroid in the 3 x 3 pixels around a peak (positive pixels only)."""
    sub = np.clip(np.nan_to_num(img[y - 1 : y + 2, x - 1 : x + 2]), 0, None)
    if sub.shape != (3, 3) or sub.sum() <= 0:
        return float(x), float(y)
    yy, xx = np.mgrid[-1:2, -1:2]
    return x + float((sub * xx).sum() / sub.sum()), y + float((sub * yy).sum() / sub.sum())


def find_peaks(img: np.ndarray, scale: float, peak_snr: float, n_max: int):
    """Centroids (x, y) of the ``n_max`` brightest local maxima above ``peak_snr`` x the noise of
    the 1-px-smoothed image."""
    sm = ndimage.gaussian_filter(np.nan_to_num(img), 1.0)
    pix_sd = 1.4826 * np.nanmedian(np.abs(sm - np.nanmedian(sm)))
    size = int(round(0.5 / scale)) | 1
    peaks = (sm == ndimage.maximum_filter(sm, size=size)) & (sm > peak_snr * pix_sd)
    py, px = np.nonzero(peaks)
    order = np.argsort(-sm[py, px])[:n_max]
    return [centroid(img, int(a), int(b)) for a, b in zip(px[order], py[order], strict=True)]


def measure(
    assoc: str, uri: str, ra: float, dec: float, peak_snr: float, offset=(0.0, 0.0)
) -> dict:
    import fsspec

    socket.setdefaulttimeout(60)
    fs, path = fsspec.core.url_to_fs(uri, anon=True)
    with fs.open(path, "rb", block_size=2**20, cache_type="readahead") as fo:
        h = fits.open(fo, lazy_load_hdus=True)
        p0 = h[0].header
        sci_h = h["SCI"].header
        w = WCS(sci_h, h)
        scale = float(np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))) * 3600
        # registration: the target appears at (ra, dec) + offset (arcsec) in this visit's WCS
        ra_m = ra + offset[0] / 3600.0 / np.cos(np.radians(dec))
        dec_m = dec + offset[1] / 3600.0
        x, y = (float(v) for v in w.world_to_pixel_values(ra_m, dec_m))
        hb = int(np.ceil(BOX / scale))
        ny, nx = sci_h["NAXIS2"], sci_h["NAXIS1"]
        x0, y0 = max(int(x) - hb, 0), max(int(y) - hb, 0)
        x1, y1 = min(int(x) + hb + 1, nx), min(int(y) + hb + 1, ny)
        if x1 - x0 < hb or y1 - y0 < hb:
            return {"assoc": assoc, "error": "target outside the image"}
        sci = np.asarray(h["SCI"].section[y0:y1, x0:x1], float)
        photfnu = sci_h.get("PHOTFNU", p0.get("PHOTFNU"))
        if photfnu is None:  # ACS: from PHOTFLAM (erg/s/cm2/A per e-/s) and PHOTPLAM (A), in Jy
            flam = sci_h.get("PHOTFLAM", p0.get("PHOTFLAM"))
            plam = sci_h.get("PHOTPLAM", p0.get("PHOTPLAM"))
            photfnu = 3.33564095e4 * float(flam) * float(plam) ** 2
        rec = {
            "assoc": assoc,
            "instrument": f"{p0.get('INSTRUME')}/{p0.get('DETECTOR')}",
            "filter": p0.get("FILTER") or p0.get("FILTER1") or "",
            "mjd": float(p0.get("EXPSTART")),
            "exptime": float(p0.get("EXPTIME")),
            "pixscale": scale,
            "photfnu": photfnu,
            "bunit": sci_h.get("BUNIT"),
        }
    if "CLEAR" in rec["filter"]:
        rec["filter"] = (
            p0.get("FILTER2") if "CLEAR" in str(p0.get("FILTER1")) else p0.get("FILTER1")
        )
    k = float(photfnu) * 1e9  # e-/s -> nJy
    xt, yt = x - x0, y - y0
    r_px, ann = APER / scale, (ANNULUS[0] / scale, ANNULUS[1] / scale)
    f = aperture_photometry(sci, xt, yt, r_px, ann) * k
    rng = np.random.default_rng(zlib.crc32(assoc.encode()))
    m = ann[1] + 1
    rnd = []
    for _ in range(N_RANDOM * 3):
        xx, yy = rng.uniform(m, sci.shape[1] - m), rng.uniform(m, sci.shape[0] - m)
        if np.hypot(xx - xt, yy - yt) * scale < EXCLUDE:
            continue
        q = aperture_photometry(sci, xx, yy, r_px, ann)
        if np.isfinite(q):
            rnd.append(q * k)
        if len(rnd) >= N_RANDOM:
            break
    rnd = np.asarray(rnd)
    err = float(1.4826 * np.median(np.abs(rnd - np.median(rnd)))) if len(rnd) > 30 else np.nan
    rec.update(flux_njy=f, err_njy=err, n_random=int(len(rnd)))
    rec["offset_arcsec"] = list(offset)
    pk = [
        p
        for p in find_peaks(sci, scale, peak_snr, MAX_FIELD + 5)
        if np.hypot(p[0] - xt, p[1] - yt) * scale >= 1.0
    ]
    px = np.array([p[0] for p in pk[:MAX_FIELD]])
    py = np.array([p[1] for p in pk[:MAX_FIELD]])
    rd = w.pixel_to_world_values(px + x0, py + y0)
    rec["field"] = [
        {
            "ra": float(r),
            "dec": float(d),
            "flux_njy": float(aperture_photometry(sci, float(a), float(b), r_px, ann) * k),
        }
        for a, b, r, d in zip(px, py, rd[0], rd[1], strict=True)
    ]
    h0 = int(np.ceil(1.2 / scale))
    xi, yi = int(round(xt)), int(round(yt))
    rec["stamp"] = (sci[yi - h0 : yi + h0 + 1, xi - h0 : xi + h0 + 1] * k).tolist()
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ra", type=float, required=True)
    ap.add_argument("--dec", type=float, required=True)
    ap.add_argument("--assoc", nargs="+", required=True)
    ap.add_argument("--peak-snr", type=float, default=20.0)
    ap.add_argument(
        "--offsets", type=Path, default=None, help="JSON {assoc: [dRA cos(dec), dDec] arcsec}"
    )
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    uris = resolve(a.assoc)
    offs = json.loads(a.offsets.read_text()) if a.offsets else {}
    missing = sorted(set(a.assoc) - set(uris))

    def job(assoc):
        try:
            off = tuple(offs.get(assoc, (0.0, 0.0)))
            return measure(assoc, uris[assoc], a.ra, a.dec, a.peak_snr, off)
        except Exception as exc:  # unreadable or missing product: recorded, never a pass
            return {"assoc": assoc, "error": str(exc)[:200]}

    with ThreadPoolExecutor(a.threads) as ex:
        rows = list(ex.map(job, sorted(uris)))
    rows += [{"assoc": m, "error": "no drizzled product found"} for m in missing]
    a.out.mkdir(parents=True, exist_ok=True)
    stamps = {r["assoc"]: r.pop("stamp") for r in rows if "stamp" in r}
    np.savez_compressed(
        a.out / "stamps.npz", **{k: np.asarray(v, np.float32) for k, v in stamps.items()}
    )
    (a.out / "hst_forced.json").write_text(
        json.dumps({"label": "derived", "uris": uris, "rows": rows}, indent=1, default=float) + "\n"
    )
    for r in rows:
        print({k: v for k, v in r.items() if k != "field"}, len(r.get("field", [])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
