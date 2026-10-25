"""SN Eos 101.2 / 101.1 difference-image flux ratio vs the VENUS AstroLensPy v2.1 model.

Frozen plan: ``results/transients/sneos_ratio/PLAN.md``. This module holds the plan's rules as
functions (resampling, registration check, aperture ratio, injection, decision); the run driver
(empirical-PSF star, S3 stamps of the 9 mosaics, contact sheet) is the next step. R_MODEL is a
model_prediction from ``results/transients/sneos/summary.json``.
"""

from __future__ import annotations

import numpy as np
from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry
from scipy.ndimage import map_coordinates
from scipy.signal import fftconvolve

R_MODEL = 1.162  # |μ_101.2 / μ_101.1|, sneos/summary.json (PLAN.md)
APER_ARCSEC = 0.20  # ASSUMPTION (PLAN.md rule 2)
ANNULUS_ARCSEC = (0.35, 0.60)  # PLAN.md rule 2
MAX_SHIFT_ARCSEC = 0.03  # ASSUMPTION (PLAN.md rule 1)
MIN_SNR = 5.0  # ASSUMPTION (PLAN.md rule 3)
INJ_TOL = 0.10  # ASSUMPTION (PLAN.md rule 4)
INJ_RATIOS = (1.16, 0.70)  # PLAN.md rule 4
FLOOR_LN = 0.30  # ASSUMPTION (PLAN.md rule 6)
SN_MASK_ARCSEC = 0.4  # PLAN.md rule 1


def resample(src: np.ndarray, src_wcs, dst_wcs, shape: tuple[int, int]) -> np.ndarray:
    """``src`` (surface brightness) bilinearly resampled onto the ``dst_wcs`` pixel grid."""
    jj, ii = np.indices(shape, dtype=float)
    ra, dec = dst_wcs.pixel_to_world_values(ii, jj)
    si, sj = src_wcs.world_to_pixel_values(ra, dec)
    return map_coordinates(src, [sj, si], order=1, cval=np.nan)


def register_offset(a: np.ndarray, b: np.ndarray, mask: np.ndarray) -> tuple[float, float]:
    """Sub-pixel (dx, dy) of ``b`` relative to ``a`` from the cross-correlation peak (parabolic
    refinement); pixels in ``mask`` and non-finite pixels are set to zero after mean removal."""

    def prep(x):
        x = np.where(mask | ~np.isfinite(x), np.nan, x)
        x = x - np.nanmean(x)
        return np.nan_to_num(x)

    pa, pb = prep(a), prep(b)
    cc = fftconvolve(pb, pa[::-1, ::-1], mode="same")
    j, i = np.unravel_index(np.argmax(cc), cc.shape)

    def refine(m, c, p):
        den = m - 2 * c + p
        return 0.0 if den == 0 else 0.5 * (m - p) / den

    di = refine(cc[j, i - 1], cc[j, i], cc[j, i + 1]) if 0 < i < cc.shape[1] - 1 else 0.0
    dj = refine(cc[j - 1, i], cc[j, i], cc[j + 1, i]) if 0 < j < cc.shape[0] - 1 else 0.0
    cy, cx = cc.shape[0] // 2, cc.shape[1] // 2  # zero lag of a "same" correlation
    return float(i + di - cx), float(j + dj - cy)


def aperture_pair(
    diff: np.ndarray, err: np.ndarray, p1: tuple[float, float], p2: tuple[float, float], pix: float
) -> dict:
    """PLAN.md rule 2: background-subtracted aperture fluxes at ``p1`` (101.1) and ``p2``
    (101.2), pixel positions (x, y); the annulus excludes the other image's aperture."""
    r = APER_ARCSEC / pix
    rin, rout = (a / pix for a in ANNULUS_ARCSEC)
    yy, xx = np.indices(diff.shape)
    out = {}
    for name, p, q in (("f1", p1, p2), ("f2", p2, p1)):
        other = np.hypot(xx - q[0], yy - q[1]) <= r
        ann = CircularAnnulus(p, rin, rout).to_mask(method="center").to_image(diff.shape) > 0
        sel = ann & ~other & np.isfinite(diff)
        bkg = float(np.median(diff[sel]))
        ap = CircularAperture(p, r)
        bad = ~np.isfinite(diff) | ~np.isfinite(err)
        t = aperture_photometry(diff - bkg, ap, error=np.where(bad, 0, err), mask=bad)
        out[name] = float(t["aperture_sum"][0])
        out["s" + name] = float(t["aperture_sum_err"][0])
    f1, f2, s1, s2 = out["f1"], out["f2"], out["sf1"], out["sf2"]
    ratio = f2 / f1
    out["snr1"], out["snr2"] = f1 / s1 if s1 else np.inf, f2 / s2 if s2 else np.inf
    out["ratio"] = ratio
    out["sratio"] = abs(ratio) * float(np.hypot(s1 / f1, s2 / f2))
    out["detected"] = bool(out["snr1"] >= MIN_SNR and out["snr2"] >= MIN_SNR)
    return out


def inject_pair(
    img: np.ndarray, psf: np.ndarray, p1: tuple[float, float], sep: tuple[float, float], f1, ratio
) -> np.ndarray:
    """``img`` plus a unit-sum ``psf`` scaled to ``f1`` at ``p1`` and ``ratio * f1`` at
    ``p1 + sep`` (pixels; integer-pixel placement after a sub-pixel bilinear shift)."""
    out = np.array(img, float, copy=True)
    k = psf / psf.sum()
    h0, h1 = k.shape[0] // 2, k.shape[1] // 2
    for (x, y), f in (
        (p1, f1),
        ((p1[0] + sep[0], p1[1] + sep[1]), ratio * f1),
    ):
        ix, iy = int(np.floor(x)), int(np.floor(y))
        fx, fy = x - ix, y - iy
        jj, ii = np.indices(k.shape, dtype=float)
        ks = map_coordinates(k, [jj - fy, ii - fx], order=1, cval=0.0)
        y0, x0 = iy - h0, ix - h1  # clip the PSF where it runs off the stamp
        ya, xa = max(y0, 0), max(x0, 0)
        yb, xb = min(y0 + k.shape[0], out.shape[0]), min(x0 + k.shape[1], out.shape[1])
        out[ya:yb, xa:xb] += f * ks[ya - y0 : yb - y0, xa - x0 : xb - x0]
    return out


def recovered(measured: float, injected: float) -> bool:
    """PLAN.md rule 4: ratio recovered within INJ_TOL (relative)."""
    return bool(abs(measured / injected - 1.0) <= INJ_TOL)


def decide(ln_r: float, sig_ln: float, r_model: float = R_MODEL) -> str:
    """PLAN.md rule 6: 'differs' when |Δ| > 3σ and > FLOOR_LN; otherwise 'agrees'."""
    if not np.isfinite(ln_r) or not np.isfinite(sig_ln):
        return "undetermined"
    d = abs(ln_r - np.log(r_model))
    return "differs" if (d > 3.0 * sig_ln and d > FLOOR_LN) else "agrees"


def combine(ln_r: np.ndarray, sig: np.ndarray) -> tuple[float, float]:
    """PLAN.md rule 5: inverse-variance mean of ln R; σ = max(propagated, scatter / √n)."""
    ln_r, sig = np.asarray(ln_r, float), np.asarray(sig, float)
    if ln_r.size == 0:
        return float("nan"), float("nan")
    w = 1.0 / sig**2
    mean = float(np.sum(w * ln_r) / np.sum(w))
    prop = float(1.0 / np.sqrt(np.sum(w)))
    scat = float(np.std(ln_r, ddof=1) / np.sqrt(ln_r.size)) if ln_r.size > 1 else 0.0
    return mean, max(prop, scat)


# --- run driver -------------------------------------------------------------------------------

EARLY = "jw06882-o038_t041"
LATER = ("jw09493-o001_t001", "jw12299-o001_t001")
BANDS = ("f277w", "f356w", "f444w")
URI = "s3://stpubdata/jwst/public/jw{p}/L3/t/o{o}/{stem}_nircam_clear-{band}_i2d.fits"
STAMP_HALF_ARCSEC = 3.0  # rule 1: 6" x 6"
PSF_HALF_ARCSEC = 1.0
INJ_RADIUS_ARCSEC = 1.8  # implementation of rule 4: 3 positions at 0/120/240 deg from the midpoint
GAIA_RADIUS_ARCMIN = 2.0
GAIA_ISOLATION_ARCSEC = 3.0  # implementation of rule 4 ("isolated")
GAIA_PM_SNR = 3.0  # implementation of rule 4: a star has a significant parallax or proper motion
OUT = "results/transients/sneos_ratio"
IMAGES = "results/lens_consistency/macs1931_inputs/m1931_images.ecsv"


def uri(stem: str, band: str) -> str:
    p, o = stem[2:7], stem[9:12]
    return URI.format(p=p, o=o, stem=stem, band=band)


class Mosaic:
    """SCI / ERR sections of one ``_i2d`` by S3 byte range."""

    def __init__(self, path: str):
        import warnings

        import fsspec
        from astropy.io import fits
        from astropy.wcs import WCS, FITSFixedWarning

        fs, p = fsspec.core.url_to_fs(path, anon=True)
        self.fo = fs.open(p, "rb", block_size=2**20, cache_type="readahead")
        self.hdul = fits.open(self.fo, lazy_load_hdus=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FITSFixedWarning)
            self.wcs = WCS(self.hdul["SCI"].header)
        self.pix = float(abs(self.wcs.proj_plane_pixel_scales()[0].to_value("arcsec")))

    def section(self, ra: float, dec: float, half_arcsec: float):
        """(sci, err, wcs of the section) centred on (ra, dec)."""
        x, y = (int(round(float(v))) for v in self.wcs.world_to_pixel_values(ra, dec))
        h = int(np.ceil(half_arcsec / self.pix))
        ny, nx = self.hdul["SCI"].header["NAXIS2"], self.hdul["SCI"].header["NAXIS1"]
        if not (h <= x < nx - h and h <= y < ny - h):
            raise ValueError(f"section at ({ra}, {dec}) runs off the mosaic")
        sl = (slice(y - h, y + h + 1), slice(x - h, x + h + 1))
        sci = np.asarray(self.hdul["SCI"].section[sl], float)
        err = np.asarray(self.hdul["ERR"].section[sl], float)
        return sci, err, self.wcs[sl]

    def close(self):
        self.hdul.close()
        self.fo.close()


def gaia_psf_star(ra: float, dec: float) -> dict:
    """Rule 4 (implementation): nearest Gaia DR3 source with a > 3σ parallax or proper motion and no
    other Gaia source within GAIA_ISOLATION_ARCSEC (VizieR I/355/gaiadr3)."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.vizier import Vizier

    c = SkyCoord(ra, dec, unit="deg")
    v = Vizier(
        columns=[
            "RA_ICRS",
            "DE_ICRS",
            "Source",
            "Gmag",
            "Plx",
            "e_Plx",
            "pmRA",
            "e_pmRA",
            "pmDE",
            "e_pmDE",
        ],
        row_limit=-1,
    )
    t = v.query_region(c, radius=GAIA_RADIUS_ARCMIN * u.arcmin, catalog="I/355/gaiadr3")[0]
    sc = SkyCoord(t["RA_ICRS"], t["DE_ICRS"], unit="deg")
    sig = np.nan_to_num(np.abs(t["Plx"] / t["e_Plx"])) > GAIA_PM_SNR
    sig |= np.nan_to_num(np.abs(t["pmRA"] / t["e_pmRA"])) > GAIA_PM_SNR
    sig |= np.nan_to_num(np.abs(t["pmDE"] / t["e_pmDE"])) > GAIA_PM_SNR
    order = np.argsort(sc.separation(c).arcsec)
    for k in order:
        if not sig[k]:
            continue
        sep = sc.separation(sc[k]).arcsec
        if np.sum(sep < GAIA_ISOLATION_ARCSEC) == 1:
            return {
                "source_id": int(t["Source"][k]),
                "ra": float(t["RA_ICRS"][k]),
                "dec": float(t["DE_ICRS"][k]),
                "gmag": float(t["Gmag"][k]),
                "dist_arcsec": float(sc[k].separation(c).arcsec),
            }
    raise RuntimeError("no isolated Gaia star in the cone")


def psf_stamp(m: Mosaic, star: dict) -> np.ndarray | None:
    """Border-median-subtracted star stamp, or None when its core is not finite (saturated)."""
    sci, _, _ = m.section(star["ra"], star["dec"], PSF_HALF_ARCSEC)
    c = sci.shape[0] // 2
    if not np.all(np.isfinite(sci[c - 2 : c + 3, c - 2 : c + 3])):
        return None
    border = np.concatenate([sci[0], sci[-1], sci[:, 0], sci[:, -1]])
    s = np.nan_to_num(sci - np.nanmedian(border))
    return s if s.sum() > 0 else None


def injection_trials(diff, err, psf, p1, p2, f1, pix) -> list[dict]:
    """Rule 4: both ratios at 3 positions INJ_RADIUS_ARCSEC from the pair midpoint."""
    mid = ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)
    sep = (p2[0] - p1[0], p2[1] - p1[1])
    r = INJ_RADIUS_ARCSEC / pix
    out = []
    for ang in (0.0, 120.0, 240.0):
        q1 = (
            mid[0] + r * np.cos(np.radians(ang)) - sep[0] / 2,
            mid[1] + r * np.sin(np.radians(ang)) - sep[1] / 2,
        )
        q2 = (q1[0] + sep[0], q1[1] + sep[1])
        for ratio in INJ_RATIOS:
            img = inject_pair(diff, psf, q1, sep, f1, ratio)
            base = aperture_pair(diff, err, q1, q2, pix)
            got = aperture_pair(img, err, q1, q2, pix)
            meas = (got["f2"] - base["f2"]) / (got["f1"] - base["f1"])
            out.append(
                {
                    "angle_deg": ang,
                    "injected": ratio,
                    "measured": float(meas),
                    "recovered": recovered(meas, ratio),
                }
            )
    return out


def run(out: str = OUT) -> dict:
    import json
    from pathlib import Path

    from astropy.table import Table

    t = Table.read(IMAGES, format="ascii.ecsv")
    t = t[np.isin(np.asarray(t["image_id"], str), ["101.1", "101.2"])]
    t.sort("image_id")
    (ra1, ra2), (de1, de2) = list(map(float, t["ra"])), list(map(float, t["dec"]))
    ram, dem = (ra1 + ra2) / 2, (de1 + de2) / 2
    star = gaia_psf_star(ram, dem)
    rows, sheets = [], {}
    for band in BANDS:
        e = Mosaic(uri(EARLY, band))
        sci_e, err_e, w = e.section(ram, dem, STAMP_HALF_ARCSEC)
        pix = e.pix
        psf = psf_stamp(e, star)
        e.close()
        p1 = tuple(float(v) for v in w.world_to_pixel_values(ra1, de1))
        p2 = tuple(float(v) for v in w.world_to_pixel_values(ra2, de2))
        yy, xx = np.indices(sci_e.shape)
        sn = (np.hypot(xx - p1[0], yy - p1[1]) <= SN_MASK_ARCSEC / pix) | (
            np.hypot(xx - p2[0], yy - p2[1]) <= SN_MASK_ARCSEC / pix
        )
        for stem in LATER:
            m = Mosaic(uri(stem, band))
            sci_l, err_l, wl = m.section(ram, dem, STAMP_HALF_ARCSEC + 1.0)
            m.close()
            res = resample(sci_l, wl, w, sci_e.shape)
            rerr = resample(err_l, wl, w, sci_e.shape)
            dx, dy = register_offset(sci_e, res, sn)
            shift = float(np.hypot(dx, dy) * pix)
            diff = sci_e - res
            derr = np.hypot(err_e, rerr)
            row = {
                "band": band,
                "later": stem,
                "pix_arcsec": pix,
                "shift_arcsec": shift,
                "registered": shift <= MAX_SHIFT_ARCSEC,
            }
            ph = aperture_pair(diff, derr, p1, p2, pix)
            if ph["f1"] < 0 and ph["f2"] < 0:  # rule 1: sign taken from the data
                diff = -diff
                ph = aperture_pair(diff, derr, p1, p2, pix)
                row["sign_flipped"] = True
            row.update(ph)
            if psf is not None and row["registered"]:
                f_inj = ph["f1"] if ph["f1"] > 0 else MIN_SNR * ph["sf1"]
                tr = injection_trials(diff, derr, psf, p1, p2, f_inj, pix)
                row["injections"] = tr
                row["gate_I"] = sum(x["recovered"] for x in tr) >= 5
            else:
                row["gate_I"] = False
            row["used"] = bool(row["registered"] and row["detected"] and row["gate_I"])
            rows.append(row)
            sheets[(band, stem)] = (sci_e, res, diff, p1, p2)
    used = [r for r in rows if r["used"] and r["ratio"] > 0]
    ln_r = np.log([r["ratio"] for r in used]) if used else np.array([])
    sig = np.array([r["sratio"] / r["ratio"] for r in used]) if used else np.array([])
    mean, s = combine(ln_r, sig)
    summary = {
        "plan": f"{OUT}/PLAN.md",
        "r_model": R_MODEL,
        "psf_star": star,
        "rows": rows,
        "n_used": len(used),
        "ln_r_obs": mean,
        "sigma_ln": s,
        "r_obs": float(np.exp(mean)) if np.isfinite(mean) else None,
        "decision": decide(mean, s),
        "provenance": {"fluxes": "observed (difference images)", "r_model": "model_prediction"},
    }
    Path(out).mkdir(parents=True, exist_ok=True)
    Path(out, "summary.json").write_text(json.dumps(summary, indent=1, default=float) + "\n")
    contact_sheet(sheets, Path(out, "contact_sheet.png"))
    return summary


def contact_sheet(sheets: dict, path) -> None:
    """Early, resampled later and difference stamps per band and pair, SN apertures circled."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(len(sheets), 3, figsize=(6, 2 * len(sheets)))
    for row, ((band, stem), (*ims, p1, p2)) in zip(
        np.atleast_2d(axes), sheets.items(), strict=True
    ):
        for ax, im, lab in zip(row, ims, ("6882 o038", stem[:8], "diff"), strict=True):
            v = np.nanpercentile(np.abs(ims[2]), 99.5)
            ax.imshow(im, origin="lower", cmap="gray", vmin=-v, vmax=v if lab == "diff" else 3 * v)
            for p in (p1, p2):
                ax.plot(*p, "o", mfc="none", mec="r", ms=8)
            ax.set_title(f"{band} {lab}", fontsize=7)
            ax.set_xticks([])
            ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(path, dpi=60)
    plt.close(fig)


if __name__ == "__main__":
    import json as _json

    s = run()
    print(_json.dumps({k: s[k] for k in ("n_used", "r_obs", "sigma_ln", "decision")}))
    for r in s["rows"]:
        print(
            {
                k: r[k]
                for k in (
                    "band",
                    "later",
                    "shift_arcsec",
                    "snr1",
                    "snr2",
                    "ratio",
                    "sratio",
                    "gate_I",
                    "used",
                )
            }
        )
