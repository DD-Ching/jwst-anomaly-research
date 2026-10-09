"""Archival HST check of the D-064 untestable rejected quasar pairs (PSF-subtracted residuals).

For each pair with HAP imaging (MAST hapcut, F814W preferred), fit two point sources with a shared
Moffat profile plus a constant background, then measure the residual flux inside the circle that
has the two images as its diameter (beyond ``Params.image_mask_arcsec`` from each image, after
removing each image's halo residual). The D-064 control lenses with
HST imaging go through the same chain first: the test is trusted only if it shows their lens
galaxies. Thresholds are ASSUMPTIONs (``Params``); outputs are ``derived``.

    python scripts/w12_niq_hst.py --out outputs/w12_niq_hst
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.modeling import fitting, models
from astropy.table import Table
from astropy.wcs import WCS

from jwst_anomaly import paths, schema

UNTESTABLE = [
    "J0130+0725",
    "J0728+2607",
    "J0941-2443",
    "J1428+0500",
    "J2355-4553",
    "J092718.37+211357.4",
    "J124257.32+254303.0",
    "J090437.41+113426.2",
    "J094235.04+231028.9",
    "J132405.19+282331.9",
    "J171101.70+292950.9",
]
CONTROLS = ["J0628-7448", "J1550+0221", "J2308+3201", "J151538.59+151135.8", "J132236.41+105239.4"]


@dataclass(frozen=True)
class Params:
    size_px: int = 250  # hapcut cutout (0.04"/px skycells: 10")
    search_arcsec: float = 3.5  # the companion lies within the catalogued separation + margin
    # residual excludes the image cores (PSF-core mismatch; 0.15" left core residuals in J0728)
    image_mask_arcsec: float = 0.3
    detect_sigma: float = 5.0  # residual flux / its noise for "residual >= 5 sigma"
    sep_match_arcsec: float = 0.5  # the fitted pair must be the catalogued pair (as D-064)


def pick_image(files, c: SkyCoord, radius_arcsec: float) -> Path | None:
    """The single-skycell F814W cutout with the most valid pixels within ``radius_arcsec`` of the
    target (a target near a skycell edge appears in two cutouts). The "combined_skycells" product
    is not used: its WCS does not describe its pixels (separations of 10^4" in a first run)."""
    best, best_n = None, 0
    for f in map(Path, files):
        if "f814w" not in f.name or "combined" in f.name or "coarse" in f.name:
            continue
        with fits.open(f) as h:
            hdu = next(x for x in h if x.data is not None and x.data.ndim == 2)
            w, data = WCS(hdu.header), np.nan_to_num(hdu.data.astype(float))
        pix = float(np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix))) * 3600)
        x0, y0 = w.world_to_pixel(c)
        yy, xx = np.indices(data.shape)
        n = int(((np.hypot(xx - x0, yy - y0) <= radius_arcsec / pix) & (data != 0)).sum())
        if n > best_n:
            best, best_n = f, n
    return best


def two_psf_fit(img: np.ndarray, peaks: list[tuple[float, float]]):
    """Two Moffat2D point sources with a shared shape (gamma, alpha) and a constant background."""
    yy, xx = np.indices(img.shape)
    bg = float(np.median(img))
    m = models.Const2D(bg)
    for x0, y0 in peaks:
        amp = float(img[int(round(y0)), int(round(x0))] - bg)
        m = m + models.Moffat2D(
            amp,
            x0,
            y0,
            gamma=2.0,
            alpha=2.5,
            bounds={"gamma": (0.5, 10.0), "alpha": (1.2, 8.0), "amplitude": (0, None)},
        )
    m.gamma_2.tied = lambda mm: mm.gamma_1
    m.alpha_2.tied = lambda mm: mm.alpha_1
    fit = fitting.LevMarLSQFitter()(m, xx, yy, img, maxiter=500)
    return fit, img - fit(xx, yy)


def find_two_peaks(
    img: np.ndarray, centre: tuple[float, float], radius_px: float, min_sep_px: float
):
    """Brightest pixel within ``radius_px`` of ``centre`` and the brightest one at least
    ``min_sep_px`` from it (the two quasar images)."""
    yy, xx = np.indices(img.shape)
    sm = img.copy()
    within = np.hypot(xx - centre[0], yy - centre[1]) <= radius_px
    sm[~within] = -np.inf
    y1, x1 = np.unravel_index(np.argmax(sm), sm.shape)
    sm[np.hypot(xx - x1, yy - y1) < min_sep_px] = -np.inf
    y2, x2 = np.unravel_index(np.argmax(sm), sm.shape)
    return [(float(x1), float(y1)), (float(x2), float(y2))]


def remove_halos(res: np.ndarray, fit, nbin: float = 1.0) -> np.ndarray:
    """Subtract, around each image, the radial profile of the residual measured on the side facing
    away from the other image (PSF-model mismatch is radially symmetric about each image; light
    between the images is not). Profiles out to the pair separation, in ``nbin``-pixel rings."""
    out = res.copy()
    pos = [(fit.x_0_1.value, fit.y_0_1.value), (fit.x_0_2.value, fit.y_0_2.value)]
    yy, xx = np.indices(res.shape)
    sep = np.hypot(pos[1][0] - pos[0][0], pos[1][1] - pos[0][1])
    for k, (x, y) in enumerate(pos):
        ox, oy = pos[1 - k]
        r = np.hypot(xx - x, yy - y)
        cosang = ((xx - x) * (ox - x) + (yy - y) * (oy - y)) / np.maximum(r * sep, 1e-9)
        away = cosang < 0  # the half-plane facing away from the other image
        near = r < sep
        rbin = (r / nbin).astype(int)
        prof = np.zeros(rbin[near].max() + 2)
        for b in np.unique(rbin[near]):
            sel = near & away & (rbin == b)
            prof[b] = np.median(res[sel]) if sel.any() else 0.0
        out[near] -= prof[rbin[near]]
    return out


def residual_between(
    res: np.ndarray, fit, pix_arcsec: float, p: Params, noise: float, valid=None
) -> dict:
    """Residual flux in the pair-diameter circle, excluding the image cores, and its S/N."""
    (x1, y1), (x2, y2) = (fit.x_0_1.value, fit.y_0_1.value), (fit.x_0_2.value, fit.y_0_2.value)
    xc, yc = (x1 + x2) / 2, (y1 + y2) / 2
    r = np.hypot(x2 - x1, y2 - y1) / 2
    yy, xx = np.indices(res.shape)
    core = p.image_mask_arcsec / pix_arcsec
    m = (np.hypot(xx - xc, yy - yc) <= r) & (np.hypot(xx - x1, yy - y1) > core)
    m &= np.hypot(xx - x2, yy - y2) > core
    if valid is not None:
        m &= valid  # missing or edge pixels are not data
    flux = float(res[m].sum())
    err = float(noise * np.sqrt(m.sum()))
    return {
        "sep_arcsec": 2 * r * pix_arcsec,
        "resid_flux": flux,
        "resid_err": err,
        "resid_snr": flux / err if err > 0 else np.nan,
        "n_pix": int(m.sum()),
        "flux_ratio_images": float(fit.amplitude_1.value / fit.amplitude_2.value),
    }


def run(args) -> None:
    from astroquery.mast import Hapcut

    p = Params()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    systems = Table.read(paths.repo_root() / "results" / "w12_niq" / "systems.ecsv")
    rows, panels, inj = [], [], []
    for name in CONTROLS + UNTESTABLE:
        r = systems[systems["name"] == name][0]
        d = out / name
        d.mkdir(exist_ok=True)
        c = SkyCoord(float(r["ra"]), float(r["dec"]), unit="deg")
        try:
            files = list(d.glob("*.fits")) or list(
                Hapcut.download_cutouts(c, size=[p.size_px, p.size_px], path=str(d))["Local Path"]
            )
            f = pick_image(files, c, p.search_arcsec)
        except Exception as e:  # no HAP footprint, service error: recorded, never a pass
            files, f = [], None
            print(f"{name}: no HAP cutout ({type(e).__name__})")
        row = {
            "name": name,
            "group": str(r["group"]),
            "sep_cat": float(r["sep_cat"]),
            "hst_image": f.name if f else "",
            "status": "no HST F814W HAP cutout",
        }
        if f is not None:
            with fits.open(f) as h:
                hdu = next(x for x in h if x.data is not None and x.data.ndim == 2)
                data = hdu.data.astype(float)
                w = WCS(hdu.header)
                zp = ab_zeropoint(hdu.header, h[0].header)
            valid_all = np.isfinite(data) & (data != 0)  # 0 = outside the drizzled footprint
            corr = noise_correlation(data, valid_all)
            img = np.where(valid_all, data, np.median(data[valid_all]))
            pix = float(np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix))) * 3600)
            x0, y0 = w.world_to_pixel(c)
            peaks = find_two_peaks(img, (x0, y0), p.search_arcsec / pix, 0.5 * r["sep_cat"] / pix)
            # fit on a stamp around the pair (the full cutout holds unrelated sources)
            xm, ym = np.mean([q[0] for q in peaks]), np.mean([q[1] for q in peaks])
            half = int(1.5 * r["sep_cat"] / pix) + 10
            y0s, x0s = max(int(ym) - half, 0), max(int(xm) - half, 0)
            stamp = (slice(y0s, int(ym) + half), slice(x0s, int(xm) + half))
            img, valid = img[stamp], valid_all[stamp]
            peaks = [(x - x0s, y - y0s) for x, y in peaks]
            fit, res = two_psf_fit(img, peaks)
            # aperture noise: pixel MAD x the drizzle correlation factor
            noise = 1.4826 * np.median(np.abs(res[valid] - np.median(res[valid]))) * corr
            clean = remove_halos(res, fit)
            row.update(residual_between(clean, fit, pix, p, noise, valid), pixel_arcsec=pix)
            row["noise_corr"] = corr
            row["resid_snr_raw"] = residual_between(res, fit, pix, p, noise, valid)["resid_snr"]
            snr = row["resid_snr"]
            # the S/N states a residual only; lens light or core mismatch needs the stamp
            if abs(row["sep_arcsec"] - row["sep_cat"]) > p.sep_match_arcsec:
                row["status"] = "pair mismatch"
            elif not np.isfinite(snr) or snr <= -p.detect_sigma:
                row["status"] = "fit inconclusive (over-subtraction)"
            else:
                row["status"] = "residual >= 5 sigma" if snr >= p.detect_sigma else "no residual"
            row["zp_ab"] = zp
            row["mag_resid"] = (
                zp - 2.5 * np.log10(row["resid_flux"]) if row["resid_flux"] > 0 else np.nan
            )
            row["mag_5sigma"] = zp - 2.5 * np.log10(p.detect_sigma * row["resid_err"])
            row["mag_limit_inject"] = np.nan
            if row["status"] == "no residual":
                rec = inject_recovery(img, peaks, fit, pix, p, zp, noise, valid)
                for mag, s_inj, frac in rec:
                    inj.append(
                        {"name": name, "inject_mag": mag, "recovered_snr": s_inj,
                         "snr_above_baseline": s_inj - snr, "flux_frac": frac}
                    )  # fmt: skip
                ok = [mag for mag, s_inj, _ in rec if s_inj - snr >= p.detect_sigma]
                row["mag_limit_inject"] = max(ok) if ok else np.nan
            panels.append((name, row, img, res, clean, fit))
        rows.append(row)
        print(row)
    t = Table(rows=rows)
    t.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source="scripts/w12_niq_hst.py on MAST HAP cutouts (hapcut)",
        params=asdict(p),
    )
    res_dir = paths.repo_root() / "results" / "w12_niq"
    t.write(res_dir / "hst_residuals.ecsv", overwrite=True)
    ti = Table(rows=inj) if inj else Table(names=("name", "inject_mag"))
    ti.meta.update(
        provenance=schema.Provenance.SIMULATED.value,
        source="Sersic n=4, r_eff=0.3 arcsec, ellip 0.2, at the SIS-predicted lens position, "
        "injected into the real stamp and run through the same fit (scripts/w12_niq_hst.py)",
    )
    ti.write(res_dir / "hst_injections.ecsv", overwrite=True)
    (res_dir / "hst_params.json").write_text(json.dumps(asdict(p), indent=2) + "\n")
    sheet(panels, res_dir / "hst_residuals.jpg")


def noise_correlation(img: np.ndarray, valid: np.ndarray, block: int = 5) -> float:
    """Factor by which drizzle-correlated noise inflates an aperture sum: std of block sums of
    source-free background pixels over sqrt(block²) x the pixel std (1 for white noise)."""
    good = valid & np.isfinite(img)
    med = np.median(img[good])
    sd = 1.4826 * np.median(np.abs(img[good] - med))
    bg = good & (np.abs(img - med) < 3 * sd)
    n0, n1 = (s // block * block for s in img.shape)
    a = np.where(bg, img - med, np.nan)[:n0, :n1].reshape(n0 // block, block, n1 // block, block)
    full = np.isfinite(a).all(axis=(1, 3))
    sums = np.nansum(a, axis=(1, 3))[full]
    if len(sums) < 20:
        return np.nan
    sd_sum = 1.4826 * np.median(np.abs(sums - np.median(sums)))
    return float(sd_sum / (block * sd))


def ab_zeropoint(hdr, primary) -> float:
    """AB zero point for an image in electrons/s from PHOTFLAM and PHOTPLAM (STScI convention)."""
    pf = hdr.get("PHOTFLAM", primary.get("PHOTFLAM"))
    pl = hdr.get("PHOTPLAM", primary.get("PHOTPLAM"))
    return float(-2.5 * np.log10(pf) - 5 * np.log10(pl) - 2.408) if pf and pl else np.nan


def inject_recovery(
    img, peaks, fit, pix, p: Params, zp: float, noise: float, valid, mags=(21, 22, 23, 24, 25)
):
    """Inject an early-type lens galaxy (Sersic n = 4, r_eff 0.3") at the SIS-predicted position
    (on the line between the images, at sep * f / (1 + f) from the fainter one, f = faint/bright)
    into the real stamp and re-measure with the same chain. Returns (mag, S/N, flux fraction)."""
    (x1, y1, a1), (x2, y2, a2) = (
        (fit.x_0_1.value, fit.y_0_1.value, fit.amplitude_1.value),
        (fit.x_0_2.value, fit.y_0_2.value, fit.amplitude_2.value),
    )
    if a1 < a2:
        (x1, y1, a1), (x2, y2, a2) = (x2, y2, a2), (x1, y1, a1)
    q = a2 / a1
    xl, yl = x2 + (x1 - x2) * q / (1 + q), y2 + (y1 - y2) * q / (1 + q)
    yy, xx = np.indices(img.shape)
    g = models.Sersic2D(1, 0.3 / pix, 4, x_0=xl, y_0=yl, ellip=0.2, theta=0.5)(xx, yy)
    g /= g.sum()
    out = []
    for mag in mags:
        flux = 10 ** (-0.4 * (mag - zp))
        f2, res = two_psf_fit(img + flux * g, peaks)
        m = residual_between(remove_halos(res, f2), f2, pix, p, noise, valid)
        out.append((mag, m["resid_snr"], m["resid_flux"] / flux))
    return out


def sheet(panels, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image

    n = len(panels)
    fig, axes = plt.subplots(n, 3, figsize=(9, 3 * n), squeeze=False)
    for i, (name, row, img, raw, clean, _fit) in enumerate(panels):
        vmax = np.percentile(img, 99.5)
        rs = np.percentile(np.abs(raw), 99)
        div = {"vmin": -rs, "vmax": rs, "cmap": "RdBu_r"}
        for ax, im, title, kw in (
            (
                axes[i, 0],
                img,
                f"{name} ({row['group']}) F814W",
                {"vmin": -0.05 * vmax, "vmax": vmax},
            ),
            (axes[i, 1], raw, f"two-PSF resid S/N={row['resid_snr_raw']:.1f}", div),
            (axes[i, 2], clean, f"halos removed S/N={row['resid_snr']:.1f}", div),
        ):
            ax.imshow(im, origin="lower", **({"cmap": "gray_r"} | kw))
            ax.set_title(title, fontsize=8)
            ax.set_xticks([])
            ax.set_yticks([])
    fig.tight_layout()
    png = path.with_suffix(".png")
    fig.savefig(png, dpi=70)
    Image.open(png).convert("RGB").save(path, quality=80)
    png.unlink()


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(paths.outputs_dir() / "w12_niq_hst"))
    run(ap.parse_args(argv))


if __name__ == "__main__":
    main()
