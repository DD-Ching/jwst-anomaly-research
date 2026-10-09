"""Empirical-PSF version of the D-064 archival HST test (`w12_niq_hst.py` uses a Moffat model).

The PSF is the median of normalised, subpixel-recentred stamps of unsaturated Gaia DR3 stars in
the same HAP F814W skycell cutout (same visit, within ~45" of the pair). The two quasar images are
fitted as shifted, scaled copies of it plus a constant background. The residual flux between the
images (after the halo correction of `w12_niq_hst.py`), its empirical sky-aperture error and the
whole-chain injection test follow that script. The control lens J2308+3201 runs first.
Outputs are ``derived``.

    python scripts/w12_niq_epsf.py --out outputs/w12_niq_epsf
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from scipy import ndimage, optimize

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w12_niq_hst as hst  # noqa: E402

from jwst_anomaly import paths, schema  # noqa: E402

TARGETS = ["J2308+3201", "J0130+0725", "J0728+2607"]  # control first, then the two pairs


@dataclass(frozen=True)
class Params:
    size_px: int = 2200  # hapcut cutout (88" at 0.04"/px): stars within ~44"
    psf_half: int = 25  # PSF stamp half-size in pixels (1")
    star_gmag: tuple = (16.5, 21.0)  # Gaia G range: unsaturated in 674 s ACS/F814W, enough S/N
    star_min_dist_arcsec: float = 4.0  # stars at least this far from the pair
    star_isolation_arcsec: float = 2.0  # no other Gaia source this close to a PSF star
    mirror_ap_arcsec: float = 0.3  # lens-position vs mirror aperture radius


def gaia_stars(c: SkyCoord, radius_deg: float):
    from astroquery.gaia import Gaia

    q = (
        "SELECT source_id, ra, dec, phot_g_mean_mag FROM gaiadr3.gaia_source WHERE "
        f"1=CONTAINS(POINT(ra,dec),CIRCLE({c.ra.deg},{c.dec.deg},{radius_deg}))"
    )
    return Gaia.launch_job(q).get_results()


def stamp(img, x, y, half):
    xi, yi = int(round(x)), int(round(y))
    if (
        xi - half < 0
        or yi - half < 0
        or xi + half + 1 > img.shape[1]
        or yi + half + 1 > img.shape[0]
    ):
        return None, 0.0, 0.0
    return img[yi - half : yi + half + 1, xi - half : xi + half + 1], x - xi, y - yi


def build_epsf(img, valid, xy, half: int):
    """Median of normalised stars, each recentred by its flux-weighted centroid (subpixel shift).
    Returns (psf normalised to unit sum, number of stars used)."""
    stamps = []
    for x, y in xy:
        st, _, _ = stamp(img, x, y, half)
        vs, _, _ = stamp(valid.astype(float), x, y, half)
        if st is None or vs.min() < 1:
            continue
        st = st - np.median(np.r_[st[0], st[-1], st[:, 0], st[:, -1]])  # edge background
        cy, cx = ndimage.center_of_mass(np.clip(st, 0, None) * (st > 0.1 * st.max()))
        st = ndimage.shift(st, (half - cy, half - cx), order=3, mode="nearest")
        if st.sum() > 0:
            stamps.append(st / st.sum())
    if not stamps:
        return None, 0
    psf = np.median(stamps, axis=0)
    psf = np.clip(psf, 0, None)
    return psf / psf.sum(), len(stamps)


def place(psf, shape, x, y):
    """PSF (unit sum) centred at (x, y) on an image of ``shape``."""
    out = np.zeros(shape)
    half = psf.shape[0] // 2
    xi, yi = int(np.floor(x)), int(np.floor(y))
    sh = ndimage.shift(psf, (y - yi, x - xi), order=3, mode="constant")
    y0, x0 = yi - half, xi - half
    ys = slice(max(y0, 0), min(y0 + psf.shape[0], shape[0]))
    xs = slice(max(x0, 0), min(x0 + psf.shape[1], shape[1]))
    out[ys, xs] = sh[ys.start - y0 : ys.stop - y0, xs.start - x0 : xs.stop - x0]
    return out


def fit_two_epsf(img, peaks, psf):
    """Background + two scaled, shifted copies of ``psf``. Returns (fit namespace, residual)."""
    (x1, y1), (x2, y2) = peaks
    bg = float(np.median(img))
    f0 = [float(img[int(y), int(x)] - bg) / psf.max() for x, y in peaks]

    def model(q):
        b, a1, xa, ya, a2, xb, yb = q
        return b + a1 * place(psf, img.shape, xa, ya) + a2 * place(psf, img.shape, xb, yb)

    q0 = [bg, f0[0], x1, y1, f0[1], x2, y2]
    sol = optimize.least_squares(lambda q: (model(q) - img).ravel(), q0, x_scale="jac")
    b, a1, xa, ya, a2, xb, yb = sol.x
    fit = SimpleNamespace(
        x_0_1=SimpleNamespace(value=xa), y_0_1=SimpleNamespace(value=ya),
        x_0_2=SimpleNamespace(value=xb), y_0_2=SimpleNamespace(value=yb),
        amplitude_1=SimpleNamespace(value=a1), amplitude_2=SimpleNamespace(value=a2),
    )  # fmt: skip
    return fit, img - model(sol.x)


def lens_vs_mirror(res, fit, pix, r_ap_arcsec, noise_px):
    """Residual flux in an aperture at the SIS-predicted lens position (on the line between the
    images, at sep * f / (1 + f) from the fainter one) minus the flux in the mirror aperture on
    the far side of the fainter image. A halo centred on that image (host galaxy, PSF colour
    mismatch) gives ~0; a lens galaxy gives a positive excess. Returns (excess, error, S/N)."""
    a1, a2 = fit.amplitude_1.value, fit.amplitude_2.value
    pb = (fit.x_0_1.value, fit.y_0_1.value)  # bright
    pf = (fit.x_0_2.value, fit.y_0_2.value)  # faint
    if a1 < a2:
        pb, pf, a1, a2 = pf, pb, a2, a1
    q = a2 / a1
    lx, ly = pf[0] + (pb[0] - pf[0]) * q / (1 + q), pf[1] + (pb[1] - pf[1]) * q / (1 + q)
    mx, my = 2 * pf[0] - lx, 2 * pf[1] - ly
    yy, xx = np.indices(res.shape)
    r = r_ap_arcsec / pix
    al = np.hypot(xx - lx, yy - ly) <= r
    am = np.hypot(xx - mx, yy - my) <= r
    exc = float(res[al].sum() - res[am].sum())
    err = float(noise_px * np.sqrt(al.sum() + am.sum()))
    return exc, err, exc / err if err > 0 else np.nan


def inject_epsf(img, val, centre, sep_px, tol_px, fit, psf, pix, hp, zp, noise,
                mags=(21, 22, 23, 24, 25)):  # fmt: skip
    """Whole-chain injection with the empirical PSF: an early-type lens galaxy (Sersic n = 4,
    r_eff 0.3") at the SIS-predicted position, then peak finding, the two-ePSF fit, the pair
    check, the halo correction and the residual. Returns [(mag, halo-removed S/N or NaN)]."""
    from astropy.modeling import models

    a1, a2 = fit.amplitude_1.value, fit.amplitude_2.value
    pb = (fit.x_0_1.value, fit.y_0_1.value)
    pf = (fit.x_0_2.value, fit.y_0_2.value)
    if a1 < a2:
        pb, pf, a1, a2 = pf, pb, a2, a1
    q = a2 / a1
    lx, ly = pf[0] + (pb[0] - pf[0]) * q / (1 + q), pf[1] + (pb[1] - pf[1]) * q / (1 + q)
    yy, xx = np.indices(img.shape)
    g = models.Sersic2D(1, 0.3 / pix, 4, x_0=lx, y_0=ly, ellip=0.2, theta=0.5)(xx, yy)
    g /= g.sum()
    out = []
    for mag in mags:
        im2 = img + 10 ** (-0.4 * (mag - zp)) * g
        pk = hst.find_two_peaks(im2, centre, hp.search_arcsec / pix, 0.5 * sep_px, sep_px, tol_px)
        f2, res = fit_two_epsf(im2, pk, psf)
        m = hst.residual_between(hst.remove_halos(res, f2), f2, pix, hp, noise, val)
        good = abs(m["sep_arcsec"] - sep_px * pix) <= hp.sep_match_arcsec
        out.append((mag, m["resid_snr"] if good else np.nan))
    return out


def run(args) -> None:
    p = Params()
    hp = hst.Params()
    out = Path(args.out)
    systems = Table.read(paths.repo_root() / "results" / "w12_niq" / "systems.ecsv")
    rows, inj, panels = [], [], []
    for name in TARGETS:
        r = systems[systems["name"] == name][0]
        c = SkyCoord(float(r["ra"]), float(r["dec"]), unit="deg")
        f = hst.pick_image(list((out / name).glob("*.fits")), c, hp.search_arcsec)
        with fits.open(f) as h:
            hdu = next(x for x in h if x.data is not None and x.data.ndim == 2)
            data = hdu.data.astype(float)
            w = WCS(hdu.header)
            zp = hst.ab_zeropoint(hdu.header, h[0].header)
        valid = np.isfinite(data) & (data != 0)
        img_full = np.where(valid, data, np.median(data[valid]))
        pix = float(np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix))) * 3600)
        g = gaia_stars(c, p.size_px * pix / 3600 / 2)
        gc = SkyCoord(g["ra"], g["dec"], unit="deg")
        _, d2, _ = gc.match_to_catalog_sky(gc, nthneighbor=2)
        ok = np.asarray(g["phot_g_mean_mag"]) > p.star_gmag[0]
        ok &= np.asarray(g["phot_g_mean_mag"]) < p.star_gmag[1]
        ok &= gc.separation(c).arcsec > p.star_min_dist_arcsec
        ok &= d2.arcsec > p.star_isolation_arcsec
        xy = [tuple(map(float, w.world_to_pixel(s))) for s in gc[ok]]
        psf, nstar = build_epsf(img_full, valid, xy, p.psf_half)
        row = {"name": name, "group": str(r["group"]), "n_psf_stars": nstar, "status": "no PSF"}
        if psf is not None:
            x0, y0 = w.world_to_pixel(c)
            sep_px, tol_px = float(r["sep_cat"]) / pix, hp.sep_match_arcsec / pix
            peaks = hst.find_two_peaks(
                img_full, (x0, y0), hp.search_arcsec / pix, 0.5 * sep_px, sep_px, tol_px
            )
            xm, ym = np.mean([q[0] for q in peaks]), np.mean([q[1] for q in peaks])
            half = int(1.5 * r["sep_cat"] / pix) + 10
            y0s, x0s = max(int(ym) - half, 0), max(int(xm) - half, 0)
            sl = (slice(y0s, int(ym) + half), slice(x0s, int(xm) + half))
            img, val = img_full[sl], valid[sl]
            pk = [(x - x0s, y - y0s) for x, y in peaks]
            fit, res = fit_two_epsf(img, pk, psf)
            npix = hst.residual_between(res, fit, pix, hp, 1.0, val)["n_pix"]
            emp = hst.empirical_aperture_noise(data, valid, npix)
            noise = emp / np.sqrt(npix)
            m = hst.residual_between(res, fit, pix, hp, noise, val)
            row.update(m, zp_ab=zp, aperture_noise_empirical=emp)
            # mirror test: same-size apertures (0.3") at the lens position and its mirror image
            n_ap = int(np.pi * (p.mirror_ap_arcsec / pix) ** 2)
            noise_ap = hst.empirical_aperture_noise(data, valid, n_ap) / np.sqrt(n_ap)
            exc, err, snr_m = lens_vs_mirror(res, fit, pix, p.mirror_ap_arcsec, noise_ap)
            row.update(lens_minus_mirror=exc, lens_minus_mirror_err=err, mirror_snr=snr_m)
            # halo test: radially symmetric residual about each image (host galaxy, PSF colour
            # mismatch) removed as in w12_niq_hst; light between the images remains
            mh = hst.residual_between(hst.remove_halos(res, fit), fit, pix, hp, noise, val)
            row.update(halo_removed_snr=mh["resid_snr"], halo_removed_flux=mh["resid_flux"])
            row["status_halo"] = (
                "light between images" if mh["resid_snr"] >= hp.detect_sigma else "none between"
            )
            if row["status_halo"] == "none between":
                rec = inject_epsf(
                    img, val, (x0 - x0s, y0 - y0s), sep_px, tol_px, fit, psf, pix, hp, zp, noise
                )
                base = mh["resid_snr"]
                for mag, s_inj in rec:
                    inj.append(
                        {"name": name, "inject_mag": mag, "snr_above_baseline": s_inj - base}
                    )
                row["mag_limit_inject"] = hst.contiguous_limit(
                    [(mag, bool(s_inj - base >= hp.detect_sigma)) for mag, s_inj in rec]
                )
            snr = m["resid_snr"]
            row["status"] = (
                "pair mismatch"
                if abs(m["sep_arcsec"] - float(r["sep_cat"])) > hp.sep_match_arcsec
                else "residual >= 5 sigma"
                if snr >= hp.detect_sigma
                else "no residual"
            )
            row["mag_resid"] = (
                zp - 2.5 * np.log10(m["resid_flux"]) if m["resid_flux"] > 0 else np.nan
            )
            panels.append((name, row, img, res))
        rows.append(row)
        keys = ("name", "status_halo", "n_psf_stars", "halo_removed_snr", "mag_limit_inject")
        print({k: row.get(k) for k in keys})
    t = Table(rows=rows)
    t.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source="scripts/w12_niq_epsf.py: empirical Gaia-star PSF on MAST HAP F814W cutouts",
        params=asdict(p),
    )
    res_dir = paths.repo_root() / "results" / "w12_niq"
    t.write(res_dir / "hst_epsf_residuals.ecsv", overwrite=True)
    ti = Table(rows=inj) if inj else Table(names=("name", "inject_mag"))
    ti.meta.update(
        provenance=schema.Provenance.SIMULATED.value,
        source="Sersic n=4, r_eff 0.3 arcsec at the SIS-predicted position; whole ePSF chain",
    )
    ti.write(res_dir / "hst_epsf_injections.ecsv", overwrite=True)
    sheet(panels, res_dir / "hst_epsf_residuals.jpg")


def sheet(panels, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image

    fig, axes = plt.subplots(len(panels), 3, figsize=(9, 3 * len(panels)), squeeze=False)
    for i, (name, row, img, res) in enumerate(panels):
        vmax = np.percentile(img, 99.5)
        rs = np.percentile(np.abs(res), 99)
        sm = ndimage.gaussian_filter(res, 2.0)
        ss = np.percentile(np.abs(sm), 99.5)
        gray = {"cmap": "gray_r", "vmin": -0.05 * vmax, "vmax": vmax}

        def div(v):
            return {"cmap": "RdBu_r", "vmin": -v, "vmax": v}

        for ax, im, title, kw in (
            (axes[i, 0], img, f"{name} F814W", gray),
            (axes[i, 1], res, f"ePSF residual S/N={row['resid_snr']:.1f}", div(rs)),
            (
                axes[i, 2],
                sm,
                f"residual smoothed 2 px (halo-removed S/N={row['halo_removed_snr']:.1f})",
                div(ss),
            ),
        ):
            ax.imshow(im, origin="lower", **kw)
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
    ap.add_argument("--out", default=str(paths.outputs_dir() / "w12_niq_epsf"))
    run(ap.parse_args(argv))


if __name__ == "__main__":
    main()
