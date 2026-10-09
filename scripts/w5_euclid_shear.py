"""W5/W1 radial-shear screen on Euclid Q1 MER shapes (D-065 next; D-TBD).

A negative-mass lens (W5 / W1, D-047) shears background galaxies *radially*; ordinary mass shears
them tangentially. This screen measures the catalogue aperture-mass S/N of radial alignment
(``exotic_screens.ApertureMass``, D-050) on a grid of trial centres, from the Euclid Q1 MER
catalogue's SExtractor moments (IRSA TAP; no PSF-corrected shear catalogue is public in Q1).

Steps (each writes into ``results/w5_shear/``):

``pacheck``  Moments measured on MER VIS image cutouts against the catalogue ``position_angle``
    (lensing-independent check of the angle convention).
``validate``  Known massive clusters (Planck PSZ2 / ACT DR5, ``CLUSTERS``) must show *tangential*
    shear (S < 0). This checks the angle convention (``pa_east_of_north``): a 90° error flips
    tangential and radial, the sign under test. Also measures the stars' mean ellipticity (PSF
    anisotropy; a constant one cancels over a full annulus).
``screen``  Trial-centre grid over pilot discs; the null rotates every shape (positions kept); the
    field maximum of S over the grid is compared with the rotation null. Injections add the
    reduced shear of a negative point mass, g = +R (θ_E/θ)² exp(2iφ) (radial), to the real shapes at
    random centres and run the same statistic.

Everything fetched is cached under ``$JWST_ANOMALY_DATA/euclid_q1_shear/`` (gitignored); outputs are
``derived``; θ_E grids, cuts and the responsivity R are ASSUMPTIONs named below.

    python scripts/w5_euclid_shear.py pacheck
    python scripts/w5_euclid_shear.py validate
    python scripts/w5_euclid_shear.py screen
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import requests
from astropy.table import Table
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exotic_screens as es  # noqa: E402

from jwst_anomaly import paths  # noqa: E402

TAP = "https://irsa.ipac.caltech.edu/TAP/sync"
TABLE = "euclid_q1_mer_catalogue"
# "position_angle" must be quoted on IRSA (unquoted, the ADQL parser rejects the query).
COLUMNS = (
    'ra, dec, ellipticity, "position_angle", semimajor_axis, flux_detection_total, point_like_prob'
)
GAL_SEL = "vis_det = 1 AND spurious_flag = 0 AND det_quality_flag = 0 AND point_like_prob < 0.1"
STAR_SEL = "vis_det = 1 AND spurious_flag = 0 AND det_quality_flag = 0 AND point_like_prob > 0.9"
VIS_LIMIT = 24.5  # galaxies brighter than this in VIS (ASSUMPTION, as D-065)
STAR_MAGS = (19.0, 22.0)  # unsaturated, high-S/N stars for the PSF (ASSUMPTION)
MIN_SIZE_OVER_PSF = 1.2  # resolved: semimajor sigma > this x the stars' median sigma (ASSUMPTION)
# Isophotal-moment shear responsivity for injections (ASSUMPTION; D-053 measured 0.4-0.6).
RESPONSIVITY = 0.5
# Known clusters inside Q1 coverage (VizieR J/A+A/594/A27 PSZ2 and J/ApJS/253/3 ACT DR5, queried
# 2026-10-09): name -> (RA, Dec, z, M500 / 1e14 Msun). Observed catalogue values.
CLUSTERS = {
    "PSZ2 G255.60-46.18": (62.7892169, -48.3030180, 0.4235, 6.32),
    "ACT-CL J0405.9-4915": (61.4923159, -49.2626922, 0.3250, 3.61),
    "ACT-CL J0405.1-4648": (61.2836742, -46.8110099, 0.3789, 2.91),
    "ACT-CL J0402.2-4611": (60.5738725, -46.1881411, 0.3606, 2.67),
}
CLUSTER_RADIUS_DEG = 0.2
CLUSTER_APERTURE = (60.0, 600.0)  # arcsec, validation annulus radii (ASSUMPTION)
# Pilot discs (the D-065 field centres), away from the validation clusters (ASSUMPTION).
PILOTS = {"EDF-F": (52.93, -28.09), "EDF-S": (61.24, -48.42), "EDF-N": (269.73, 66.02)}
# The EDF-N 0.3° row query did not return in 15 min (2026-10-09).
DEFAULT_PILOTS = ("EDF-F", "EDF-S")
PILOT_RADIUS_DEG = 0.3
THETA_E_ARCSEC = (30.0, 60.0, 120.0)  # trial Einstein radii; aperture 1.5-3 θ_E (D-065 forecast)
GRID_STEP_OVER_THETA_E = 1.0  # grid spacing in θ_E (ASSUMPTION)
N_NULL = 200  # rotation draws (ASSUMPTION)
N_INJ = 40  # injections per θ_E and field (ASSUMPTION)
MIN_N = 20  # sources needed in an aperture (ASSUMPTION)
MIN_EFFICIENCY = 0.5  # quote a limit only if every field's efficiency reaches this (ASSUMPTION)


def vis_flux_ujy(mag: float) -> float:
    return 10 ** ((23.9 - mag) / 2.5)


def disc_query(ra: float, dec: float, radius: float, star: bool) -> str:
    sel = STAR_SEL if star else GAL_SEL
    lo = STAR_MAGS[1] if star else VIS_LIMIT
    flux = f"flux_detection_total > {vis_flux_ujy(lo):.4f}"
    if star:
        flux += f" AND flux_detection_total < {vis_flux_ujy(STAR_MAGS[0]):.4f}"
    return (
        f"SELECT {COLUMNS} FROM {TABLE} WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), "
        f"CIRCLE('ICRS', {ra}, {dec}, {radius})) AND {sel} AND {flux}"
    )


def fetch(ra: float, dec: float, radius: float, star: bool, cache: Path) -> Table:
    """Rows in a disc (CONTAINS uses IRSA's spatial index, D-065), cached as ECSV."""
    q = disc_query(ra, dec, radius, star)
    key = hashlib.sha1(q.encode()).hexdigest()[:10]  # a changed selection never reuses old rows
    path = cache / f"{'star' if star else 'gal'}_{ra:.4f}_{dec:+.4f}_{radius:.3f}_{key}.ecsv"
    if path.exists():
        return Table.read(path)
    for attempt in range(4):
        try:
            r = requests.post(TAP, data={"QUERY": q, "FORMAT": "csv"}, timeout=900)
            r.raise_for_status()
            # IRSA reports query errors as a VOTable with HTTP 200.
            if r.text.lstrip().startswith("<"):
                raise RuntimeError(f"IRSA TAP error: {r.text[:400]}")
            t = Table.read(r.text, format="ascii.csv")
            break
        except (requests.RequestException, ValueError, RuntimeError):
            if attempt == 3:
                raise
            time.sleep(2 ** (attempt + 1))
    cache.mkdir(parents=True, exist_ok=True)
    t.write(path, overwrite=True)
    return t


def tangent_plane(ra, dec, ra0: float, dec0: float) -> tuple[np.ndarray, np.ndarray]:
    """Gnomonic offsets in arcsec, x West and y North (the ``ApertureMass`` frame)."""
    a, d = np.deg2rad(np.asarray(ra, float)), np.deg2rad(np.asarray(dec, float))
    a0, d0 = math.radians(ra0), math.radians(dec0)
    cosc = np.sin(d0) * np.sin(d) + np.cos(d0) * np.cos(d) * np.cos(a - a0)
    xi = np.cos(d) * np.sin(a - a0) / cosc  # East
    eta = (np.cos(d0) * np.sin(d) - np.sin(d0) * np.cos(d) * np.cos(a - a0)) / cosc
    k = 180.0 / math.pi * 3600.0
    return -xi * k, eta * k


def pa_east_of_north(position_angle_deg) -> np.ndarray:
    """Catalogue ``position_angle`` -> PA east of north (mod 180°).

    The TAP column description says "CCW/x (THETA_IMAGE)", but moments measured on 25 MER VIS
    cutouts (north-up, East left) give PA east of north = ``position_angle`` to 0.8° (median), and
    ``position_angle`` + 90° is 89° off (2026-10-09, ``validate`` notes). Known clusters then show
    tangential shear (S < 0), as positive mass must."""
    return np.mod(np.asarray(position_angle_deg, float), 180.0)


def star_psf(stars: Table) -> dict:
    """Median sigma (px) and mean complex ellipticity of stars (PSF size and anisotropy)."""
    a = np.asarray(stars["semimajor_axis"], float)
    ell = np.asarray(stars["ellipticity"], float)
    eps = es.intrinsic_ellipticity(
        a, a * (1.0 - ell), pa_east_of_north(stars["position_angle"]), 0.0
    )
    ok = np.isfinite(eps) & np.isfinite(a)
    e = eps[ok]
    return {
        "n_stars": int(ok.sum()),
        "sigma_px_median": float(np.median(a[ok] * np.sqrt(1.0 - ell[ok]))),
        "mean_e1": float(e.real.mean()),
        "mean_e2": float(e.imag.mean()),
        # error of the mean of one component (e1 or e2): rms|e| / sqrt(2 N)
        "mean_e_err": float(np.sqrt(np.mean(np.abs(e) ** 2) / (2 * len(e)))),
        "rms_abs_e": float(np.sqrt(np.mean(np.abs(e) ** 2))),
    }


def galaxy_shapes(gals: Table, psf_sigma: float) -> np.ndarray:
    """PSF-deconvolved complex ellipticity (PA east of north), NaN for unresolved rows."""
    a = np.asarray(gals["semimajor_axis"], float)
    ell = np.asarray(gals["ellipticity"], float)
    eps = es.intrinsic_ellipticity(
        a, a * (1.0 - ell), pa_east_of_north(gals["position_angle"]), psf_sigma
    )
    resolved = a * np.sqrt(np.clip(1.0 - ell, 0, 1)) > MIN_SIZE_OVER_PSF * psf_sigma
    return np.where(resolved & (np.abs(eps) < 1.0), eps, np.nan + 0j)


def sheared_catalogue(gals: Table, g: np.ndarray, rng: np.random.Generator) -> Table:
    """Catalogue copy whose *observed* moments carry the reduced shear ``g`` (for injections).

    Each row's observed ellipticity (no PSF deconvolution) is rotated by a random angle (removing
    any real signal, as the null does) and sheared; ``ellipticity``, ``position_angle`` and
    ``semimajor_axis`` (area kept) are rewritten, so the resolved cut and the PSF deconvolution of
    :func:`galaxy_shapes` run on the injected shapes (whole chain)."""
    a = np.asarray(gals["semimajor_axis"], float)
    ell = np.clip(np.asarray(gals["ellipticity"], float), 0.0, 0.999)
    pa = pa_east_of_north(gals["position_angle"])
    q = 1.0 - ell
    e_obs = (1.0 - q) / (1.0 + q) * np.exp(2j * np.deg2rad(pa))
    e_new = apply_shear(e_obs * np.exp(1j * rng.uniform(0, 2 * np.pi, len(a))), g)
    m = np.minimum(np.abs(e_new), 0.999)
    q_new = (1.0 - m) / (1.0 + m)
    return Table(
        {
            "semimajor_axis": a * np.sqrt(q / q_new),
            "ellipticity": 1.0 - q_new,
            "position_angle": np.rad2deg(np.angle(e_new) / 2.0),
        }
    )


def radial_shear(x, y, cx: float, cy: float, theta_e: float, r_min: float) -> np.ndarray:
    """Complex reduced shear of a negative point mass of Einstein radius θ_E at (cx, cy), radial:
    g = (θ_E/θ)² exp(2iφ), φ the PA of the source seen from the centre (model_prediction).
    Zero inside ``r_min`` (κ = 0 outside θ_E; inner sources are excluded by the aperture anyway)."""
    de, dn = -(np.asarray(x) - cx), np.asarray(y) - cy
    th = np.hypot(de, dn)
    phi = np.arctan2(de, dn)
    with np.errstate(divide="ignore"):
        amp = np.where(th >= r_min, (theta_e / np.maximum(th, 1e-9)) ** 2, 0.0)
    return np.minimum(amp, 0.9) * np.exp(2j * phi)


def apply_shear(e: np.ndarray, g: np.ndarray) -> np.ndarray:
    """Lensed ellipticity (ε + g)/(1 + g* ε) (Seitz & Schneider 1997), |g| < 1."""
    return (e + g) / (1.0 + np.conj(g) * e)


def cmd_validate(args) -> dict:
    cache = paths.data_root() / "euclid_q1_shear"
    jobs = [(n, c) for n, c in CLUSTERS.items()]
    with ThreadPoolExecutor(8) as ex:  # galaxy and star queries all in flight at once
        fg = [ex.submit(fetch, c[0], c[1], CLUSTER_RADIUS_DEG, False, cache) for _, c in jobs]
        fs = [ex.submit(fetch, c[0], c[1], CLUSTER_RADIUS_DEG, True, cache) for _, c in jobs]
        gal, star = [f.result() for f in fg], [f.result() for f in fs]
    out = {}
    for (name, (ra, dec, z, m)), g, s in zip(jobs, gal, star, strict=True):
        psf = star_psf(s)
        e = galaxy_shapes(g, psf["sigma_px_median"])
        x, y = tangent_plane(g["ra"], g["dec"], ra, dec)
        ok = np.isfinite(e)
        res = {"z": z, "m500_1e14": m, "n_gal": len(g), "n_resolved": int(ok.sum()), "psf": psf}
        for kind in ("tophat", "pointmass"):
            r_in, r_out = CLUSTER_APERTURE
            am = es.ApertureMass(x[ok], y[ok], [0.0], [0.0], r_out, kind, r_in, MIN_N)
            sv, sx = am.snr(e[ok])
            res[kind] = {"S": float(sv[0]), "S_cross": float(sx[0]), "n": int(am.n[0])}
            # Same statistic with PA rotated by 90° (the TAP description's reading) for the record.
            sw, _ = am.snr(-e[ok])
            res[kind]["S_if_pa_plus_90"] = float(sw[0])
        # Faint galaxies only (VIS > 23): fewer cluster members, same sign expected.
        mag = 23.9 - 2.5 * np.log10(np.asarray(g["flux_detection_total"], float))
        fk = ok & (mag > 23.0)
        r_in, r_out = CLUSTER_APERTURE
        am = es.ApertureMass(x[fk], y[fk], [0.0], [0.0], r_out, "tophat", r_in, MIN_N)
        res["tophat_faint_vis_gt_23"] = {"S": float(am.snr(e[fk])[0][0]), "n": int(am.n[0])}
        # Tangential-shear profile in annuli (derived), for the record.
        prof = []
        for r_in, r_out in ((60, 150), (150, 300), (300, 600)):
            am = es.ApertureMass(x[ok], y[ok], [0.0], [0.0], r_out, "tophat", r_in, MIN_N)
            proj = am.w @ e[ok]
            n = int(am.n[0])
            prof.append(
                {
                    "r_arcsec": [r_in, r_out],
                    "n": n,
                    "mean_e_t": float(-proj[0].real / max(n, 1)),
                    "mean_e_x": float(-proj[0].imag / max(n, 1)),
                    "err": float(np.sqrt(np.mean(np.abs(e[ok]) ** 2) / 2 / max(n, 1))),
                }
            )
        res["profile"] = prof
        out[name] = res
        print(name, json.dumps({k: res[k] for k in ("n_resolved", "tophat", "pointmass")}))
    return out


def screen_field(name: str, cache: Path, rng: np.random.Generator) -> dict:
    ra0, dec0 = PILOTS[name]
    g = fetch(ra0, dec0, PILOT_RADIUS_DEG, False, cache)
    s = fetch(ra0, dec0, PILOT_RADIUS_DEG, True, cache)
    psf = star_psf(s)
    e_all = galaxy_shapes(g, psf["sigma_px_median"])
    x_all, y_all = tangent_plane(g["ra"], g["dec"], ra0, dec0)
    ok = np.isfinite(e_all)
    x, y, e = x_all[ok], y_all[ok], e_all[ok]
    rmax = PILOT_RADIUS_DEG * 3600.0
    res = {"n_gal": len(g), "n_resolved": int(ok.sum()), "psf": psf, "theta_e": {}}
    for te in THETA_E_ARCSEC:
        r_in, r_out = 1.5 * te, 3.0 * te
        step = GRID_STEP_OVER_THETA_E * te
        gx, gy = np.meshgrid(np.arange(-rmax, rmax + 1, step), np.arange(-rmax, rmax + 1, step))
        inside = np.hypot(gx, gy) <= rmax - r_out  # full aperture inside the disc
        cx, cy = gx[inside], gy[inside]
        am = es.ApertureMass(x, y, cx, cy, r_out, "pointmass", r_in, MIN_N)
        sv, sx = am.snr(e)
        null = am.null_max(e, N_NULL, rng)
        smax = float(np.nanmax(sv))
        thr = float(np.quantile(null, 0.99))  # field-wise 1 % false-alarm threshold (ASSUMPTION)
        # Injections: one negative point mass at a uniformly random position inside the centre
        # region (off-grid); detected if any centre within one grid step has S > thr.
        tree = cKDTree(np.c_[cx, cy])
        det = []
        for _ in range(N_INJ):
            rr = (rmax - r_out) * math.sqrt(rng.uniform())
            aa = rng.uniform(0, 2 * np.pi)
            ix, iy = rr * math.cos(aa), rr * math.sin(aa)
            gi = RESPONSIVITY * radial_shear(x_all, y_all, ix, iy, te, r_in)
            ei = galaxy_shapes(sheared_catalogue(g, gi, rng), psf["sigma_px_median"])
            oi = np.isfinite(ei)
            near = tree.query_ball_point([ix, iy], step)
            am_i = es.ApertureMass(
                x_all[oi], y_all[oi], cx[near], cy[near], r_out, "pointmass", r_in, MIN_N
            )
            si, _ = am_i.snr(ei[oi])
            det.append(bool(np.any(si > thr)))  # NaN (too few sources) counts as missed
        res["theta_e"][f"{te:g}"] = {
            "n_centres": int(len(cx)),
            "area_deg2": float(len(cx) * step**2 / 3600.0**2),
            "S_max": smax,
            "S_min": float(np.nanmin(sv)),
            "S_cross_absmax": float(np.nanmax(np.abs(sx))),
            "null_max_q50": float(np.median(null)),
            "threshold_p01": thr,
            "p_random_of_S_max": float(np.mean(null >= smax)),
            "flags": [
                {"ra_off_arcsec": float(cx[i]), "dec_off_arcsec": float(cy[i]), "S": float(sv[i])}
                for i in np.flatnonzero(sv > thr)
            ],
            "injection_efficiency": float(np.mean(det)),
        }
        print(
            name,
            te,
            json.dumps({k: v for k, v in res["theta_e"][f"{te:g}"].items() if k != "flags"}),
        )
    return res


def cmd_screen(args) -> dict:
    cache = paths.data_root() / "euclid_q1_shear"
    names = list(DEFAULT_PILOTS) if not args.field else [args.field]
    with ThreadPoolExecutor(8) as ex:  # prefetch every disc, galaxies and stars, concurrently
        futs = [
            ex.submit(fetch, *PILOTS[n], PILOT_RADIUS_DEG, st, cache)
            for n in names
            for st in (False, True)
        ]
        [f.result() for f in futs]
    # one generator per field (seed, field index) so a single-field run reproduces its numbers
    out = {
        n: screen_field(n, cache, np.random.default_rng([args.seed, list(PILOTS).index(n)]))
        for n in names
    }
    out["limits"] = density_limits(out)
    return out


def density_limits(fields: dict) -> dict:
    """95 % upper limit on the sky density of negative point masses per θ_E, with zero flags:
    n95 = 3 / Σ ε_i A_i over fields (derived; Poisson, efficiency from the injections)."""
    out = {}
    for te in THETA_E_ARCSEC:
        k = f"{te:g}"
        rows = [f["theta_e"][k] for f in fields.values()]
        n_flags = sum(len(r["flags"]) for r in rows)
        eff_area = sum(r["injection_efficiency"] * r["area_deg2"] for r in rows)
        out[k] = {
            "n_flags": n_flags,
            "effective_area_deg2": eff_area,
            # no limit where injections show the test is blind (scripts/CLAUDE.md)
            "n95_deg2": 3.0 / eff_area
            if n_flags == 0 and min(r["injection_efficiency"] for r in rows) >= MIN_EFFICIENCY
            else None,
        }
    return out


def moment_pa(img: np.ndarray, wcs, ra: float, dec: float, radius_px: float) -> float:
    """PA east of north (deg, mod 180) of the flux-weighted second moments of ``img`` within
    ``radius_px`` of (ra, dec), pixels above 3x the edge-strip rms; the axis is mapped to the sky
    through the WCS, so no orientation of the image is assumed."""
    yy, xx = np.indices(img.shape)
    px, py = wcs.world_to_pixel_values(ra, dec)
    edge = np.r_[img[:4].ravel(), img[-4:].ravel()]
    m = (np.hypot(xx - px, yy - py) < radius_px) & (img > 3.0 * np.nanstd(edge))
    f = np.where(m, img, 0.0)
    dx, dy = xx - px, yy - py
    qxx, qyy, qxy = ((f * dx * dx).sum(), (f * dy * dy).sum(), (f * dx * dy).sum())
    t = 0.5 * math.atan2(2.0 * qxy, qxx - qyy)  # axis angle in the pixel frame
    r2, d2 = wcs.pixel_to_world_values(px + 5.0 * math.cos(t), py + 5.0 * math.sin(t))
    return float(np.degrees(np.arctan2((r2 - ra) * math.cos(math.radians(dec)), d2 - dec)) % 180.0)


def axis_diff(a, b) -> np.ndarray:
    """Signed difference of two axis angles (deg), wrapped to [-90, 90)."""
    d = np.asarray(a, float) - np.asarray(b, float)
    return np.mod(d + 90.0, 180.0) - 90.0


# One MER VIS tile with ACT-CL J0405.1-4648 (IRSA SIA, collection euclid_DpdMerBksMosaic).
PACHECK_TILE = (
    "https://irsa.ipac.caltech.edu/ibe/data/euclid/q1/MER/102022477/VIS/"
    "EUC_MER_BGSUB-MOSAIC-VIS_TILE102022477-C38626_20241018T201336.984902Z_00.00.fits"
)
PACHECK_CENTRE = (61.2837, -46.8110, 0.12)  # RA, Dec, radius (deg) inside that tile


def cmd_pacheck(args) -> dict:
    """Catalogue ``position_angle`` against moments measured on 8″ VIS cutouts (IRSA IBE cutouts)
    of bright, elongated galaxies (VIS < 21.5, ellipticity > 0.4, semimajor > 4 px)."""
    import gzip
    import io

    from astropy.io import fits
    from astropy.wcs import WCS

    ra0, dec0, rad = PACHECK_CENTRE
    ra_c, dec_c = CLUSTERS["ACT-CL J0405.1-4648"][:2]
    g = fetch(ra_c, dec_c, CLUSTER_RADIUS_DEG, False, paths.data_root() / "euclid_q1_shear")
    mag = 23.9 - 2.5 * np.log10(np.asarray(g["flux_detection_total"], float))
    x, y = tangent_plane(g["ra"], g["dec"], ra0, dec0)
    sel = np.flatnonzero(
        (mag < 21.5)
        & (np.asarray(g["ellipticity"], float) > 0.4)
        & (np.asarray(g["semimajor_axis"], float) > 4.0)
        & (np.hypot(x, y) < rad * 3600.0)
    )[: args.n]

    def one(i):
        ra, dec = float(g["ra"][i]), float(g["dec"][i])
        r = requests.get(
            PACHECK_TILE, params={"center": f"{ra},{dec}", "size": "8arcsec"}, timeout=120
        )
        r.raise_for_status()
        b = gzip.decompress(r.content) if r.content[:2] == b"\x1f\x8b" else r.content
        with fits.open(io.BytesIO(b)) as h:
            img, wcs = h[0].data.astype(float), WCS(h[0].header)
        return moment_pa(img, wcs, ra, dec, 3.0 * float(g["semimajor_axis"][i]))

    with ThreadPoolExecutor(8) as ex:
        pa_img = np.array(list(ex.map(one, sel)))
    cat = np.asarray(g["position_angle"], float)[sel]
    res = {"n": len(sel), "pa_image_deg": pa_img.round(2).tolist(), "position_angle": cat.tolist()}
    for lab, v in (("as_is", cat), ("plus_90", cat + 90.0), ("mirror", -cat)):
        res[f"median_abs_diff_{lab}"] = float(np.median(np.abs(axis_diff(v, pa_img))))
    print(json.dumps({k: v for k, v in res.items() if k.startswith("median") or k == "n"}))
    return res


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("validate")
    pc = sub.add_parser("pacheck")
    pc.add_argument("-n", type=int, default=25)
    s = sub.add_parser("screen")
    s.add_argument("--field", choices=list(PILOTS))
    s.add_argument("--seed", type=int, default=20261009)
    args = p.parse_args(argv)
    out_dir = paths.repo_root() / "results" / "w5_shear"
    out_dir.mkdir(parents=True, exist_ok=True)
    res = {"validate": cmd_validate, "screen": cmd_screen, "pacheck": cmd_pacheck}[args.cmd](args)
    meta = {
        "provenance": {
            "shapes": "observed (Euclid Q1 MER catalogue, IRSA TAP; SExtractor moments)",
            "statistics": "derived",
            "injections": "model_prediction (negative point mass, radial reduced shear)",
            "service": f"{TAP} table {TABLE}",
            # run time; IRSA rows may come from an earlier cache (file mtimes under the data root)
            "run_utc": time.strftime("%Y-%m-%dT%H:%MZ", time.gmtime()),
            "first_irsa_access_utc": "2026-10-09",
            "assumptions": {
                "vis_limit": VIS_LIMIT,
                "star_mags": STAR_MAGS,
                "min_size_over_psf": MIN_SIZE_OVER_PSF,
                "responsivity": RESPONSIVITY,
                "galaxy_selection": GAL_SEL,
                "star_selection": STAR_SEL,
                "theta_e_arcsec": THETA_E_ARCSEC,
                "grid_step_over_theta_e": GRID_STEP_OVER_THETA_E,
                "aperture_over_theta_e": [1.5, 3.0],
                "threshold": "99th percentile of the field maximum under shape rotations",
                "n_null": N_NULL,
                "n_inj": N_INJ,
                "min_n": MIN_N,
                "min_efficiency": MIN_EFFICIENCY,
                "cluster_aperture_arcsec": CLUSTER_APERTURE,
                "pilot_radius_deg": PILOT_RADIUS_DEG,
                "seed": getattr(args, "seed", None),
                "pa_mapping": "PA E of N = position_angle (25 MER VIS cutouts, 4 clusters)",
            },
        }
    }
    field = getattr(args, "field", None)
    path = out_dir / (f"{args.cmd}_{field}.json" if field else f"{args.cmd}.json")
    path.write_text(json.dumps({**meta, "results": res}, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
