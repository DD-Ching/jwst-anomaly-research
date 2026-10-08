"""Galaxy count maps for the W5 count-deficit screen (D-063).

A repulsive (negative-mass) point lens demagnifies everything behind it inside about one Einstein
radius, so bright background galaxies thin out there (``exotic_sim.count_ratio``). Stars in the
Milky Way lie in front of any extragalactic lens and are unaffected, while masks and survey holes
remove both. This module bins a Legacy Surveys DR10 Tractor extract onto a flat tangent-plane
grid, measures the galaxy count in a disk at every cell against the coverage-weighted mean, and
injects predicted holes for efficiency. Exotic physics is a hypothesis; a flag is an anomaly to
vet, and a null becomes an injection-calibrated limit.

Data access: the Astro Data Lab TAP service rejects expressions in GROUP BY and sub-selects in
FROM (tested 2026-10-08), so counts are binned client-side from row queries.
"""

from __future__ import annotations

import numpy as np
from astropy.table import Table
from scipy import signal, stats

from . import exotic_sim, schema

DATALAB_TAP = "https://datalab.noirlab.edu/tap"
TRACTOR_COLUMNS = ("ra", "dec", "type", "mag_r", "maskbits")
# DR10 maskbits used as the mask: NPRIMARY, BRIGHT, SATUR_G/R/Z, ALLMASK_G/R/Z, MEDIUM, GALAXY,
# CLUSTER. WISE bits (8, 9) and the i-band bits (14, 15) are ignored: the screen is optical r-band.
MASK_BITS = sum(1 << b for b in (0, 1, 2, 3, 4, 5, 6, 7, 11, 12, 13))


def fetch_tractor(ra_range, dec_range, mag_max: float, tap_url: str = DATALAB_TAP) -> Table:
    """All DR10 Tractor rows with ``mag_r < mag_max`` in a RA/Dec box (one TAP query, observed)."""
    import pyvo

    (ra0, ra1), (de0, de1) = ra_range, dec_range
    query = (
        f"SELECT {', '.join(TRACTOR_COLUMNS)} FROM ls_dr10.tractor "
        f"WHERE ra BETWEEN {ra0} AND {ra1} AND dec BETWEEN {de0} AND {de1} AND mag_r < {mag_max}"
    )
    tab = pyvo.dal.TAPService(tap_url).run_sync(query, maxrec=5_000_000).to_table()
    tab.meta.update(
        provenance=schema.Provenance.OBSERVED.value,
        source=f"Legacy Surveys DR10 ls_dr10.tractor via {tap_url}: {query}",
    )
    return tab


def tangent_plane(ra, dec, ra0: float, dec0: float) -> tuple[np.ndarray, np.ndarray]:
    """Gnomonic projection in arcmin about (ra0, dec0); x grows with RA."""
    ra, dec = np.radians(ra), np.radians(dec)
    a0, d0 = np.radians(ra0), np.radians(dec0)
    cosc = np.sin(d0) * np.sin(dec) + np.cos(d0) * np.cos(dec) * np.cos(ra - a0)
    x = np.cos(dec) * np.sin(ra - a0) / cosc
    y = (np.cos(d0) * np.sin(dec) - np.sin(d0) * np.cos(dec) * np.cos(ra - a0)) / cosc
    return np.degrees(x) * 60.0, np.degrees(y) * 60.0


def count_grids(
    tab: Table,
    ra0: float,
    dec0: float,
    half_size_arcmin: float,
    cell_arcmin: float,
    mag_gal: float,
    mask_bits: int = MASK_BITS,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per-cell counts (all ``derived``): galaxies, stars, unmasked rows and all rows.

    Galaxies: ``type != PSF``, ``maskbits & mask_bits == 0``, ``mag_r < mag_gal``. Stars: unmasked
    ``PSF`` rows. The unmasked fraction of *all* rows traces masks independently of a lens, which
    thins galaxies but not the fraction; with about one row per cell it is only used summed over a
    disk (``disk_deficit``).
    """
    clean = (np.asarray(tab["maskbits"]) & mask_bits) == 0
    psf = np.char.strip(np.asarray(tab["type"]).astype(str)) == "PSF"
    gal = clean & ~psf & (np.asarray(tab["mag_r"]) < mag_gal)
    sels = (gal, clean & psf, clean, np.ones(len(tab), bool))
    return tuple(grid(tab, sel, ra0, dec0, half_size_arcmin, cell_arcmin) for sel in sels)


def grid(tab: Table, sel, ra0: float, dec0: float, half_size_arcmin: float, cell_arcmin: float):
    """Counts of the rows ``sel`` per cell of the square tangent-plane grid (rows = +Dec)."""
    x, y = tangent_plane(np.asarray(tab["ra"])[sel], np.asarray(tab["dec"])[sel], ra0, dec0)
    n = int(round(2 * half_size_arcmin / cell_arcmin))
    edges = np.linspace(-half_size_arcmin, half_size_arcmin, n + 1)
    return np.histogram2d(y, x, bins=(edges, edges))[0]


def disk_kernel(radius_cells: float) -> np.ndarray:
    r = int(np.ceil(radius_cells))
    yy, xx = np.mgrid[-r : r + 1, -r : r + 1]
    return (np.hypot(xx, yy) <= radius_cells).astype(float)


def disk_deficit(
    gal: np.ndarray, clean: np.ndarray, every: np.ndarray, radius_cells: float, min_coverage: float
):
    """Observed and expected galaxy counts in a disk at every cell, and the Poisson tail.

    Unmasked fraction f = sum(clean) / sum(every) over the disk; E = n_bar * A * f with A the disk
    area inside the grid and n_bar the field's galaxies per cell per unit f. Disks with f <
    ``min_coverage`` (ASSUMPTION) or no rows get NaN. Returns (O, E, log10 p) with
    p = P(N <= O | Poisson(E)); clustering widens the true null, so p ranks and the empirical
    distribution over the field (plus injections) calibrates.
    """
    k = disk_kernel(radius_cells)

    def conv(a):
        return signal.fftconvolve(a, k, mode="same")

    obs = np.rint(conv(gal))
    n_every = conv(every)
    frac = np.divide(conv(clean), n_every, out=np.zeros_like(n_every), where=n_every > 0.5)
    n_bar = gal.sum() / (gal.size * clean.sum() / every.sum())
    exp = n_bar * conv(np.ones_like(gal)) * frac
    ok = frac >= min_coverage
    logp = np.where(ok, stats.poisson.logcdf(obs, np.where(ok, exp, 1.0)) / np.log(10), np.nan)
    return np.where(ok, obs, np.nan), np.where(ok, exp, np.nan), logp


def clustering_k(obs: np.ndarray, exp: np.ndarray, radius_cells: float) -> tuple[float, int]:
    """Negative-binomial shape k of disk counts, from disks on a grid spaced one diameter apart.

    Var(O) = E + E²/k (method of moments on independent disks); galaxy clustering makes k finite.
    Returns (k, number of independent disks). ``derived``; k = inf if no excess variance.
    """
    s = max(int(np.ceil(2 * radius_cells)), 1)
    o, e = obs[::s, ::s].ravel(), exp[::s, ::s].ravel()
    ok = np.isfinite(o) & np.isfinite(e) & (e > 0)
    o, e = o[ok], e[ok]
    excess = ((o - e) ** 2 - e).sum() / (e**2).sum()
    return (1.0 / excess if excess > 0 else np.inf), int(ok.sum())


def nb_log10_tail(obs, exp, k: float) -> np.ndarray:
    """log10 P(N <= obs), negative binomial with mean ``exp`` and shape ``k`` (Poisson if inf)."""
    obs, exp = np.asarray(obs, float), np.asarray(exp, float)
    ok = np.isfinite(obs) & np.isfinite(exp) & (exp > 0)
    o, e = np.where(ok, obs, 0), np.where(ok, exp, 1.0)
    if np.isinf(k):
        lt = stats.poisson.logcdf(o, e)
    else:
        lt = stats.nbinom.logcdf(o, k, k / (k + e))
    return np.where(ok, lt / np.log(10), np.nan)


def power_law_counts(mags, cum_counts, bright_slope: float = 0.6):
    """N(>S) from cumulative counts N(<m) measured at ``mags``, as a function of S = 10^(-0.4 m).

    Interpolated in log N between the measured magnitudes; brighter than the first one it uses
    ``bright_slope`` = dlog10 N / dm (ASSUMPTION: 0.6 is the Euclidean slope that bright galaxy
    counts approach), fainter than the last one the last measured slope.
    """
    m = np.asarray(mags, float)
    lg = np.log10(np.asarray(cum_counts, float))
    faint_slope = (lg[-1] - lg[-2]) / (m[-1] - m[-2])

    def n_brighter(flux):
        mm = -2.5 * np.log10(np.asarray(flux, float))
        out = np.interp(mm, m, lg)
        out = np.where(mm < m[0], lg[0] + bright_slope * (mm - m[0]), out)
        out = np.where(mm > m[-1], lg[-1] + faint_slope * (mm - m[-1]), out)
        return 10.0**out

    return n_brighter


def predicted_ratio(x, n_brighter, mag_lim: float) -> np.ndarray:
    """Image-plane count ratio around a repulsive point lens (``model_prediction``, n = 1)."""
    return exotic_sim.count_ratio(x, n_brighter, 10 ** (-0.4 * mag_lim), n=1.0, sign=-1)


def inject_hole(gal, n_bar, centre, theta_e_cells, ratio_fn, rng) -> np.ndarray:
    """Counts with a predicted hole at ``centre`` (row, col): binomial thinning where the ratio is
    below 1, Poisson additions of n_bar * (ratio - 1) per cell where it is above 1 (simulation)."""
    yy, xx = np.indices(gal.shape)
    x = np.hypot(yy - centre[0], xx - centre[1]) / theta_e_cells
    ratio = np.ones(gal.shape)
    near = x < 5
    with np.errstate(all="ignore"):
        ratio[near] = ratio_fn(np.maximum(x[near], 1e-3))
    ratio = np.nan_to_num(ratio, nan=1.0)  # |mu| -> inf on the critical circle x = 1 (zero measure)
    out = rng.binomial(gal.astype(int), np.clip(ratio, 0, 1)).astype(float)
    out += rng.poisson(n_bar * np.clip(ratio - 1, 0, None))
    return out
