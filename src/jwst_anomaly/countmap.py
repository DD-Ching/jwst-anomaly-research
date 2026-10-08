"""Galaxy count maps on HEALPix pixels for the W5 count-deficit screen (D-063).

A negative-mass lens (n = 1, ε < 0) would leave a deficit of background galaxies inside about its
Einstein radius θ_E (``exotic_sim.count_ratio``; docs/exotic_lensing.md, W5). Testing this needs
galaxy counts with the survey's mask and depth, not individual galaxies, so the survey enters as a
count map (:class:`signatures.CountMapSurvey`): one row per HEALPix ``nest4096`` pixel (0.86′ on a
side, 2.05 × 10⁻⁴ deg²) with the number of selected galaxies and the per-pixel mask and depth.

:class:`LegacySurveysCountMap` aggregates the Legacy Surveys DR10 Tractor catalogue
(``ls_dr10.tractor``) server-side on the NOIRLab Astro Data Lab TAP (``GROUP BY nest4096``), so a
few hundred deg² cost tens of MB instead of a multi-GB catalogue download. The ADQL front end
rejects sub-selects, CASE, SIGN and GROUP BY on expressions (2026-10-08), so each 2° × 2° chunk
runs three queries: selected galaxies, all primary sources (depth, nobs, E(B−V)) and sources with
any ``maskbits`` set. DR10 random catalogues are not on Data Lab; the unmasked fraction of a pixel
is estimated as 1 − n_bad / n_all (**assumption**: the masked fraction of catalogue sources tracks
the masked fraction of area).

The screen helpers (:func:`pixel_table`, :func:`rasterise`, :func:`matched_filter`,
:func:`inject_deficit`) are pure numpy and run on synthetic maps in the tests. Exotic physics is a
hypothesis; a flag is an anomaly to vet, and a null becomes an injection-calibrated limit.
"""

from __future__ import annotations

import math
import os
import threading
import time
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack

from jwst_anomaly import exotic_sim, paths, schema

TAP_URL = "https://datalab.noirlab.edu/tap"
TABLE = "ls_dr10.tractor"
NSIDE = 4096
PIX_AREA_DEG2 = 4 * math.pi * (180 / math.pi) ** 2 / (12 * NSIDE**2)
REFERENCE = (
    "Legacy Surveys DR10 Tractor catalogue (Dey et al. 2019, AJ 157, 168) via NOIRLab Astro Data "
    "Lab TAP " + TAP_URL
)

# Map columns (one row per nest4096 pixel). observed counts of catalogue rows; ``w`` derived.
MAP_COLUMNS = ("pix", "n_gal", "n_all", "n_bad", "depth_r", "nobs_min", "ebv", "w")


def depth_mag(galdepth_ivar) -> np.ndarray:
    """5σ galaxy depth (AB mag) from the Tractor ``galdepth`` inverse variance (nanomaggies⁻²)."""
    iv = np.clip(np.asarray(galdepth_ivar, float), 1e-12, None)
    return 22.5 - 2.5 * np.log10(5.0 / np.sqrt(iv))


@dataclass(frozen=True)
class Region:
    """An RA/Dec box (deg), fetched in ``step`` × ``step`` chunks."""

    name: str
    ra_min: float
    ra_max: float
    dec_min: float
    dec_max: float
    step: float = 2.0

    def chunks(self) -> list[tuple[float, float, float, float]]:
        out = []
        for ra in np.arange(self.ra_min, self.ra_max, self.step):
            for dec in np.arange(self.dec_min, self.dec_max, self.step):
                out.append(
                    (
                        float(ra),
                        float(min(ra + self.step, self.ra_max)),
                        float(dec),
                        float(min(dec + self.step, self.dec_max)),
                    )
                )
        return out

    def area_deg2(self) -> float:
        """Spherical area of the box (deg²), before masking."""
        d = math.radians
        return (
            (180 / math.pi) ** 2
            * d(self.ra_max - self.ra_min)
            * (math.sin(d(self.dec_max)) - math.sin(d(self.dec_min)))
        )


def galaxy_selection(mag_lim: float) -> str:
    """ADQL WHERE clause of the galaxy sample (extended Tractor models, unmasked, dereddened r)."""
    return (
        f"maskbits=0 AND type<>'PSF' AND type<>'DUP' AND dered_mag_r<{mag_lim:.2f} "
        "AND dered_mag_r>0"
    )


def chunk_queries(box: tuple[float, float, float, float], mag_lim: float) -> dict[str, str]:
    """The three ADQL aggregations of one chunk (half-open in RA/Dec, so chunks do not overlap)."""
    r0, r1, d0, d1 = box
    where = f"ra>={r0} AND ra<{r1} AND dec>={d0} AND dec<{d1} AND brick_primary=1"
    return {
        "gal": f"SELECT nest4096, COUNT(*) AS n_gal FROM {TABLE} WHERE {where} AND "
        f"{galaxy_selection(mag_lim)} GROUP BY nest4096",
        "all": f"SELECT nest4096, COUNT(*) AS n_all, AVG(galdepth_r) AS galdepth_r, "
        f"MIN(nobs_g) AS nobs_g, MIN(nobs_r) AS nobs_r, MIN(nobs_z) AS nobs_z, AVG(ebv) AS ebv "
        f"FROM {TABLE} WHERE {where} GROUP BY nest4096",
        "bad": f"SELECT nest4096, COUNT(*) AS n_bad FROM {TABLE} WHERE {where} AND maskbits<>0 "
        "GROUP BY nest4096",
    }


def _run_tap(query: str, retries: int = 3) -> Table:
    import pyvo

    service = pyvo.dal.TAPService(TAP_URL)
    for attempt in range(retries):
        try:
            return service.run_sync(query, maxrec=1_000_000).to_table()
        except Exception:  # network or service error: back off and retry
            if attempt == retries - 1:
                raise
            time.sleep(30 * (attempt + 1))
    raise RuntimeError("unreachable")


def merge_chunk(gal: Table, all_: Table, bad: Table) -> Table:
    """Join the three per-pixel aggregations of one chunk into map rows (``observed`` counts)."""
    pix = np.asarray(all_["nest4096"], np.int64)
    order = np.argsort(pix)
    pix = pix[order]

    def lookup(t: Table, col: str) -> np.ndarray:
        out = np.zeros(pix.size, np.int64)
        if len(t):
            p = np.asarray(t["nest4096"], np.int64)
            i = np.searchsorted(pix, p)
            ok = (i < pix.size) & (pix[np.minimum(i, pix.size - 1)] == p)
            out[i[ok]] = np.asarray(t[col], np.int64)[ok]
        return out

    nobs = np.min(
        [np.asarray(all_[c], float)[order] for c in ("nobs_g", "nobs_r", "nobs_z")], axis=0
    )
    return Table(
        {
            "pix": pix,
            "n_gal": lookup(gal, "n_gal"),
            "n_all": np.asarray(all_["n_all"], np.int64)[order],
            "n_bad": lookup(bad, "n_bad"),
            "depth_r": depth_mag(np.asarray(all_["galdepth_r"], float)[order]),
            "nobs_min": nobs,
            "ebv": np.asarray(all_["ebv"], float)[order],
        },
        meta={
            "provenance": schema.Provenance.OBSERVED.value,
            "source": f"{REFERENCE}: one chunk's per-pixel aggregations (galaxies, all, masked)",
        },
    )


class LegacySurveysCountMap:
    """Legacy Surveys DR10 galaxy counts per nest4096 pixel in one region (CountMapSurvey).

    Chunks are cached as FITS under ``$JWST_ANOMALY_DATA/cache/w5_counts/<region>/``.
    """

    nside = NSIDE

    def __init__(self, region: Region, mag_lim: float = 23.5, cache: Path | None = None):
        self.region = region
        self.mag_lim = mag_lim
        self.name = f"LS-DR10 {region.name} r<{mag_lim}"
        self.cache = cache or paths.cache_dir() / "w5_counts" / f"{region.name}_r{mag_lim:.1f}"

    def chunk_path(self, box) -> Path:
        r0, _, d0, _ = box
        return self.cache / f"chunk_ra{r0:+07.2f}_dec{d0:+06.2f}.fits"

    def fetch_chunk(self, box, run: Callable[[str], Table] = _run_tap) -> Path:
        path = self.chunk_path(box)
        if path.exists():
            return path
        q = chunk_queries(box, self.mag_lim)
        t = merge_chunk(run(q["gal"]), run(q["all"]), run(q["bad"]))
        t.meta["query_gal"] = q["gal"]
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(f".{os.getpid()}.{threading.get_ident()}.tmp.fits")
        t.write(tmp, overwrite=True)
        tmp.replace(path)
        return path

    def missing_chunks(self) -> list:
        return [b for b in self.region.chunks() if not self.chunk_path(b).exists()]

    def count_map(self) -> Table:
        """All cached chunks as one map table (``observed`` counts, ``derived`` w)."""
        missing = self.missing_chunks()
        if missing:
            raise FileNotFoundError(f"{len(missing)} chunks not fetched (run fetch first)")
        parts = []
        for b in self.region.chunks():
            p = Table.read(self.chunk_path(b))
            # the cache key is coarse (region, mag_lim to 0.1, lower-left corner): refuse a chunk
            # fetched with another selection or geometry instead of silently reusing it
            if p.meta.get("query_gal") != chunk_queries(b, self.mag_lim)["gal"]:
                raise ValueError(f"{self.chunk_path(b)} was fetched with another query; refetch it")
            p.meta.clear()
            parts.append(p)
        t = combine_duplicates(vstack(parts))
        return finalize_map(t, source=f"{REFERENCE}; {self.name}; {galaxy_selection(self.mag_lim)}")

    def area_deg2(self) -> float:
        return float(np.sum(self.count_map()["w"]) * PIX_AREA_DEG2)


def combine_duplicates(t: Table) -> Table:
    """Merge rows of one pixel split between chunks (a pixel straddling a chunk edge).

    Counts add; depth and E(B−V) are n_all-weighted means; nobs_min is the minimum. Without
    this, edge pixels keep only one chunk's part and every chunk boundary looks like a deficit.
    """
    pix = np.asarray(t["pix"], np.int64)
    uniq, inv = np.unique(pix, return_inverse=True)
    if uniq.size == pix.size:
        out = t.copy()
        out.sort("pix")
        out.meta.update(
            provenance=schema.Provenance.OBSERVED.value,
            source=f"{REFERENCE}: chunk rows (no pixel split at chunk edges)",
        )
        return out
    n_all = np.asarray(t["n_all"], float)

    def add(col):
        return np.bincount(inv, weights=np.asarray(t[col], float), minlength=uniq.size)

    wsum = np.maximum(add("n_all"), 1e-300)
    nobs = np.full(uniq.size, np.inf)
    np.minimum.at(nobs, inv, np.asarray(t["nobs_min"], float))
    depth_ivar = 25.0 * 10 ** (0.8 * (np.asarray(t["depth_r"], float) - 22.5))
    ivar = np.bincount(inv, weights=n_all * depth_ivar, minlength=uniq.size) / wsum
    return Table(
        {
            "pix": uniq,
            "n_gal": add("n_gal").astype(np.int64),
            "n_all": add("n_all").astype(np.int64),
            "n_bad": add("n_bad").astype(np.int64),
            "depth_r": depth_mag(ivar),
            "nobs_min": nobs,
            "ebv": np.bincount(inv, weights=n_all * np.asarray(t["ebv"], float)) / wsum,
        },
        meta={
            "provenance": schema.Provenance.OBSERVED.value,
            "source": f"{REFERENCE}: chunk rows with pixels split at chunk edges summed",
        },
    )


def finalize_map(t: Table, source: str) -> Table:
    """Add the unmasked fraction ``w`` = 1 − n_bad/n_all and the provenance meta."""
    n_all = np.asarray(t["n_all"], float)
    t["w"] = np.where(n_all > 0, 1.0 - np.asarray(t["n_bad"], float) / np.maximum(n_all, 1), 0.0)
    t.meta["provenance"] = schema.Provenance.OBSERVED.value
    t.meta["source"] = source
    t.meta["nside"] = NSIDE
    t.meta["note"] = "n_* observed catalogue-row counts per nest4096 pixel; w derived (assumption)"
    return t


# ----------------------------------------------------------------------------- prediction


def power_counts(alpha: float) -> Callable[[np.ndarray], np.ndarray]:
    """N(>S) ∝ S^-alpha (ASSUMPTION; alpha = 2.5 × d log10 N / dm)."""
    return lambda s: np.asarray(s, float) ** (-alpha)


def tabulated_counts(
    mag: np.ndarray, n_brighter: np.ndarray, bright_slope: float = 0.6
) -> Callable:
    """N(>S) from a measured cumulative count table N(< m) (per deg², observed or derived).

    S is flux relative to the limit (S = 10^(-0.4 (m - m_lim))) — only ratios matter. Brighter
    than the table, d log N / dm = ``bright_slope`` (default 0.6, Euclidean: an ASSUMPTION, needed
    because a few bright galaxies per 100 deg² cannot be counted in a pilot area); fainter, the
    last measured slope continues (ASSUMPTION).
    """
    m = np.asarray(mag, float)
    lg = np.log10(np.asarray(n_brighter, float))
    if np.any(np.diff(m) <= 0) or np.any(~np.isfinite(lg)):
        raise ValueError("mag must increase and counts be positive")
    s_lo, s_hi = bright_slope, (lg[-1] - lg[-2]) / (m[-1] - m[-2])

    def n_of_s(s, m_lim: float = float(m[-1])):
        mm = m_lim - 2.5 * np.log10(np.asarray(s, float))
        out = np.interp(mm, m, lg)
        out = np.where(mm < m[0], lg[0] + s_lo * (mm - m[0]), out)
        out = np.where(mm > m[-1], lg[-1] + s_hi * (mm - m[-1]), out)
        return 10.0**out

    return n_of_s


def deficit_profile(
    x, counts: Callable, flux_limit: float = 1.0, mu_max: float = 30.0
) -> np.ndarray:
    """Predicted N_obs/N̄ at image radius x (θ_E units) for n = 1, ε < 0 (``model_prediction``).

    Both images lie on the source's side, so every image radius x > 0 is reached:
    the inner image (x < 1) of sources far outside, demagnified as |μ| ≈ x⁴. This is
    ``exotic_sim.count_ratio`` with |μ| capped at ``mu_max`` near the radial critical curve x = 1,
    where extended galaxies (not point sources) limit the magnification (ASSUMPTION; it matters
    only for a narrow ring around x = 1).
    """
    x = np.maximum(np.asarray(x, float), 1e-4)
    mu = np.minimum(np.abs(exotic_sim.image_plane_magnification(x, 1.0, -1)), mu_max)
    return counts(flux_limit / mu) / mu / counts(flux_limit)


def expected_missing(theta_e_arcmin: float, density_per_arcmin2: float, profile, x_max=3.0):
    """Expected net missing galaxies within x_max θ_E (positive = deficit), and the N̄ there."""
    x = np.linspace(1e-3, x_max, 3000)
    dx = x[1] - x[0]
    ring = 2 * np.pi * x * dx * theta_e_arcmin**2 * density_per_arcmin2
    return float(np.sum((1 - profile(x)) * ring)), float(np.sum(ring))


# ----------------------------------------------------------------------------- flat raster


@dataclass
class Raster:
    """A flat grid (sinusoidal, equal-area) of cell size ``cell`` arcmin around (ra0, dec0)."""

    ra0: float
    dec0: float
    cell: float
    nx: int
    ny: int
    x0: float  # arcmin of the grid's first column edge
    y0: float

    def xy(self, ra, dec) -> tuple[np.ndarray, np.ndarray]:
        ra = np.asarray(ra, float)
        dec = np.asarray(dec, float)
        dra = (ra - self.ra0 + 180.0) % 360.0 - 180.0
        return dra * np.cos(np.radians(dec)) * 60.0, (dec - self.dec0) * 60.0

    def radec(self, x, y) -> tuple[np.ndarray, np.ndarray]:
        dec = self.dec0 + np.asarray(y, float) / 60.0
        ra = self.ra0 + np.asarray(x, float) / 60.0 / np.cos(np.radians(dec))
        return ra % 360.0, dec

    def centres(self) -> tuple[np.ndarray, np.ndarray]:
        x = self.x0 + (np.arange(self.nx) + 0.5) * self.cell
        y = self.y0 + (np.arange(self.ny) + 0.5) * self.cell
        return np.meshgrid(x, y)


def make_raster(region: Region, cell: float = 0.25) -> Raster:
    ra0 = 0.5 * (region.ra_min + region.ra_max)
    dec0 = 0.5 * (region.dec_min + region.dec_max)
    half_w = 0.5 * (region.ra_max - region.ra_min) * 60.0 + cell
    half_h = 0.5 * (region.dec_max - region.dec_min) * 60.0 + cell
    nx = int(math.ceil(2 * half_w / cell))
    ny = int(math.ceil(2 * half_h / cell))
    return Raster(ra0, dec0, cell, nx, ny, -half_w, -half_h)


def pixel_index_image(raster: Raster, pix_sorted: np.ndarray) -> np.ndarray:
    """Row index into the (sorted) map for every raster cell; −1 where the pixel is absent."""
    from astropy import units as u
    from astropy_healpix import HEALPix

    hp = HEALPix(nside=NSIDE, order="nested")
    xs, ys = raster.centres()
    idx = np.full(xs.shape, -1, np.int64)
    for j in range(0, raster.ny, 256):
        ra, dec = raster.radec(xs[j : j + 256], ys[j : j + 256])
        p = hp.lonlat_to_healpix(ra * u.deg, dec * u.deg)
        i = np.searchsorted(pix_sorted, p)
        i = np.minimum(i, pix_sorted.size - 1)
        idx[j : j + 256] = np.where(pix_sorted[i] == p, i, -1)
    return idx


def rasterise(values: np.ndarray, idx: np.ndarray, per_cell: bool = False) -> np.ndarray:
    """Paint per-pixel values on the raster; counts become per-cell density when ``per_cell``.

    With ``per_cell`` a pixel's count is spread evenly over the raster cells inside it (the
    number of cells per pixel is counted, so totals are conserved).
    """
    v = np.asarray(values, float)
    out = np.zeros(idx.shape)
    ok = idx >= 0
    if per_cell:
        ncell = np.bincount(idx[ok], minlength=v.size).astype(float)
        out[ok] = v[idx[ok]] / np.maximum(ncell[idx[ok]], 1)
    else:
        out[ok] = v[idx[ok]]
    return out


# ----------------------------------------------------------------------------- matched filter


def _fft_conv(img: np.ndarray, ker: np.ndarray) -> np.ndarray:
    from scipy.signal import fftconvolve

    return fftconvolve(img, ker, mode="same")


def block_sum(img: np.ndarray, f: int) -> np.ndarray:
    """Sum ``f`` × ``f`` blocks (the trailing partial rows/columns are dropped)."""
    if f == 1:
        return img
    ny, nx = (img.shape[0] // f) * f, (img.shape[1] // f) * f
    return img[:ny, :nx].reshape(ny // f, f, nx // f, f).sum(axis=(1, 3))


def kernels(theta_e_cells: float, profile_fn, x_max: float = 2.5):
    """Top-hat W (r ≤ x_max θ_E), template T = profile − 1 inside it, and the inner disc x ≤ 1.

    T is averaged over 3 × 3 sub-samples per cell to follow the steep core.
    """
    r_max = x_max * theta_e_cells
    h = int(math.ceil(r_max)) + 1
    yy, xx = np.mgrid[-h : h + 1, -h : h + 1].astype(float)
    sub = (np.arange(3) - 1) / 3.0
    t = np.zeros_like(xx)
    for dy in sub:
        for dx in sub:
            t += profile_fn(np.hypot(xx + dx, yy + dy) / theta_e_cells) - 1.0
    t /= 9.0
    r = np.hypot(xx, yy)
    w = (r <= r_max).astype(float)
    inner = (r <= theta_e_cells).astype(float)
    return w, t * w, inner


class MatchedFilter:
    """Local least-squares fit D = m a (1 + A T) around every cell of a masked count image.

    ``weight`` m is the unmasked fraction per cell (0 outside the footprint); ``density`` D the
    galaxies per cell, counted in the unmasked part only. A = 1 is the predicted deficit
    (``model_prediction``), A = 0 no lens. The weight convolutions are computed once, so
    re-filtering an injected map costs two FFT convolutions. ``sigma`` is the Poisson error of A for
    independent cells (an ASSUMPTION: cells inside one HEALPix pixel are correlated and the galaxy
    field is clustered, so significances are recalibrated on an independent null region).
    ``cover`` and ``full`` are the unmasked fractions of the inner disc and of the whole aperture.
    """

    def __init__(self, weight: np.ndarray, theta_e_cells: float, profile_fn, x_max: float = 2.5):
        self.w_k, self.t_k, inner = kernels(theta_e_cells, profile_fn, x_max)
        m2 = weight**2
        self.weight = weight
        self.s11 = _fft_conv(m2, self.w_k)
        self.s12 = _fft_conv(m2, self.t_k)
        self.s22 = _fft_conv(m2, self.t_k**2)
        self.det = self.s11 * self.s22 - self.s12**2
        self.cover = _fft_conv(weight, inner) / inner.sum()
        self.full = _fft_conv(weight, self.w_k) / self.w_k.sum()

    def __call__(self, density: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        md = self.weight * density
        y1 = _fft_conv(md, self.w_k)
        y2 = _fft_conv(md, self.t_k)
        with np.errstate(divide="ignore", invalid="ignore"):
            a = (y1 * self.s22 - y2 * self.s12) / self.det
            b = (self.s11 * y2 - self.s12 * y1) / self.det
            amp = b / a
            sigma = np.sqrt(np.abs(a) * self.s11 / self.det) / np.abs(a)
        return amp, sigma


def find_peaks(
    z: np.ndarray, ok: np.ndarray, size: int, z_min: float
) -> tuple[np.ndarray, np.ndarray]:
    """(row, col) of local maxima of z (within ``size`` cells) where ``ok`` and z ≥ z_min."""
    from scipy.ndimage import maximum_filter

    zz = np.where(ok & np.isfinite(z), z, -np.inf)
    peak = (
        (zz == maximum_filter(zz, size=max(3, size), mode="constant", cval=-np.inf))
        & np.isfinite(zz)
        & (zz >= z_min)
    )
    return np.nonzero(peak)


# ----------------------------------------------------------------------------- injection


def inject_deficit(
    counts: np.ndarray,
    pix_xy: tuple[np.ndarray, np.ndarray],
    centre_xy: tuple[float, float],
    theta_e: float,
    profile_fn,
    rng: np.random.Generator,
    x_max: float = 3.0,
) -> np.ndarray:
    """Apply the predicted count ratio around one centre to per-pixel counts (``simulated``).

    Where the ratio is < 1, each galaxy survives with probability ratio (binomial thinning, i.e.
    sources demagnified below the limit or moved out); where it is > 1 the pixel gains
    Poisson((ratio − 1) n) galaxies. The ratio is averaged over the pixel (2 × 2 sub-samples at
    0.86′ / 4 offsets). Positions and θ_E in the same units (arcmin).
    """
    x, y = pix_xy
    r = np.hypot(x - centre_xy[0], y - centre_xy[1])
    near = r <= x_max * theta_e + 1.0
    out = np.array(counts, dtype=np.int64, copy=True)
    if not near.any():
        return out
    side = math.sqrt(PIX_AREA_DEG2) * 60.0
    ratio = np.zeros(int(near.sum()))
    for dx in (-0.25, 0.25):
        for dy in (-0.25, 0.25):
            rr = np.hypot(x[near] + dx * side - centre_xy[0], y[near] + dy * side - centre_xy[1])
            ratio += np.where(rr / theta_e <= x_max, profile_fn(rr / theta_e), 1.0)
    ratio /= 4.0
    n = out[near]
    thin = rng.binomial(n, np.clip(ratio, 0.0, 1.0))
    extra = rng.poisson(np.clip(ratio - 1.0, 0.0, None) * n)
    out[near] = thin + extra
    return out
