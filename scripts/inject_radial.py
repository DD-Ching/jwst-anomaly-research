"""Injection-recovery of W1 negative-mass lenses through the ``radial`` screen (D-049).

W1 (docs/exotic_lensing.md, D-047) is a compact lens with ε < 0 and n = 1: no images of sources
within 2 θ_E (the umbra), and two radially stretched images on the source's side of every source
beyond. This script paints such lenses into a field's real catalogue, runs the screen's own arc
selection, grid, null draws and convergence search (``exotic_screens``: ``radial_candidates``,
``radial_grid``, ``anti_window_draw``, ``line_counts``, ``convergence_peaks``; same defaults), and
measures how often a lens is recovered. With the real screens null, the recovery efficiencies give
95 % upper limits on the surface density of such lenses, per lens mass.

Everything injected is ``simulated``; efficiencies and limits are ``derived`` from simulated signals
in real catalogues. A limit says how many such lenses the screen would have seen; it is not a
detection of anything. The lens mass |M| is the parameter; θ_E per source follows from the
point-lens formula (``model_prediction``).

Injection (all ASSUMPTIONs, documented in docs/exotic_limits.md):
- lens of mass |M| at the cluster redshift, at a uniform random position in the screened footprint
  (1" grid points inside ``--max-radius`` with a catalogue source within 4");
- lensed sources are the catalogue's own rows behind the lens position: no photo-z, or a photo-z
  the screen calls background (z160 > z_l + 0.1); not bright point sources (stars). Each row's
  θ_E uses its photo-z (z_phot), or z_s = 2 without one;
- rows with beta < 2 θ_E (umbra) are removed; rows with 2 <= beta <= 4 θ_E are replaced by their
  images (``exotic_sim.inject_images``); rows beyond 4 θ_E are left as they are (D-047 range);
- each image keeps its source's columns; its shape is the source's second-moment ellipse, PSF
  deconvolved, mapped by the signed lens Jacobian and re-convolved; magnitude - 2.5 log10 |mu|,
  isophotal area x |mu| and S/N x sqrt(|mu|) (background-limited);
- the two images of a source are painted as one blended row when they are closer than the sum of
  their isophotal equivalent radii sqrt(area / pi): fluxes add, second moments combine
  flux-weighted about the joint centroid;
- rows with S/N < 5 are not detected (dropped);
- the PSF sigma is the 1st percentile of semimajor x (1 - ellipticity) among S/N > 50 sources.

Recovery: a convergence peak with p_random < 0.05 (the screen's definition) within 2" of the
injected centre. The null is the screen's: 200 draws of each arc's angle inside its ``anti``
window. Trials run in batches of 10; each batch draws its own independent 200-draw null for the
real arcs, and each trial removes the replaced arcs from it and adds fresh draws for the injected
ones (only the touched grid blocks are recomputed).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from astropy.table import Table, vstack

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lens_consistency as lc  # noqa: E402

from jwst_anomaly import exotic_sim, paths, schema  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "exotic_screens", Path(__file__).resolve().parent / "exotic_screens.py"
)
es = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(es)

BETA_MIN, BETA_MAX = 2.0, 4.0  # lensed-source range in theta_E (D-047 recommendation; ASSUMPTION)
P_RECOVER = 0.05  # the screen's p_random threshold
RECOVER_TOL = 2.0  # arcsec between the peak and the injected centre (ASSUMPTION)
FOOTPRINT_STEP = 1.0  # arcsec grid for the screened footprint (ASSUMPTION)
FOOTPRINT_RADIUS = 4.0  # a grid point is screened if a catalogue source lies within (ASSUMPTION)
POISSON_UL_95 = 2.996  # -ln(0.05): 95 % upper limit on a Poisson mean with zero events
Z_SOURCE = 2.0  # source redshift of rows without photo-z, and of the theta_E label (ASSUMPTION)
SNR_DETECT = 5.0  # painted rows below this isophotal S/N are not detected (ASSUMPTION)
BATCH = 10  # trials per independent null realisation
BLOCK = 32  # grid cells per block for the incremental null maximum
MASSES = (2e10, 2e11, 2e12, 8e12, 2e13)  # |M| in M_sun (theta_E ~ 0.3"-10" at z_l 0.4, z_s 2)
SCREEN_DEFAULTS = es.radial_defaults()

# The inputs each field's null screen used (docs/fields/*.md, docs/exotic_lensing.md). Paths are
# relative to the data root. Photo-z kinds: "eazy" (DJA zout) or "canucs" (CANUCS DR1 catalogue,
# converted as in docs/fields/macs1149.md).
_MAST = "cache/mast/{0}/{0}_cat.ecsv"
_CANUCS = "cache/external/canucs/hlsp_canucs_jwst-hst_multi_{0}-clu_multi_v1_photometry-cat.fits.gz"
FIELDS = {
    "smacs0723": {
        "model": "smacs0723-iclv2",
        "catalog": _MAST.format("jw02736-o001_t001_nircam_clear-f200w"),
        "photoz": ("eazy", "cache/dja/photoz-v7.4/smacs0723-grizli-v7.4-fix.eazypy.zout.fits"),
        "max_radius": 60.0,
    },
    "elgordo": {
        "model": "elgordo-caminha23",
        "catalog": _MAST.format("jw01176-o241_t012_nircam_clear-f200w"),
        "photoz": (
            "eazy",
            "cache/dja/photoz-elgordo-v7.0/elgordo-grizli-v7.0-fix.eazypy.zout.fits",
        ),
        "max_radius": 110.0,
    },
    "abell2744": {
        "model": "abell2744-bergamini23",
        "catalog": _MAST.format("jw02561-o001_t003_nircam_clear-f200w"),
        "photoz": (
            "eazy",
            "cache/dja/photoz-abell2744-v7.2/abell2744clu-grizli-v7.2-fix.eazypy.zout.fits",
        ),
        "max_radius": 120.0,
    },
    "macs0416": {
        "model": "macs0416-cats",
        "catalog": _MAST.format("jw01208-o004_t002_nircam_clear-f200w"),
        "photoz": ("canucs", _CANUCS.format("macs0416")),
        "max_radius": 100.0,
    },
    "macs1149": {
        "model": "macs1149-cats",
        "catalog": _MAST.format("jw01208-o008_t004_nircam_clear-f200w"),
        "photoz": ("canucs", _CANUCS.format("macs1149")),
        "max_radius": 100.0,
    },
    "macs0717": {
        "model": "macs0717-cats",
        "catalog": _MAST.format("jw06882-o029_t063_nircam_clear-f200w"),
        "photoz": None,
        "max_radius": 120.0,
    },
    "abell370": {
        "model": "abell370-cats",
        "catalog": _MAST.format("jw01208-o002_t001_nircam_clear-f200w"),
        "photoz": ("canucs", _CANUCS.format("a370")),
        "max_radius": 100.0,
        "spike_stars": True,  # Gaia DR3 (scripts/gaia_stars.py), as in docs/fields/abell370.md
        "centre": (39.97134, -1.58226),
    },
    "abells1063": {
        "model": "abells1063-cats",
        "catalog": _MAST.format("jw03293-o001_t001_nircam_clear-f200w"),
        "photoz": None,  # the DJA v7.5 zout is a 75 MB tarball; the null holds with or without it
        "max_radius": 100.0,
    },
}


# ----------------------------------------------------------------------------- physics helpers


def _cosmology(cosmology):
    if cosmology is None:
        from astropy.cosmology import Planck18 as cosmology
    return cosmology


def theta_e_arcsec(mass_msun: float, z_l: float, z_s, cosmology=None) -> np.ndarray:
    """Einstein (scale) radius in arcsec of an n = 1 lens of mass |M| at ``z_l`` for sources at
    ``z_s`` (array; NaN where z_s <= z_l). ``model_prediction``; Planck18 by default."""
    cosmo = _cosmology(cosmology)
    z_s = np.atleast_1d(np.asarray(z_s, float))
    out = np.full(z_s.shape, np.nan)
    ok = np.isfinite(z_s) & (z_s > z_l)
    if ok.any():
        d_l = cosmo.angular_diameter_distance(z_l).si.value
        d_s = cosmo.angular_diameter_distance(z_s[ok]).si.value
        d_ls = cosmo.angular_diameter_distance(z_l, z_s[ok]).si.value
        rad = exotic_sim.einstein_radius(
            exotic_sim.eps_bar_point_mass(mass_msun), 1, d_l, d_s, d_ls
        )
        out[ok] = np.degrees(rad) * 3600.0
    return out


def theta_e_to_mass(theta_e, z_l: float, z_s: float = Z_SOURCE, cosmology=None):
    """|M| in M_sun of an n = 1 lens with Einstein (scale) radius ``theta_e`` arcsec (inverse of
    :func:`theta_e_arcsec`)."""
    ref = theta_e_arcsec(1.0, z_l, z_s, cosmology)[0]
    return (np.asarray(theta_e, float) / ref) ** 2


def surface_density_limit(efficiency, area_deg2, n_ul: float = POISSON_UL_95) -> float:
    """95 % upper limit (deg^-2) with zero detections: n_ul / sum(efficiency x area)."""
    exposure = float(np.sum(np.asarray(efficiency, float) * np.asarray(area_deg2, float)))
    return n_ul / exposure if exposure > 0 else float("inf")


def _ellipse_cov(a, b, pa_deg):
    """Second-moment matrices (East, North basis) of ellipses with sigmas ``a`` >= ``b``."""
    pa = np.deg2rad(np.asarray(pa_deg, float))
    u = np.stack([np.sin(pa), np.cos(pa)], axis=-1)  # major axis, PA east of north
    v = np.stack([np.cos(pa), -np.sin(pa)], axis=-1)
    a2, b2 = np.asarray(a, float) ** 2, np.asarray(b, float) ** 2
    return (
        a2[:, None, None] * u[:, :, None] * u[:, None, :]
        + b2[:, None, None] * v[:, :, None] * v[:, None, :]
    )


def cov_shape(cov) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(semimajor sigma, ellipticity 1 - b/a, PA east of north) of second-moment matrices."""
    w, vec = np.linalg.eigh(cov)
    major = vec[:, :, 1]  # eigh sorts ascending
    pa = np.mod(np.degrees(np.arctan2(major[:, 0], major[:, 1])), 180.0)
    return np.sqrt(w[:, 1]), 1.0 - np.sqrt(np.maximum(w[:, 0], 0.0) / w[:, 1]), pa


def lensed_cov(a, b, pa_deg, dx, dy, lam_r, lam_t, psf_sigma: float) -> np.ndarray:
    """Image second moments (px^2, East/North) of sources (sigmas ``a``, ``b``, ``pa_deg``).

    ``dx``, ``dy``: image offsets from the lens (East, North); ``lam_r``, ``lam_t``: signed Jacobian
    eigenvalues at the image. The PSF (isotropic, ``psf_sigma`` px) is removed from the source
    moments, the intrinsic ellipse mapped by A^-1, and the PSF added back (ASSUMPTION: Gaussian
    moments, isolated lens)."""
    cov = _ellipse_cov(a, b, pa_deg)
    eye = np.eye(2)
    floor = (0.05 * psf_sigma) ** 2
    w, vec = np.linalg.eigh(cov - psf_sigma**2 * eye)
    cov_int = vec @ (np.maximum(w, floor)[..., None] * np.swapaxes(vec, -1, -2))
    r = np.stack([np.asarray(dx, float), np.asarray(dy, float)], axis=-1)
    r = r / np.linalg.norm(r, axis=-1, keepdims=True)
    t = np.stack([-r[:, 1], r[:, 0]], axis=-1)
    inv = (1.0 / np.asarray(lam_r, float))[:, None, None] * r[:, :, None] * r[:, None, :] + (
        1.0 / np.asarray(lam_t, float)
    )[:, None, None] * t[:, :, None] * t[:, None, :]
    return inv @ cov_int @ np.swapaxes(inv, -1, -2) + psf_sigma**2 * eye


def lensed_shapes(a, b, pa_deg, dx, dy, lam_r, lam_t, psf_sigma: float):
    """:func:`lensed_cov` as (semimajor, ellipticity, pa_deg)."""
    return cov_shape(lensed_cov(a, b, pa_deg, dx, dy, lam_r, lam_t, psf_sigma))


def background_mask(shapes: Table, z_lens: float, z_margin: float = 0.1) -> np.ndarray:
    """Rows a lens at ``z_lens`` can lens (``exotic_screens.lensable_mask``)."""
    return es.lensable_mask(shapes, z_lens, z_margin)


def source_redshift(shapes: Table, z_default: float = Z_SOURCE) -> np.ndarray:
    """Each row's source redshift for θ_E: its photo-z, or ``z_default`` without one."""
    if "z_phot" not in shapes.colnames:
        return np.full(len(shapes), z_default)
    z = np.asarray(shapes["z_phot"], float)
    return np.where(np.isfinite(z), z, z_default)


def _merge_pairs(img: Table, area_src: np.ndarray, pixel_scale: float):
    """Indices of image pairs (outer, inner) of one source whose isophotes overlap: separation
    below the sum of the equivalent radii sqrt(area |mu| / pi)."""
    sid = np.asarray(img["source_id"], int)
    k = np.asarray(img["image"], int)
    outer = {s: i for i, s in enumerate(sid) if k[i] == 0}
    pairs = []
    for i, s in enumerate(sid):
        if k[i] != 1 or s not in outer:
            continue
        j = outer[s]
        sep = np.hypot(img["dx"][i] - img["dx"][j], img["dy"][i] - img["dy"][j])
        mu = np.abs(np.asarray(img["mu"], float)[[j, i]])
        r_eq = np.sqrt(area_src[[j, i]] * mu / np.pi) * pixel_scale
        if sep < r_eq.sum():
            pairs.append((j, i))
    return pairs


def paint_lens(
    shapes: Table,
    xs: np.ndarray,
    ys: np.ndarray,
    background: np.ndarray,
    theta_row: np.ndarray,
    model,
    x_l: float,
    y_l: float,
    psf_sigma: float,
    first_label: int,
    pixel_scale: float,
    snr_floor: float = SNR_DETECT,
) -> tuple[np.ndarray, Table, dict]:
    """Paint one W1 lens at model-frame (``x_l``, ``y_l``) into ``shapes`` (``simulated`` rows).

    ``xs``, ``ys``: model-frame positions of ``shapes`` (x West, y North); ``theta_row``: each row's
    θ_E (arcsec) for this lens mass (NaN: not lensable). Returns the mask of rows to keep, the image
    rows (same columns, new labels from ``first_label``) and counts."""
    dx = -(xs - x_l)  # East
    dy = ys - y_l
    with np.errstate(invalid="ignore", divide="ignore"):
        beta = np.hypot(dx, dy) / theta_row
    hit = background & np.isfinite(beta) & (beta < BETA_MAX)
    lensed = np.flatnonzero(hit & (beta >= BETA_MIN))
    keep = ~hit
    info = {"n_umbra": int(np.sum(hit & (beta < BETA_MIN))), "n_lensed": len(lensed)}
    meta = {
        "provenance": schema.Provenance.SIMULATED.value,
        "source": "exotic_sim.inject_images(n=1, sign=-1) painted on "
        f"{shapes.meta.get('source', 'catalogue')}",
    }
    empty = shapes[:0].copy()
    empty.meta = meta
    if not len(lensed):
        return keep, empty, info | {"n_images": 0, "n_merged": 0}
    th = theta_row[lensed]
    img = exotic_sim.inject_images(
        dx[lensed] / th, dy[lensed] / th, 1.0, n=1.0, sign=-1, source_id=lensed
    )
    if not len(img):
        return keep, empty, info | {"n_images": 0, "n_merged": 0}
    k = np.asarray(img["source_id"], int)
    scale = theta_row[k]
    img["dx"] = np.asarray(img["dx"], float) * scale
    img["dy"] = np.asarray(img["dy"], float) * scale
    lam_r = 1.0 / np.asarray(img["stretch_r"], float)
    lam_t = 1.0 / np.asarray(img["stretch_t"], float)
    # signs: lam_t > 0 always for eps < 0; lam_r < 0 for the inner image (x < 1)
    lam_r = np.where(np.asarray(img["image"]) == 1, -lam_r, lam_r)
    a = np.asarray(shapes["semimajor_px"], float)[k]
    b = a * (1.0 - np.asarray(shapes["ellipticity"], float)[k])
    pa_src = np.asarray(shapes["pa_obs"], float)[k]
    dxi, dyi = np.asarray(img["dx"], float), np.asarray(img["dy"], float)
    cov = lensed_cov(a, b, pa_src, dxi, dyi, lam_r, lam_t, psf_sigma)
    mu = np.abs(np.asarray(img["mu"], float))
    area_src = np.asarray(shapes["area_px"], float)[k]
    # blended pairs: one row at the flux-weighted centroid with the combined second moments
    pairs = _merge_pairs(img, area_src, pixel_scale)
    drop = []
    for j, i in pairs:
        w = mu[[j, i]] / mu[[j, i]].sum()
        cx, cy = w @ dxi[[j, i]], w @ dyi[[j, i]]
        c = np.zeros((2, 2))
        for m, wm in zip((j, i), w, strict=True):
            d = np.array([dxi[m] - cx, dyi[m] - cy]) / pixel_scale
            c += wm * (cov[m] + np.outer(d, d))
        cov[j], dxi[j], dyi[j], mu[j] = c, cx, cy, mu[[j, i]].sum()
        drop.append(i)
    sel = np.setdiff1d(np.arange(len(img)), drop)
    rows = shapes[k[sel]]
    semi, ell, pa = cov_shape(cov[sel])
    mu = mu[sel]
    rows["semimajor_px"] = semi
    rows["ellipticity"] = ell
    rows["pa_obs"] = pa
    rows["mag"] = np.asarray(rows["mag"], float) - 2.5 * np.log10(mu)
    rows["area_px"] = np.asarray(rows["area_px"], float) * mu
    rows["snr"] = np.asarray(rows["snr"], float) * np.sqrt(mu)
    if "is_extended" in rows.colnames:  # a lensed image is never a star for the spike veto
        rows["is_extended"] = True
    ra, dec = model.to_sky(x_l - dxi[sel], y_l + dyi[sel])
    rows["ra"], rows["dec"] = ra, dec
    rows["label"] = first_label + np.arange(len(rows))
    rows = rows[np.asarray(rows["snr"], float) >= snr_floor]
    rows.meta = meta
    return keep, rows, info | {"n_images": len(rows), "n_merged": len(pairs)}


def screened_footprint(xs, ys, max_radius: float, radius: float = FOOTPRINT_RADIUS):
    """Model-frame grid points (``FOOTPRINT_STEP`` arcsec) inside ``max_radius`` with a catalogue
    source within ``radius``: the area the screen looked at. Returns (gx, gy, area_arcsec2)."""
    from scipy.spatial import cKDTree

    step = FOOTPRINT_STEP
    g = np.arange(-max_radius + step / 2, max_radius, step)
    gx, gy = np.meshgrid(g, g)
    gx, gy = gx.ravel(), gy.ravel()
    inside = np.hypot(gx, gy) <= max_radius
    gx, gy = gx[inside], gy[inside]
    d, _ = cKDTree(np.column_stack([xs, ys])).query(np.column_stack([gx, gy]))
    ok = d <= radius
    return gx[ok], gy[ok], float(ok.sum() * step**2)


def footprint_border(xs, ys, max_radius: float, radii=(4.0, 5.0, 6.0, 8.0)) -> dict:
    """Measured border excess of the footprint: area(r) grows by about perimeter x dr once the
    interior holes are filled, so a linear fit over ``radii`` extrapolated to r = 0 gives the area
    enclosed by the outermost sources; the excess at ``FOOTPRINT_RADIUS`` is an upper bound on the
    overestimate (the true edge lies beyond the outermost sources)."""
    radii = tuple(sorted(set(radii) | {FOOTPRINT_RADIUS}))
    areas = [screened_footprint(xs, ys, max_radius, r)[2] for r in radii]
    slope, a0 = np.polyfit(radii, areas, 1)
    a4 = areas[radii.index(FOOTPRINT_RADIUS)]
    return {
        "radii": list(radii),
        "areas_arcsec2": areas,
        "area_r0_arcsec2": float(a0),
        # clipped to [0, 0.99]: a flat or concave area(r) (footprint cut by max_radius) must not
        # give a negative or vanishing corrected area
        "excess_frac": float(np.clip((a4 - a0) / a4, 0.0, 0.99)),
    }


# ----------------------------------------------------------------------------- the screen


class RadialInjector:
    """The radial screen on one field: real arcs' lines cached, null drawn per batch."""

    def __init__(self, model, shapes: Table, args, extra: Table | None = None):
        self.model, self.shapes, self.args, self.extra = model, shapes, args, extra
        cand, self.base_counts_info = es.radial_candidates(model, shapes, args, extra)
        self.base = cand
        self.gx, self.gy = es.radial_grid(args.max_radius, args.grid_arcsec)
        self.g = self.gx[0]
        cx, cy = model.to_frame(cand["ra"], cand["dec"])
        self.base_xy = (np.atleast_1d(cx), np.atleast_1d(cy))
        self.pa_pred = np.asarray(cand["pa_pred_z2"], float)
        self.counts = es.line_counts(
            cx, cy, cand["pa_obs"], self.gx, self.gy, args.line_tol_arcsec, args.max_len_arcsec
        ).astype(np.int16)
        self.label_index = {int(v): i for i, v in enumerate(cand["label"])}
        xs, ys = model.to_frame(shapes["ra"], shapes["dec"])
        self.xs, self.ys = np.asarray(xs), np.asarray(ys)
        self.background = background_mask(shapes, model.z_lens)
        self.next_label = int(np.max(shapes["label"])) + 1
        # the screen's own null (its seed and draw order) for the real-field numbers
        self.draw_null(np.random.default_rng(args.seed))
        self.screen_max_rand = self.max_rand.copy()

    # --- grid windows ---------------------------------------------------------------------
    def _window(self, x, y) -> tuple[int, int, int, int]:
        r = self.args.max_len_arcsec + self.args.line_tol_arcsec
        step, g0, n = self.args.grid_arcsec, self.g[0], len(self.g)
        i0 = max(int(np.floor((y - r - g0) / step)), 0)
        i1 = min(int(np.ceil((y + r - g0) / step)) + 1, n)
        j0 = max(int(np.floor((x - r - g0) / step)), 0)
        j1 = min(int(np.ceil((x + r - g0) / step)) + 1, n)
        return i0, i1, j0, j1

    def _add(self, grid, x, y, pa, sign: int, origin=(0, 0)) -> None:
        """Add (``sign`` +1) or remove (-1) the lines of arcs at (x, y, pa) to ``grid``, a block of
        the full grid starting at cell ``origin``. Exactly ``line_counts`` restricted to each arc's
        window (a line never reaches beyond max_len + tol)."""
        a = self.args
        oi, oj = origin
        ni, nj = grid.shape
        for xi, yi, pi in zip(x, y, pa, strict=True):
            i0, i1, j0, j1 = self._window(xi, yi)
            i0, i1, j0, j1 = max(i0, oi), min(i1, oi + ni), max(j0, oj), min(j1, oj + nj)
            if i0 >= i1 or j0 >= j1:
                continue
            sub = es.line_counts(
                [xi],
                [yi],
                [pi],
                self.gx[i0:i1, j0:j1],
                self.gy[i0:i1, j0:j1],
                a.line_tol_arcsec,
                a.max_len_arcsec,
            )
            grid[i0 - oi : i1 - oi, j0 - oj : j1 - oj] += sign * sub.astype(np.int16)

    def _block_max(self, grid) -> np.ndarray:
        n = grid.shape[0]
        nb = -(-n // BLOCK)
        pad = np.full((nb * BLOCK, nb * BLOCK), np.iinfo(np.int16).min, np.int16)
        pad[:n, :n] = grid
        return pad.reshape(nb, BLOCK, nb, BLOCK).max(axis=(1, 3))

    def draw_null(self, rng) -> None:
        """A fresh 200-draw null of the real arcs (``exotic_screens.anti_window_draw``)."""
        n_rand = self.args.n_random
        bx, by = self.base_xy
        self.rand_pa = np.empty((n_rand, len(bx)))
        self.rand = np.zeros((n_rand, *self.gx.shape), np.int16)
        self.rand_block = []
        for d in range(n_rand):
            self.rand_pa[d] = es.anti_window_draw(self.pa_pred, rng)
            self._add(self.rand[d], bx, by, self.rand_pa[d], +1)
            self.rand_block.append(self._block_max(self.rand[d]))
        self.max_rand = np.array([b.max() for b in self.rand_block], int)

    def p_random(self, n_lines, max_rand=None) -> float:
        mr = self.screen_max_rand if max_rand is None else max_rand
        return float(np.mean(mr >= n_lines))

    def min_lines_for(self, p: float = P_RECOVER) -> int:
        """Smallest line count the real-field null puts below p_random ``p``."""
        k = int(self.screen_max_rand.max()) + 1
        while k > 1 and self.p_random(k - 1) < p:
            k -= 1
        return k

    def _null_max(self, rem, add_xy, pa_add_pred, rng) -> np.ndarray:
        """Per-draw grid maximum after removing arcs ``rem`` and adding new arcs, recomputing only
        the grid blocks their windows touch."""
        bx, by = self.base_xy
        ax, ay = add_xy
        wins = [self._window(x, y) for x, y in zip(bx[rem], by[rem], strict=True)]
        wins += [self._window(x, y) for x, y in zip(ax, ay, strict=True)]
        i0 = min(w[0] for w in wins) // BLOCK
        i1 = (max(w[1] for w in wins) - 1) // BLOCK + 1
        j0 = min(w[2] for w in wins) // BLOCK
        j1 = (max(w[3] for w in wins) - 1) // BLOCK + 1
        out = np.empty(self.args.n_random, int)
        lo = np.iinfo(np.int16).min
        for d in range(self.args.n_random):
            bm = self.rand_block[d].copy()
            bm[i0:i1, j0:j1] = lo
            sub = self.rand[d][i0 * BLOCK : i1 * BLOCK, j0 * BLOCK : j1 * BLOCK].copy()
            org = (i0 * BLOCK, j0 * BLOCK)
            self._add(sub, bx[rem], by[rem], self.rand_pa[d][rem], -1, org)
            pa_r = es.anti_window_draw(pa_add_pred, rng)
            self._add(sub, ax, ay, pa_r, +1, org)
            out[d] = max(int(bm.max()), int(sub.max()))
        return out

    def trial(self, x_l, y_l, theta_row, rng, psf_sigma, pixel_scale) -> dict:
        """Inject one lens, run the screen, and report whether it recovers the centre."""
        keep, img, info = paint_lens(
            self.shapes,
            self.xs,
            self.ys,
            self.background,
            theta_row,
            self.model,
            x_l,
            y_l,
            psf_sigma,
            self.next_label,
            pixel_scale,
        )
        inj = self.shapes[keep]
        if len(img):
            inj = vstack([inj, img], metadata_conflicts="silent")
        inj.meta.update(
            provenance=schema.Provenance.SIMULATED.value,
            source=f"{self.shapes.meta.get('source', 'catalogue')} with one simulated W1 lens",
        )
        cand, _ = es.radial_candidates(self.model, inj, self.args, self.extra)
        labels = {int(v) for v in cand["label"]}
        removed = np.array([i for lab, i in self.label_index.items() if lab not in labels], int)
        added = cand[np.array([int(v) not in self.label_index for v in cand["label"]], bool)]
        counts = self.counts.copy()
        bx, by = self.base_xy
        pa_obs = np.asarray(self.base["pa_obs"], float)
        self._add(counts, bx[removed], by[removed], pa_obs[removed], -1)
        ax, ay = self.model.to_frame(added["ra"], added["dec"])
        ax, ay = np.atleast_1d(ax), np.atleast_1d(ay)
        self._add(counts, ax, ay, np.asarray(added["pa_obs"], float), +1)
        if len(removed) or len(added):
            max_rand = self._null_max(removed, (ax, ay), np.asarray(added["pa_pred_z2"]), rng)
        else:
            max_rand = self.max_rand
        peaks = es.convergence_peaks(counts, self.args.min_lines)
        best_n, best_p, best_sep = 0, 1.0, float("nan")
        for iy, ix in peaks:
            sep = float(np.hypot(self.gx[iy, ix] - x_l, self.gy[iy, ix] - y_l))
            if sep <= self.args.recover_tol and counts[iy, ix] > best_n:
                best_n, best_sep = int(counts[iy, ix]), sep
                best_p = self.p_random(best_n, max_rand)
        n_inj_cand = int(np.sum(np.asarray(added["label"]) >= self.next_label))
        return info | {
            "x": float(x_l),
            "y": float(y_l),
            "n_image_arcs": n_inj_cand,
            "n_lines": best_n,
            "p_random": best_p,
            "sep": best_sep,
            "recovered": bool(best_n and best_p < P_RECOVER),
        }


# ----------------------------------------------------------------------------- field set-up


def canucs_zout(path: Path) -> Table:
    """CANUCS DR1 photometry catalogue -> zout columns (docs/fields/macs1149.md): ``Z_SPEC`` where
    > 0 for all three, otherwise ``Z_ML``, ``Z160``, ``Z840``."""
    t = Table.read(path, hdu=1)
    zs = np.asarray(t["Z_SPEC"], float)
    spec = np.isfinite(zs) & (zs > 0)
    out = Table({"ra": np.asarray(t["RA"], float), "dec": np.asarray(t["DEC"], float)})
    for col, src in (("z_phot", "Z_ML"), ("z160", "Z160"), ("z840", "Z840")):
        out[col] = np.where(spec, zs, np.asarray(t[src], float))
    out.meta.update(provenance=schema.Provenance.DERIVED.value, source=str(path))
    return out


psf_sigma_px = es.psf_sigma_px


def pixel_scale_arcsec(cat_path: Path) -> float:
    """Catalogue pixel scale from the two most distant pixel/sky centroid pairs (``derived``)."""
    t = Table.read(cat_path)
    x, y = np.asarray(t["xcentroid"], float), np.asarray(t["ycentroid"], float)
    i, j = int(np.argmin(x + y)), int(np.argmax(x + y))
    sep = t["sky_centroid"][i].separation(t["sky_centroid"][j]).arcsec
    return float(sep / np.hypot(x[i] - x[j], y[i] - y[j]))


def load_field(name: str, out: Path):
    spec = FIELDS[name]
    root = paths.data_root()
    model, _, _ = lc.load_model(spec["model"])
    lc.apply_frame_offset(spec["model"], model)
    cat = root / spec["catalog"]
    shapes = lc.load_shapes(cat)
    photoz = None
    if spec["photoz"]:
        kind, rel = spec["photoz"]
        photoz = root / rel
        if kind == "canucs":
            zpath = out / f"{name}_canucs_zout.fits"
            if not zpath.exists():
                canucs_zout(photoz).write(zpath, overwrite=True)
            photoz = zpath
        lc.attach_photoz(shapes, photoz)
    extra = None
    if spec.get("spike_stars"):
        star_path = out / f"{name}_gaia_stars.ecsv"
        if not star_path.exists():
            import subprocess

            ra, dec = spec["centre"]
            cmd = [sys.executable, str(Path(__file__).parent / "gaia_stars.py"), str(ra), str(dec)]
            subprocess.run([*cmd, "--radius-arcmin", "4", "--out", str(star_path)], check=True)
        extra = es.read_spike_stars(star_path)
    return model, shapes, extra, cat, photoz


def run_field(name: str, args) -> dict:
    t0 = time.time()
    out = args.out / name
    out.mkdir(parents=True, exist_ok=True)
    model, shapes, extra, cat, photoz = load_field(name, out)
    sargs = SimpleNamespace(**SCREEN_DEFAULTS)
    sargs.max_radius = FIELDS[name]["max_radius"]
    sargs.recover_tol = args.recover_tol
    inj = RadialInjector(model, shapes, sargs, extra)
    t_base = time.time() - t0
    psf = psf_sigma_px(shapes)
    scale = pixel_scale_arcsec(cat)
    fx, fy, area = screened_footprint(inj.xs, inj.ys, sargs.max_radius)
    border = footprint_border(inj.xs, inj.ys, sargs.max_radius)
    bg_in = inj.background & (np.hypot(inj.xs, inj.ys) <= sargs.max_radius)
    sigma_bg = float(bg_in.sum() / area)  # lensable rows per arcsec^2 in the screened area
    z_src = source_redshift(shapes)
    seeds = np.random.SeedSequence(args.seed).spawn(len(args.mass))
    trials, eff = [], {}
    for mass, seed in zip(args.mass, seeds, strict=True):
        rng = np.random.default_rng(seed)
        theta_row = np.where(inj.background, theta_e_arcsec(mass, model.z_lens, z_src), np.nan)
        theta_2 = float(theta_e_arcsec(mass, model.z_lens, Z_SOURCE)[0])
        sub = []
        for k in range(args.n_inject):
            if k % BATCH == 0:
                inj.draw_null(rng)  # independent null realisation for each batch
            i = rng.integers(0, len(fx))
            jx, jy = rng.uniform(-FOOTPRINT_STEP / 2, FOOTPRINT_STEP / 2, 2)
            r = inj.trial(fx[i] + jx, fy[i] + jy, theta_row, rng, psf, scale)
            r["mass_msun"] = mass
            sub.append(r)
        trials.extend(sub)
        n_rec = int(sum(r["recovered"] for r in sub))
        p = n_rec / args.n_inject
        lensable = np.isfinite(theta_row) & bg_in
        eff[f"{mass:.0e}"] = {
            "mass_msun": mass,
            "theta_e_zs2_arcsec": theta_2,
            "theta_e_median_arcsec": float(np.nanmedian(theta_row[lensable])),
            "n_inject": args.n_inject,
            "n_recovered": n_rec,
            "efficiency": p,
            "efficiency_err": float(np.sqrt(p * (1 - p) / args.n_inject)),
            "mean_lensed": float(np.mean([r["n_lensed"] for r in sub])),
            "mean_merged": float(np.mean([r["n_merged"] for r in sub])),
            "mean_image_arcs": float(np.mean([r["n_image_arcs"] for r in sub])),
            "frac_ge_min_lines": float(np.mean([r["n_lines"] >= sargs.min_lines for r in sub])),
        }
    t_all = time.time() - t0
    table = Table(rows=trials)
    table.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source=f"inject_radial.py {name}: simulated W1 lenses in {cat}",
    )
    table.write(out / "trials.ecsv", overwrite=True)
    summary = {
        "field": name,
        "model": FIELDS[name]["model"],
        "z_lens": model.z_lens,
        "catalog": str(cat),
        "photoz": str(photoz) if photoz else None,
        "spike_stars": len(extra) if extra is not None else None,
        "screen": inj.base_counts_info
        | {
            "max_lines": int(inj.counts.max()),
            "p_random_max": inj.p_random(int(inj.counts.max())),
            "min_lines_p05": inj.min_lines_for(),
        },
        "screened_area_arcsec2": area,
        "screened_area_deg2": area / 3600.0**2,
        "footprint_border": border,
        "lensable_density_arcsec2": sigma_bg,
        "psf_sigma_px": psf,
        "pixel_scale_arcsec": scale,
        "efficiency": eff,
        "wall_time_s": {"base_screen": round(t_base, 1), "total": round(t_all, 1)},
        "assumptions": {
            "beta_range_theta_e": [BETA_MIN, BETA_MAX],
            "z_source_without_photoz": Z_SOURCE,
            "snr_detect": SNR_DETECT,
            "recover_tol_arcsec": args.recover_tol,
            "p_recover": P_RECOVER,
            "trials_per_null": BATCH,
            "footprint": {"step": FOOTPRINT_STEP, "radius": FOOTPRINT_RADIUS},
            "screen": SCREEN_DEFAULTS | {"max_radius": sargs.max_radius},
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    return summary


def _limits(summaries: list[dict], masses) -> dict:
    limits = {}
    for m in masses:
        key = f"{m:.0e}"
        effs = [s["efficiency"][key]["efficiency"] for s in summaries]
        areas = [s["screened_area_deg2"] for s in summaries]
        th = [s["efficiency"][key]["theta_e_zs2_arcsec"] for s in summaries]
        # border-corrected: each footprint shrunk to its r -> 0 area (conservative)
        a0 = [
            a * (1.0 - s.get("footprint_border", {}).get("excess_frac", 0.0))
            for a, s in zip(areas, summaries, strict=True)
        ]
        limits[key] = {
            "mass_msun": m,
            "theta_e_zs2_arcsec_range": [min(th), max(th)],
            "effective_area_deg2": float(np.dot(effs, areas)),
            "upper_limit_deg2": surface_density_limit(effs, areas),
            "upper_limit_border_corrected_deg2": surface_density_limit(effs, a0),
        }
    return limits


def combine(summaries: list[dict], masses) -> dict:
    """Headline limits use only fields with photo-z. Without photo-z every non-star row counts as
    lensable (at z_s = 2), so cluster members and foreground galaxies get painted as W1 images and
    the efficiency is biased high; those fields enter only the separate, optimistic set."""
    pz = [s for s in summaries if s.get("photoz")]  # as the field was actually run
    for s in summaries:
        if bool(s.get("photoz")) != (FIELDS[s["field"]]["photoz"] is not None):
            print(
                f"warning: {s['field']}: summary photo-z differs from FIELDS; re-run it", flush=True
            )
    out = {
        "fields": [s["field"] for s in pz],
        "total_area_deg2": float(sum(s["screened_area_deg2"] for s in pz)),
        "limits_95": _limits(pz, masses) if pz else {},
        "provenance": schema.Provenance.DERIVED.value,
    }
    if len(pz) < len(summaries):
        out["optimistic_all_fields"] = {
            "fields": [s["field"] for s in summaries],
            "no_photoz": [s["field"] for s in summaries if s not in pz],
            "total_area_deg2": float(sum(s["screened_area_deg2"] for s in summaries)),
            "limits_95": _limits(summaries, masses),
        }
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--fields", nargs="+", default=sorted(FIELDS), choices=sorted(FIELDS))
    ap.add_argument("--mass", nargs="+", type=float, default=list(MASSES), help="|M| in M_sun")
    ap.add_argument("--n-inject", type=int, default=200, help="lenses per field and mass")
    ap.add_argument("--recover-tol", type=float, default=RECOVER_TOL)
    ap.add_argument("--seed", type=int, default=49)
    ap.add_argument("--out", type=Path, default=paths.outputs_dir() / "inject_radial")
    ap.add_argument(
        "--combine-only", action="store_true", help="combine existing <out>/<field>/summary.json"
    )
    args = ap.parse_args(argv)
    keys = [f"{m:.0e}" for m in args.mass]  # result keys (as stored in summary.json)
    if len(set(keys)) != len(keys):
        raise SystemExit(f"error: --mass values {args.mass} collide in the result keys {keys}")
    summaries = []
    for name in args.fields:
        if args.combine_only:
            summaries.append(json.loads((args.out / name / "summary.json").read_text()))
            continue
        s = run_field(name, args)
        keys = ("field", "efficiency", "screened_area_deg2", "wall_time_s", "screen")
        print(json.dumps({k: s[k] for k in keys}, indent=1), flush=True)
        summaries.append(s)
    combined = combine(summaries, args.mass)
    (args.out / "limits.json").write_text(json.dumps(combined, indent=1))
    print(json.dumps(combined, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
