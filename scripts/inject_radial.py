"""Injection-recovery of W1 negative-mass lenses through the ``radial`` screen (D-049).

W1 (docs/exotic_lensing.md, D-047) is a compact lens with ε < 0 and n = 1: no images of sources
within 2 θ_E (the umbra), and two radially stretched images on the source's side of every source
beyond. This script paints such lenses into a field's real catalogue, runs the screen's own arc
selection and convergence search (``exotic_screens.radial_candidates``, ``line_counts``,
``convergence_peaks``, the same thresholds and null), and measures how often a lens is recovered.
With the real screens null, the recovery efficiencies give 95 % upper limits on the surface density
of such lenses.

Everything injected is ``simulated``; efficiencies and limits are ``derived`` from simulated signals
in real catalogues. A limit says how many such lenses the screen would have seen; it is not a
detection of anything, and the masses are a ``model_prediction`` of the point-lens formula.

Injection (all ASSUMPTIONs, documented in docs/exotic_limits.md):
- lens at a uniform random position in the screened footprint (catalogue sources within 4" of a
  1" grid point, inside ``--max-radius`` of the model centre), at the cluster redshift;
- lensed sources are the catalogue's own rows behind the lens position: rows with a background or
  no photo-z (the screen's rule), not bright point sources (stars). Rows with beta < 2 theta_E
  (umbra) are removed; rows with 2 <= beta <= 4 theta_E are replaced by their two images
  (``exotic_sim.inject_images``); rows beyond 4 theta_E are left as they are (D-047 range);
- each image keeps its source's columns; its shape is the source's second-moment ellipse, PSF
  deconvolved, mapped by the lens Jacobian (signed eigenvalues along the radial and tangential
  directions) and re-convolved; magnitude - 2.5 log10 |mu|, isophotal area x |mu| (surface
  brightness is conserved) and S/N x sqrt(|mu|) (background-limited). Images below the catalogue's
  own S/N floor are dropped (not detected);
- the PSF sigma is the 1st percentile of ``semiminor_sigma`` among S/N > 50 catalogue sources
  (``derived``).

Recovery: a convergence peak with p_random < 0.05 (the screen's definition) within
``--recover-tol`` (2") of the injected centre. The null is the screen's: 200 draws of each arc's
position angle inside its ``anti`` window. The real arcs' draws are computed once per field with
the screen's seed; injected or removed arcs are added to or subtracted from each draw. The result
has the same distribution as a full re-run, at a fraction of the cost.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from astropy.table import Table, vstack

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exotic_screens as es  # noqa: E402
import lens_consistency as lc  # noqa: E402

from jwst_anomaly import exotic_sim, paths, schema  # noqa: E402

BETA_MIN, BETA_MAX = 2.0, 4.0  # lensed-source range in theta_E (D-047 recommendation; ASSUMPTION)
P_RECOVER = 0.05  # the screen's p_random threshold
RECOVER_TOL = 2.0  # arcsec between the peak and the injected centre (ASSUMPTION)
FOOTPRINT_STEP = 1.0  # arcsec grid for the screened footprint (ASSUMPTION)
FOOTPRINT_RADIUS = (
    4.0  # a grid point is screened if a catalogue source lies within this (ASSUMPTION)
)
POISSON_UL_95 = 2.996  # -ln(0.05): 95 % upper limit on a Poisson mean with zero events
Z_SOURCE = 2.0  # source redshift for the theta_E -> |M| conversion (ASSUMPTION)
SCREEN_DEFAULTS = es.radial_defaults()  # the settings the field screens ran with

# The inputs each field's null screen used (docs/fields/*.md, docs/exotic_lensing.md). Paths are
# relative to the data root. Photo-z kinds: "eazy" (DJA zout) or "canucs" (CANUCS DR1 catalogue,
# converted as in docs/fields/macs1149.md).
_MAST = "cache/mast/{0}/{0}_cat.ecsv"
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
        "photoz": (
            "canucs",
            "cache/external/canucs/hlsp_canucs_jwst-hst_multi_macs0416-clu_multi_v1_photometry-cat.fits.gz",
        ),
        "max_radius": 100.0,
    },
    "macs1149": {
        "model": "macs1149-cats",
        "catalog": _MAST.format("jw01208-o008_t004_nircam_clear-f200w"),
        "photoz": (
            "canucs",
            "cache/external/canucs/hlsp_canucs_jwst-hst_multi_macs1149-clu_multi_v1_photometry-cat.fits.gz",
        ),
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
        "photoz": (
            "canucs",
            "cache/external/canucs/hlsp_canucs_jwst-hst_multi_a370-clu_multi_v1_photometry-cat.fits.gz",
        ),
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


def theta_e_to_mass(theta_e_arcsec, z_l: float, z_s: float = Z_SOURCE, cosmology=None):
    """|M| in M_sun of an n = 1 lens with Einstein (scale) radius ``theta_e_arcsec``.

    Inverts ``exotic_sim.einstein_radius`` with ``eps_bar = 4 G |M| / c^2`` (``model_prediction``;
    Planck18 unless ``cosmology`` is given)."""
    if cosmology is None:
        from astropy.cosmology import Planck18 as cosmology
    d_l = cosmology.angular_diameter_distance(z_l).si.value
    d_s = cosmology.angular_diameter_distance(z_s).si.value
    d_ls = cosmology.angular_diameter_distance(z_l, z_s).si.value
    theta = np.deg2rad(np.asarray(theta_e_arcsec, float) / 3600.0)
    eps_bar = theta**2 * d_s * d_l / d_ls
    return eps_bar * exotic_sim.C_SI**2 / (4.0 * exotic_sim.G_SI * exotic_sim.MSUN_KG)


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


def lensed_shapes(a, b, pa_deg, dx, dy, lam_r, lam_t, psf_sigma: float):
    """Image shapes of sources (sigmas ``a``, ``b`` in px, ``pa_deg``) seen through the lens.

    ``dx``, ``dy``: image offsets from the lens (East, North); ``lam_r``, ``lam_t``: signed Jacobian
    eigenvalues at the image (``exotic_sim.solve_images``). The PSF (isotropic, ``psf_sigma`` px) is
    removed from the source moments, the intrinsic ellipse mapped by A^-1, and the PSF added back
    (ASSUMPTION: Gaussian moments, isolated lens). Returns (semimajor, ellipticity, pa_deg)."""
    cov = _ellipse_cov(a, b, pa_deg)
    eye = np.eye(2)
    floor = (0.05 * psf_sigma) ** 2
    w, vec = np.linalg.eigh(cov - psf_sigma**2 * eye)
    cov_int = vec @ (np.maximum(w, floor)[..., None] * np.swapaxes(vec, -1, -2))
    r = np.stack([dx, dy], axis=-1)
    r = r / np.linalg.norm(r, axis=-1, keepdims=True)
    t = np.stack([-r[:, 1], r[:, 0]], axis=-1)
    inv = (1.0 / np.asarray(lam_r, float))[:, None, None] * r[:, :, None] * r[:, None, :] + (
        1.0 / np.asarray(lam_t, float)
    )[:, None, None] * t[:, :, None] * t[:, None, :]
    cov_img = inv @ cov_int @ np.swapaxes(inv, -1, -2) + psf_sigma**2 * eye
    w, vec = np.linalg.eigh(cov_img)
    major = vec[:, :, 1]  # eigh sorts ascending
    pa = np.mod(np.degrees(np.arctan2(major[:, 0], major[:, 1])), 180.0)
    return np.sqrt(w[:, 1]), 1.0 - np.sqrt(w[:, 0] / w[:, 1]), pa


def background_mask(shapes: Table, z_lens: float, z_margin: float = 0.1) -> np.ndarray:
    """Rows a lens at ``z_lens`` can lens: the screen's background rule (no photo-z, or one that
    puts the source behind; ``orientation_table``'s z160 > z_lens + margin), minus bright point
    sources (``exotic_screens.SPIKE_STAR_MAG``; stars)."""
    mag = np.asarray(shapes["mag"], float)
    star = ~np.asarray(shapes["is_extended"], bool) & np.isfinite(mag) & (mag < es.SPIKE_STAR_MAG)
    ok = ~star
    if "z_phot" in shapes.colnames:
        has_pz = np.isfinite(np.asarray(shapes["z_phot"], float))
        z160 = np.asarray(shapes["z160"], float)
        behind = np.isfinite(z160) & (z160 > z_lens + z_margin)
        ok &= ~has_pz | behind
    return ok


def paint_lens(
    shapes: Table,
    xs: np.ndarray,
    ys: np.ndarray,
    background: np.ndarray,
    model,
    x_l: float,
    y_l: float,
    theta_e: float,
    psf_sigma: float,
    snr_floor: float,
    first_label: int,
) -> tuple[np.ndarray, Table, dict]:
    """Paint one W1 lens at model-frame (``x_l``, ``y_l``) into ``shapes`` (``simulated`` rows).

    ``xs``, ``ys``: model-frame positions of ``shapes`` (x West, y North). Returns the mask of rows
    to keep, the image rows (same columns, new labels from ``first_label``) and counts."""
    dx = -(xs - x_l)  # East
    dy = ys - y_l
    beta = np.hypot(dx, dy) / theta_e
    hit = background & (beta < BETA_MAX)
    lensed = np.flatnonzero(hit & (beta >= BETA_MIN))
    keep = ~hit
    info = {"n_umbra": int(np.sum(hit & (beta < BETA_MIN))), "n_lensed": len(lensed)}
    if not len(lensed):
        info["n_images"] = 0
        return keep, shapes[:0].copy(), info
    img = exotic_sim.inject_images(
        dx[lensed], dy[lensed], theta_e, n=1.0, sign=-1, source_id=lensed
    )
    rows = shapes[np.asarray(img["source_id"], int)]
    k = np.asarray(img["source_id"], int)
    lam_r = 1.0 / np.asarray(img["stretch_r"], float)
    lam_t = 1.0 / np.asarray(img["stretch_t"], float)
    # signs: lam_t > 0 always for eps < 0; lam_r < 0 for the inner image (x < 1)
    lam_r = np.where(np.asarray(img["image"]) == 1, -lam_r, lam_r)
    a = np.asarray(shapes["semimajor_px"], float)[k]
    b = a * (1.0 - np.asarray(shapes["ellipticity"], float)[k])
    dxi, dyi = np.asarray(img["dx"], float), np.asarray(img["dy"], float)
    pa_src = np.asarray(shapes["pa_obs"], float)[k]
    semi, ell, pa = lensed_shapes(a, b, pa_src, dxi, dyi, lam_r, lam_t, psf_sigma)
    mu = np.abs(np.asarray(img["mu"], float))
    rows["semimajor_px"] = semi
    rows["ellipticity"] = ell
    rows["pa_obs"] = pa
    rows["mag"] = np.asarray(rows["mag"], float) - 2.5 * np.log10(mu)
    rows["area_px"] = np.asarray(rows["area_px"], float) * mu
    rows["snr"] = np.asarray(rows["snr"], float) * np.sqrt(mu)
    ra, dec = model.to_sky(x_l - np.asarray(img["dx"], float), y_l + np.asarray(img["dy"], float))
    rows["ra"], rows["dec"] = ra, dec
    rows["label"] = first_label + np.arange(len(rows))
    rows = rows[np.asarray(rows["snr"], float) >= snr_floor]
    rows.meta = {
        "provenance": schema.Provenance.SIMULATED.value,
        "source": f"exotic_sim.inject_images W1 (n=1, eps<0) theta_E={theta_e} at ({x_l}, {y_l})",
    }
    info["n_images"] = len(rows)
    return keep, rows, info


def screened_footprint(xs, ys, max_radius: float, step: float = FOOTPRINT_STEP):
    """Model-frame grid points (``step`` arcsec) inside ``max_radius`` with a catalogue source
    within ``FOOTPRINT_RADIUS``: the area the screen looked at. Returns (gx, gy, area_arcsec2)."""
    from scipy.spatial import cKDTree

    g = np.arange(-max_radius + step / 2, max_radius, step)
    gx, gy = np.meshgrid(g, g)
    gx, gy = gx.ravel(), gy.ravel()
    inside = np.hypot(gx, gy) <= max_radius
    gx, gy = gx[inside], gy[inside]
    d, _ = cKDTree(np.column_stack([xs, ys])).query(np.column_stack([gx, gy]))
    ok = d <= FOOTPRINT_RADIUS
    return gx[ok], gy[ok], float(ok.sum() * step**2)


# ----------------------------------------------------------------------------- the screen


class RadialInjector:
    """The radial screen on one field, with cached real-arc lines and null draws."""

    def __init__(self, model, shapes: Table, args, extra: Table | None = None):
        self.model, self.shapes, self.args, self.extra = model, shapes, args, extra
        cand, self.base_counts_info = es.radial_candidates(model, shapes, args, extra)
        self.base = cand
        g = np.arange(-args.max_radius, args.max_radius + 1e-9, args.grid_arcsec)
        self.g = g
        self.gx, self.gy = np.meshgrid(g, g)
        cx, cy = model.to_frame(cand["ra"], cand["dec"])
        self.base_xy = (np.asarray(cx), np.asarray(cy))
        self.counts = es.line_counts(
            cx, cy, cand["pa_obs"], self.gx, self.gy, args.line_tol_arcsec, args.max_len_arcsec
        ).astype(np.int16)
        # the screen's null, draw by draw, with its seed and draw order
        rng = np.random.default_rng(args.seed)
        pa_t = np.asarray(cand["pa_pred_z2"], float)
        self.rand_pa = np.empty((args.n_random, len(cand)))
        self.rand = np.empty((args.n_random, *self.gx.shape), np.int16)
        for d in range(args.n_random):
            pa_r = np.mod(pa_t + 90.0 + rng.uniform(-30.0, 30.0, len(cx)), 180.0)
            self.rand_pa[d] = pa_r
            self.rand[d] = es.line_counts(
                cx, cy, pa_r, self.gx, self.gy, args.line_tol_arcsec, args.max_len_arcsec
            )
        self.max_rand = self.rand.reshape(args.n_random, -1).max(axis=1)
        self.label_index = {int(v): i for i, v in enumerate(cand["label"])}
        xs, ys = model.to_frame(shapes["ra"], shapes["dec"])
        self.xs, self.ys = np.asarray(xs), np.asarray(ys)
        self.background = background_mask(shapes, model.z_lens)
        self.next_label = int(np.max(shapes["label"])) + 1

    def p_random(self, n_lines, max_rand=None) -> float:
        mr = self.max_rand if max_rand is None else max_rand
        return float(np.mean(mr >= n_lines))

    def min_lines_for(self, p: float = P_RECOVER) -> int:
        """Smallest line count the real-field null puts below p_random ``p``."""
        k = int(self.max_rand.max()) + 1
        while k > 1 and self.p_random(k - 1) < p:
            k -= 1
        return k

    def _window(self, x, y):
        r = self.args.max_len_arcsec + self.args.line_tol_arcsec
        step, g0 = self.args.grid_arcsec, self.g[0]
        i0 = max(int(np.floor((y - r - g0) / step)), 0)
        i1 = min(int(np.ceil((y + r - g0) / step)) + 1, len(self.g))
        j0 = max(int(np.floor((x - r - g0) / step)), 0)
        j1 = min(int(np.ceil((x + r - g0) / step)) + 1, len(self.g))
        return slice(i0, i1), slice(j0, j1)

    def _add(self, grid, x, y, pa, sign: int) -> None:
        a = self.args
        for xi, yi, pi in zip(x, y, pa, strict=True):
            sy, sx = self._window(xi, yi)
            grid[sy, sx] += sign * es.line_counts(
                [xi],
                [yi],
                [pi],
                self.gx[sy, sx],
                self.gy[sy, sx],
                a.line_tol_arcsec,
                a.max_len_arcsec,
            ).astype(np.int16)

    def trial(self, x_l: float, y_l: float, theta_e: float, rng, psf_sigma, snr_floor) -> dict:
        """Inject one lens, run the screen, and report whether it recovers the centre."""
        keep, img, info = paint_lens(
            self.shapes,
            self.xs,
            self.ys,
            self.background,
            self.model,
            x_l,
            y_l,
            theta_e,
            psf_sigma,
            snr_floor,
            self.next_label,
        )
        inj = (
            vstack([self.shapes[keep], img], metadata_conflicts="silent")
            if len(img)
            else (self.shapes[keep])
        )
        cand, _ = es.radial_candidates(self.model, inj, self.args, self.extra)
        labels = {int(v) for v in cand["label"]}
        removed = [i for lab, i in self.label_index.items() if lab not in labels]
        added = cand[np.array([int(v) not in self.label_index for v in cand["label"]], bool)]
        counts = self.counts.copy()
        bx, by = self.base_xy
        pa_obs = np.asarray(self.base["pa_obs"], float)
        self._add(counts, bx[removed], by[removed], pa_obs[removed], -1)
        ax, ay = self.model.to_frame(added["ra"], added["dec"])
        ax, ay = np.atleast_1d(ax), np.atleast_1d(ay)
        self._add(counts, ax, ay, np.asarray(added["pa_obs"], float), +1)
        # null: the cached draws of the real arcs, minus removed arcs, plus fresh draws of new ones
        changed = bool(removed) or len(added) > 0
        if changed:
            max_rand = np.empty(self.args.n_random, int)
            pa_t = np.asarray(added["pa_pred_z2"], float)
            for d in range(self.args.n_random):
                rc = self.rand[d].copy()
                self._add(rc, bx[removed], by[removed], self.rand_pa[d][removed], -1)
                pa_r = np.mod(pa_t + 90.0 + rng.uniform(-30.0, 30.0, len(added)), 180.0)
                self._add(rc, ax, ay, pa_r, +1)
                max_rand[d] = rc.max()
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
    return out


def psf_sigma_px(cat_path: Path) -> float:
    """1st percentile of ``semiminor_sigma`` over S/N > 50 sources: the narrowest objects are
    PSF-limited in their minor axis (``derived``)."""
    t = Table.read(cat_path)
    with np.errstate(divide="ignore", invalid="ignore"):
        snr = np.asarray(t["isophotal_flux"], float) / np.asarray(t["isophotal_flux_err"], float)
    b = np.asarray(t["semiminor_sigma"], float)
    sel = np.isfinite(snr) & (snr > 50) & np.isfinite(b)
    if not sel.any():
        raise ValueError(f"{cat_path}: no source with S/N > 50 to estimate the PSF width")
    return float(np.percentile(b[sel], 1))


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
            subprocess.run(
                [
                    sys.executable,
                    str(Path(__file__).parent / "gaia_stars.py"),
                    str(ra),
                    str(dec),
                    "--radius-arcmin",
                    "4",
                    "--out",
                    str(star_path),
                ],
                check=True,
            )
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
    psf = psf_sigma_px(cat)
    snr_floor = float(np.nanmin(np.asarray(shapes["snr"], float)))
    fx, fy, area = screened_footprint(inj.xs, inj.ys, sargs.max_radius)
    bg_in = inj.background & (np.hypot(inj.xs, inj.ys) <= sargs.max_radius)
    sigma_bg = float(bg_in.sum() / area)  # lensable rows per arcsec^2 in the screened area
    rng = np.random.default_rng(args.seed)
    trials, eff = [], {}
    for theta_e in args.theta_e:
        pick = rng.integers(0, len(fx), args.n_inject)
        jitter = rng.uniform(-FOOTPRINT_STEP / 2, FOOTPRINT_STEP / 2, (args.n_inject, 2))
        rec = []
        for k in range(args.n_inject):
            r = inj.trial(
                fx[pick[k]] + jitter[k, 0], fy[pick[k]] + jitter[k, 1], theta_e, rng, psf, snr_floor
            )
            r["theta_e"] = theta_e
            trials.append(r)
            rec.append(r["recovered"])
        n_rec = int(np.sum(rec))
        sub = [r for r in trials if r["theta_e"] == theta_e]
        eff[str(theta_e)] = {
            "n_inject": args.n_inject,
            "n_recovered": n_rec,
            "efficiency": n_rec / args.n_inject,
            "efficiency_err": float(
                np.sqrt(n_rec / args.n_inject * (1 - n_rec / args.n_inject) / args.n_inject)
            ),
            "expected_lensed": sigma_bg * np.pi * (BETA_MAX**2 - BETA_MIN**2) * theta_e**2,
            "mean_lensed": float(np.mean([r["n_lensed"] for r in sub])),
            "mean_image_arcs": float(np.mean([r["n_image_arcs"] for r in sub])),
            "frac_ge_min_lines": float(np.mean([r["n_lines"] >= sargs.min_lines for r in sub])),
            "mass_msun": float(theta_e_to_mass(theta_e, model.z_lens)),
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
        "lensable_density_arcsec2": sigma_bg,
        "psf_sigma_px": psf,
        "snr_floor": snr_floor,
        "efficiency": eff,
        "wall_time_s": {"base_screen": round(t_base, 1), "total": round(t_all, 1)},
        "assumptions": {
            "beta_range_theta_e": [BETA_MIN, BETA_MAX],
            "recover_tol_arcsec": args.recover_tol,
            "p_recover": P_RECOVER,
            "footprint": {"step": FOOTPRINT_STEP, "radius": FOOTPRINT_RADIUS},
            "z_source_for_mass": Z_SOURCE,
            "screen": SCREEN_DEFAULTS | {"max_radius": sargs.max_radius},
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    return summary


def _limits(summaries: list[dict]) -> dict:
    """95 % limits per theta_E present in every summary (keys as stored)."""
    keys = set.intersection(*(set(s["efficiency"]) for s in summaries))
    limits = {}
    for te in sorted(keys, key=float):
        effs = [s["efficiency"][te]["efficiency"] for s in summaries]
        areas = [s["screened_area_deg2"] for s in summaries]
        masses = [s["efficiency"][te]["mass_msun"] for s in summaries]
        limits[te] = {
            "effective_area_deg2": float(np.dot(effs, areas)),
            "upper_limit_deg2": surface_density_limit(effs, areas),
            "mass_msun_range": [min(masses), max(masses)],
        }
    return limits


def combine(summaries: list[dict]) -> dict:
    """Headline limits use only fields with photo-z. Without photo-z every non-star row counts as
    lensable, so cluster members and foreground galaxies get painted as W1 images and the
    efficiency is biased high; those fields are reported separately as optimistic."""
    pz = [s for s in summaries if s.get("photoz")]
    out = {
        "fields": [s["field"] for s in pz],
        "total_area_deg2": float(sum(s["screened_area_deg2"] for s in pz)),
        "limits_95": _limits(pz) if pz else {},
        "provenance": schema.Provenance.DERIVED.value,
    }
    if len(pz) < len(summaries):
        out["optimistic_all_fields"] = {
            "fields": [s["field"] for s in summaries],
            "no_photoz": [s["field"] for s in summaries if not s.get("photoz")],
            "total_area_deg2": float(sum(s["screened_area_deg2"] for s in summaries)),
            "limits_95": _limits(summaries),
        }
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument(
        "--fields", nargs="+", default=["smacs0723", "elgordo", "abell2744"], choices=sorted(FIELDS)
    )
    ap.add_argument("--theta-e", nargs="+", type=float, default=[0.3, 1.0, 3.0])
    ap.add_argument("--n-inject", type=int, default=100, help="lenses per field and theta_E")
    ap.add_argument("--recover-tol", type=float, default=RECOVER_TOL)
    ap.add_argument("--seed", type=int, default=49)
    ap.add_argument("--out", type=Path, default=paths.outputs_dir() / "inject_radial")
    ap.add_argument(
        "--combine-only", action="store_true", help="combine existing <out>/<field>/summary.json"
    )
    args = ap.parse_args(argv)
    summaries = []
    for name in args.fields:
        if args.combine_only:
            summaries.append(json.loads((args.out / name / "summary.json").read_text()))
            continue
        s = run_field(name, args)
        print(
            json.dumps(
                {
                    k: s[k]
                    for k in ("field", "efficiency", "screened_area_deg2", "wall_time_s", "screen")
                },
                indent=1,
            ),
            flush=True,
        )
        summaries.append(s)
    combined = combine(summaries)
    (args.out / "limits.json").write_text(json.dumps(combined, indent=1))
    print(json.dumps(combined, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
