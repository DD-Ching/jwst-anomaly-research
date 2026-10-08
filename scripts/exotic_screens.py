"""Exotic-lens screens on a lensing-cluster field (docs/exotic_lensing.md; D-031).

Two model-assisted screens for the patterns the literature predicts for exotic lenses. Every hit is
a *candidate for vetting* (``/vet-candidate``), never evidence: each screen names the ordinary
effects that produce the same pattern, and tests the cheap ones itself.

``fluxratio``  (demagnification; Ellis-wormhole-like lenses, Kitamura et al. 2013)
    Forced photometry of every catalogued multiple image in two bands. Each image's flux divided by
    its model |μ| estimates the unlensed source flux; an image whose estimate is below 1/3 of the
    median of its siblings is *underluminous* (demagnified relative to the model), above 3x
    *overluminous*. Only compact images (f(0.2")/f(0.4") >= 0.6) take part: a resolved arc's
    aperture flux follows its surface brightness, which lensing conserves, not |μ|. Ordinary
    explanations tested here: |μ| > 50 (μ unreliable near critical curves),
    chromatic ratios (dust or a blend: the two bands disagree by more than 0.5 mag). Left for
    vetting: microlensing by intracluster stars, substructure, time-delayed variability.

``radial``  (negative convergence: radially stretched images around a dark centre; Izumi et al.
2013)
    Elongated background sources whose major axis is perpendicular to the model's predicted stretch
    (``anti`` in ``lens_consistency.py arcs``) and that the model does not predict to be radial
    (radial magnification 1/|1 - κ + γ| < 3 at their redshift). A grid search then finds points
    where at least ``--min-lines`` such major axes pass within ``--line-tol-arcsec``; a point with
    no catalog source within ``--dark-radius-arcsec`` is a *dark-centre candidate*. The false-alarm
    rate comes from the same search with each arc's position angle redrawn inside its own anti
    window (``--n-random`` draws).
    Ordinary explanations left for vetting: blends, intrinsic alignments, voids, model error.

``shear``  (negative tangential shear; W1, D-047/D-050)
    A catalogue aperture-mass map (Schneider 1996): on a grid of centres, the S/N of the weighted
    mean *radial* alignment of background shapes (PSF-deconvolved second moments, the cluster
    model's reduced shear removed, rows with |g| >= 0.5 or κ >= 1 masked). The null rotates every
    shape by a random angle (positions kept); a centre's ``p_random`` is the fraction of draws
    whose field maximum reaches its S/N. Ordinary explanations left for vetting: cluster-model
    shear errors, the cluster's own radial arcs, PSF anisotropy, blends, intrinsic alignments.

Every threshold is an ASSUMPTION (defaults below); outputs are ``derived`` with ``model_prediction``
inputs. Results go to ``outputs/exotic_screens/<model>/``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lens_consistency as lc  # noqa: E402

from jwst_anomaly import lensmodel, paths, schema  # noqa: E402

MAX_REF_MU = 50.0  # |μ| above this is too uncertain to scale fluxes by (ASSUMPTION)
LOW, HIGH = 1 / 3, 3.0  # luminosity-ratio bounds for under- and overluminous images (ASSUMPTION)
MAX_CHROMATIC_MAG = 0.5  # band-to-band disagreement that marks a ratio as chromatic (ASSUMPTION)
MIN_CONCENTRATION = (
    0.6  # f(0.2")/f(0.4") below this: resolved, aperture flux not ∝ |μ| (ASSUMPTION)
)
MAX_RADIAL_MU = 3.0  # predicted radial magnification above which a radial arc is ordinary


# --------------------------------------------------------------------------------- fluxratio


def luminosity_ratios(flux: np.ndarray, err: np.ndarray, mu: np.ndarray, systems) -> Table:
    """Per image: ``L = flux / |μ|`` relative to the median ``L`` of its siblings.

    Only images with S/N > 5 and |μ| <= ``MAX_REF_MU`` take part; others get NaN. ``ratio_snr`` is
    the significance of ``ratio`` differing from 1 (flux errors only)."""
    flux, err, mu = (np.asarray(a, float) for a in (flux, err, mu))
    systems = np.asarray(systems).astype(str)
    ok = np.isfinite(flux) & (err > 0) & (flux / np.where(err > 0, err, np.inf) > 5)
    ok &= np.isfinite(mu) & (np.abs(mu) <= MAX_REF_MU)
    lum = np.where(ok, flux / np.abs(mu), np.nan)
    ratio = np.full(len(flux), np.nan)
    rsnr = np.full(len(flux), np.nan)
    for s in dict.fromkeys(systems):
        idx = np.flatnonzero((systems == s) & ok)
        if len(idx) < 2:
            continue
        for k in idx:
            o = idx[idx != k]
            ref = float(np.median(lum[o]))
            # error of the reference: median-of-n ~ 1.25 x mean error / sqrt(n) (ASSUMPTION)
            ref_err = 1.25 * float(np.mean(err[o] / np.abs(mu[o]))) / np.sqrt(len(o))
            ratio[k] = lum[k] / ref
            total = np.hypot(err[k], ref_err * abs(mu[k]))
            rsnr[k] = abs(flux[k] - ref * abs(mu[k])) / total
    out = Table({"system": systems, "lum_ratio": ratio, "ratio_snr": rsnr})
    out["usable"] = ok
    return out


def peak_offset(img, xx, yy, radius: float = 0.4) -> tuple[float, float]:
    """Offset (arcsec) of the brightest finite pixel within ``radius`` (as ``lc.peak_flux``)."""
    inner = np.where((np.hypot(xx, yy) <= radius) & np.isfinite(img), img, -np.inf)
    if not np.isfinite(inner).any():
        return float("nan"), float("nan")
    k = int(np.argmax(inner))
    return float(xx.flat[k]), float(yy.flat[k])


def photometry(img, err, xx, yy) -> tuple[float, float, float]:
    """Flux and error at r = 0.2" and the concentration f(0.2")/f(0.4"), at one peak and with the
    same background annulus (``lc.aperture_snr``; NaN on gaps, edges or masks)."""
    cx, cy = peak_offset(img, xx, yy)
    if not np.isfinite(cx):
        return float("nan"), float("nan"), float("nan")
    f_small, e_small = lc.aperture_snr(img, err, xx, yy, cx, cy)
    f_big, _ = lc.aperture_snr(img, err, xx, yy, cx, cy, r_ap=0.4)
    conc = f_small / f_big if np.isfinite(f_big) and f_big > 0 else float("nan")
    return f_small, e_small, conc


def flux_class(r1: float, r2: float, s1: float, s2: float, compact: bool = True) -> str:
    """Class of an image from its luminosity ratios in two bands and their significances.

    ``compact`` False (a resolved image) makes it untestable: aperture flux then follows surface
    brightness, which lensing conserves, not |μ|."""
    if not compact:
        return "resolved"
    if not (np.isfinite(r1) and np.isfinite(r2)):
        return "untestable"
    if abs(2.5 * np.log10(r1 / r2)) > MAX_CHROMATIC_MAG:
        return "chromatic"  # dust, a blend or a colour gradient, not a lensing (achromatic) effect
    if r1 < LOW and r2 < LOW and min(s1, s2) >= 5:
        return "underluminous"
    if r1 > HIGH and r2 > HIGH and min(s1, s2) >= 5:
        return "overluminous"
    return "consistent"


def cmd_fluxratio(args) -> dict:
    if not lc.has_image_list(args.model):
        raise SystemExit(
            f"error: {args.model}: no usable multiple-image list (none published, or excluded by "
            "the image-plane rms gate of D-035)"
        )
    model, files, par = lc.load_model(args.model)
    images, zml = lc.image_list(args.model, files, par)
    lc.apply_frame_offset(args.model, model, images)  # into the JWST frame (D-034)
    bt = lensmodel.backtrace_images(model, images, zml)
    out = Table(
        {
            "image_id": bt["image_id"],
            "system": bt["system"],
            "ra": bt["ra"],
            "dec": bt["dec"],
            "mu": bt["magnification"],
        }
    )
    for tag, uri in (("b1", args.band1), ("b2", args.band2)):
        stamp, close = lc.image_stamper(uri, 1.5)
        try:
            stamps = [stamp(float(r["ra"]), float(r["dec"])) for r in out]
        finally:
            close()
        ph = [photometry(*st) for st in stamps]
        out[f"flux_{tag}"] = [p[0] for p in ph]
        out[f"err_{tag}"] = [p[1] for p in ph]
        out[f"conc_{tag}"] = [p[2] for p in ph]
    # only images unresolved in both bands take part: a resolved arc's aperture flux follows
    # its surface brightness, which lensing conserves, not |μ|
    compact = (np.asarray(out["conc_b1"], float) >= MIN_CONCENTRATION) & (
        np.asarray(out["conc_b2"], float) >= MIN_CONCENTRATION
    )
    out["compact"] = compact
    for tag in ("b1", "b2"):
        flux = np.where(compact, np.asarray(out[f"flux_{tag}"], float), np.nan)
        lr = luminosity_ratios(flux, out[f"err_{tag}"], out["mu"], out["system"])
        out[f"ratio_{tag}"], out[f"snr_{tag}"] = lr["lum_ratio"], lr["ratio_snr"]
    out["flux_class"] = [
        flux_class(r["ratio_b1"], r["ratio_b2"], r["snr_b1"], r["snr_b2"], bool(r["compact"]))
        for r in out
    ]
    out.meta.update(model._meta())
    out.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        bands=[args.band1, args.band2],
        assumptions={
            "max_ref_mu": MAX_REF_MU,
            "ratio_bounds": [LOW, HIGH],
            "max_chromatic_mag": MAX_CHROMATIC_MAG,
            "min_concentration": MIN_CONCENTRATION,
            "aperture": "lens_consistency.peak_flux (r = 0.2 arcsec, recentred within 0.4 arcsec)",
        },
    )
    dest = args.out / args.model
    lc._write(out, dest / "fluxratio.ecsv")
    classes = ("underluminous", "overluminous", "chromatic", "consistent", "resolved", "untestable")
    flagged = out[np.isin(out["flux_class"], ["underluminous", "overluminous"])]
    summary = {
        "model": args.model,
        "bands": [args.band1, args.band2],
        "classes": {c: int(np.sum(out["flux_class"] == c)) for c in classes},
        "flagged": [
            {
                "image": str(r["image_id"]),
                "class": str(r["flux_class"]),
                "mu": round(float(r["mu"]), 2),
                "ratio": [round(float(r["ratio_b1"]), 3), round(float(r["ratio_b2"]), 3)],
                "snr": [round(float(r["snr_b1"]), 1), round(float(r["snr_b2"]), 1)],
            }
            for r in flagged
        ],
    }
    (dest / "fluxratio.json").write_text(json.dumps(summary, indent=1))
    return summary


# ------------------------------------------------------------------------------------ radial


def radial_magnification(model, ra, dec, z) -> np.ndarray:
    """Predicted radial magnification 1/|1 - κ + γ| (the eigen-direction a radial arc stretches)."""
    p = model.evaluate(ra, dec, z)
    lam_r = 1.0 - np.asarray(p["kappa"], float) + np.asarray(p["gamma"], float)
    with np.errstate(divide="ignore"):
        return 1.0 / np.abs(lam_r)


def line_counts(x, y, pa_deg, gx, gy, tol: float, max_len: float) -> np.ndarray:
    """For each grid point, how many sources' major-axis lines pass within ``tol`` of it.

    A source at (x, y) (model frame, arcsec; x West, y North) with position angle ``pa_deg``
    (degrees E of N) defines a line; only points within ``max_len`` of the source along it
    count."""
    pa = np.deg2rad(np.asarray(pa_deg, float))
    ux, uy = -np.sin(pa), np.cos(pa)  # E of N in a West/North frame: East = -x
    counts = np.zeros(gx.shape, int)
    for xi, yi, ui, vi in zip(x, y, ux, uy, strict=True):
        dx, dy = gx - xi, gy - yi
        along = dx * ui + dy * vi
        perp = np.abs(dx * vi - dy * ui)
        counts += (perp <= tol) & (np.abs(along) <= max_len) & (np.abs(along) >= tol)
    return counts


def convergence_peaks(counts: np.ndarray, min_lines: int) -> list[tuple[int, int]]:
    """One grid index per connected region with at least ``min_lines`` lines: its maximum."""
    from scipy import ndimage

    labels, n = ndimage.label(counts >= min_lines)
    if n == 0:
        return []
    idx = ndimage.maximum_position(counts, labels, index=np.arange(1, n + 1))
    return [(int(a), int(b)) for a, b in idx]


SPIKE_STAR_MAG = 20.0  # point sources brighter than this have long diffraction spikes (ASSUMPTION)
SPIKE_ALIGN_DEG = 7.0


def spike_radius(mag, cap: float = 20.0) -> np.ndarray:
    """Spike length (arcsec) of a point source of magnitude ``mag``: 3" x 10^(0.2 (20 - m)),
    clipped to 3" - ``cap`` (the D-027 bright-star mask form; ASSUMPTION)."""
    return np.clip(3.0 * 10 ** (0.2 * (SPIKE_STAR_MAG - np.asarray(mag, float))), 3.0, cap)


# Saturated / off-mosaic stars (D-043): their spikes reach 37" in Abell 370. Catalogued stars
# keep the 20" cap: their catalogue magnitudes are unreliable when the core is clipped, and a
# longer cap on them re-creates the El Gordo over-veto that D-034 fixed (ASSUMPTION).
EXTERNAL_SPIKE_CAP = 60.0
SPIKE_STAR_DEDUP_ARCSEC = (
    1.0  # an external star this close to a catalogued bright point source is the same star
)
_STAR_COLUMNS = {
    "ra": ("ra", "RA_ICRS", "RA"),
    "dec": ("dec", "DE_ICRS", "DEC"),
    "mag": ("mag", "phot_g_mean_mag", "Gmag"),
}


def read_spike_stars(path: Path) -> Table:
    """External bright stars (e.g. Gaia DR3) as ``ra, dec, mag`` from an ECSV/FITS table.

    The magnitude may be ``mag``, ``phot_g_mean_mag`` or ``Gmag``. Gaia G (Vega) is used directly
    in the AB spike-length law, with no colour term, so red stars' spikes are underestimated by up
    to about 2x (ASSUMPTION). Masked values become NaN and those rows are ignored."""
    t = Table.read(path)
    cols = {}
    for key, names in _STAR_COLUMNS.items():
        found = next((c for c in names if c in t.colnames), None)
        if found is None:
            raise ValueError(f"{path}: no {key} column (expected one of {', '.join(names)})")
        col = t[found]
        cols[key] = np.asarray(col.filled(np.nan) if hasattr(col, "filled") else col, float)
    out = Table(cols)
    out.meta.update(provenance=schema.Provenance.OBSERVED.value, source=str(path))
    return out


def spike_segments(src: Table, shapes: Table, extra_stars: Table | None = None) -> np.ndarray:
    """Elongated sources that are diffraction-spike segments (El Gordo, D-034).

    A segment lies within ``spike_radius`` of a bright point source (not ``is_extended``,
    brighter than ``SPIKE_STAR_MAG``), its major axis is within ``SPIKE_ALIGN_DEG`` of the
    direction to that source, and that direction is one of the field's spike axes. JWST spikes
    sit at fixed position angles for one pointing: a hexagonal set (theta, theta + 60,
    theta + 120 deg) plus the weaker axis at theta + 90 deg. theta is estimated from all aligned
    pairs (mode of PA mod 60). Genuine radial arcs around a compact source survive unless they
    happen to lie on a spike axis.

    ``extra_stars`` (``ra, dec, mag``; e.g. Gaia DR3 from :func:`read_spike_stars`) adds stars the
    catalogue misses because they are saturated or off the mosaic, with spikes up to
    ``EXTERNAL_SPIKE_CAP`` (D-043)."""
    mag = np.asarray(shapes["mag"], float)
    bright = ~np.asarray(shapes["is_extended"], bool) & np.isfinite(mag) & (mag < SPIKE_STAR_MAG)
    ra_s = list(np.asarray(shapes["ra"], float)[bright])
    dec_s = list(np.asarray(shapes["dec"], float)[bright])
    radius = list(spike_radius(mag[bright]))
    if extra_stars is not None and len(extra_stars):
        em = np.asarray(extra_stars["mag"], float)
        ok = np.isfinite(em) & (em < SPIKE_STAR_MAG)
        if ra_s and ok.any():  # a star already in the catalogue is seeded once, from the catalogue
            ce = SkyCoord(extra_stars["ra"], extra_stars["dec"], unit="deg")
            _, sep, _ = ce.match_to_catalog_sky(SkyCoord(ra_s, dec_s, unit="deg"))
            ok &= sep.arcsec > SPIKE_STAR_DEDUP_ARCSEC
        ra_s += list(np.asarray(extra_stars["ra"], float)[ok])
        dec_s += list(np.asarray(extra_stars["dec"], float)[ok])
        radius += list(spike_radius(em[ok], cap=EXTERNAL_SPIKE_CAP))
    out = np.zeros(len(src), bool)
    if not ra_s or not len(src):
        return out
    cs = SkyCoord(src["ra"], src["dec"], unit="deg")
    ct = SkyCoord(ra_s, dec_s, unit="deg")
    radius = np.asarray(radius)
    pa_src = np.asarray(src["pa_obs"], float)
    cand, cand_pa = [], []
    for k in range(len(ct)):
        sep = cs.separation(ct[k]).arcsec
        near = (sep <= radius[k]) & (sep > 0.5)
        if not near.any():
            continue
        pa_to_star = np.mod(cs[near].position_angle(ct[k]).deg, 180.0)
        aligned = lensmodel.axis_offset_deg(pa_src[near], pa_to_star) <= SPIKE_ALIGN_DEG
        cand.extend(np.flatnonzero(near)[aligned])
        cand_pa.extend(pa_to_star[aligned])
    if not cand:
        return out
    cand, cand_pa = np.asarray(cand), np.asarray(cand_pa)
    hist, edges = np.histogram(np.mod(cand_pa, 60.0), bins=60, range=(0.0, 60.0))
    smooth = hist + np.roll(hist, 1) + np.roll(hist, -1)
    theta = float(edges[int(np.argmax(smooth))] + 0.5)
    axes = np.array([theta, theta + 60.0, theta + 120.0, theta + 90.0]) % 180.0
    on_axis = (
        np.min([lensmodel.axis_offset_deg(cand_pa, np.full(len(cand_pa), a)) for a in axes], axis=0)
        <= SPIKE_ALIGN_DEG
    )
    out[cand[on_axis]] = True
    return out


def radial_candidates(model, shapes: Table, args, extra: Table | None = None) -> tuple[Table, dict]:
    """The ``radial`` screen's arc selection: elongated, not spike segments, behind the lens (when a
    photo-z says so), ``anti`` to the model's stretch and not predicted radial by the model.

    ``shapes`` is a :func:`lens_consistency.load_shapes` table (with photo-z columns if attached);
    ``args`` carries the CLI thresholds. Returns the selected arcs and the selection counts.
    Shared by ``cmd_radial`` and ``scripts/inject_radial.py`` (D-049), so injections see the same
    screen."""
    sel = (
        (np.asarray(shapes["ellipticity"], float) >= args.min_ellipticity)
        & (np.asarray(shapes["semimajor_px"], float) >= args.min_semimajor_px)
        & (np.asarray(shapes["snr"], float) >= args.min_snr)
    )
    src = shapes[sel]
    x, y = model.to_frame(src["ra"], src["dec"])
    src = src[np.hypot(x, y) <= args.max_radius]
    spike = spike_segments(src, shapes, extra)
    n_spike = int(spike.sum())
    src = src[~spike]
    ot = lc.orientation_table(model, src)
    # cluster members and foreground objects are not lensed: with a photo-z, keep only sources
    # it puts behind the lens (sources without a photo-z stay, on the redshift grid)
    if "z_phot" in ot.colnames:
        has_pz = np.isfinite(np.asarray(ot["z_phot"], float))
        not_background = has_pz & (np.asarray(ot["z_basis"]) != "photo-z")
        n_not_background = int(not_background.sum())
        ot = ot[~not_background]
    else:
        n_not_background = 0
    anti = ot[ot["orientation_class"] == "anti"]
    # the model's own radial arcs are ordinary: drop sources it predicts to be radially stretched
    # the largest predicted radial magnification over the redshifts the class used: the z grid,
    # or the photo-z range (z160, z_phot, z840) for sources with a background photo-z
    use_pz = np.asarray(anti["z_basis"]) == "photo-z"
    zsets = [np.full(len(anti), z) for z in lc.Z_GRID]
    if use_pz.any():
        zsets = [
            np.where(use_pz, np.asarray(anti[c], float), zg)
            for c, zg in zip(("z160", "z_phot", "z840"), zsets, strict=True)
        ]
    mu_r = np.max([radial_magnification(model, anti["ra"], anti["dec"], z) for z in zsets], axis=0)
    anti["mu_radial"] = mu_r
    n_mu_nan = int(np.sum(~np.isfinite(mu_r)))
    cand = anti[np.isfinite(mu_r) & (mu_r < MAX_RADIAL_MU)]
    counts = {
        "n_spike_segments_dropped": n_spike,
        "n_elongated": len(src),
        "n_not_background_dropped": n_not_background,
        "n_anti": len(anti),
        "n_anti_not_model_radial": len(cand),
        "n_mu_radial_nan": n_mu_nan,
    }
    return cand, counts


def radial_grid(max_radius: float, step: float) -> tuple[np.ndarray, np.ndarray]:
    """The screen's search grid (model frame, arcsec): a square of half-width ``max_radius``."""
    g = np.arange(-max_radius, max_radius + 1e-9, step)
    return np.meshgrid(g, g)


def anti_window_draw(pa_pred, rng) -> np.ndarray:
    """One null draw: each arc keeps its selection (anti, i.e. >= 60 deg from the predicted
    tangential direction ``pa_pred``) but takes a uniform random angle inside that window."""
    pa_pred = np.asarray(pa_pred, float)
    return np.mod(pa_pred + 90.0 + rng.uniform(-30.0, 30.0, len(pa_pred)), 180.0)


def add_radial_options(r: argparse.ArgumentParser) -> None:
    """The ``radial`` thresholds (ASSUMPTIONs) and their defaults."""
    r.add_argument("--min-ellipticity", type=float, default=0.5)
    r.add_argument("--min-semimajor-px", type=float, default=2.0)
    r.add_argument("--min-snr", type=float, default=10.0)
    r.add_argument("--max-radius", type=float, default=60.0, help="arcsec from the model centre")
    r.add_argument("--grid-arcsec", type=float, default=0.5)
    r.add_argument("--line-tol-arcsec", type=float, default=1.0)
    r.add_argument("--max-len-arcsec", type=float, default=15.0)
    r.add_argument("--min-lines", type=int, default=3)
    r.add_argument("--dark-radius-arcsec", type=float, default=1.0)
    r.add_argument("--n-random", type=int, default=200)
    r.add_argument("--seed", type=int, default=1)


def radial_defaults() -> dict:
    """The ``radial`` CLI defaults, e.g. for injection-recovery (D-049)."""
    ap = argparse.ArgumentParser()
    add_radial_options(ap)
    return vars(ap.parse_args([]))


def cmd_radial(args) -> dict:
    model, _, _ = lc.load_model(args.model)  # a Lenstool model or published deflection maps
    lc.apply_frame_offset(args.model, model)  # into the JWST frame (D-034, D-040)
    shapes = lc.load_shapes(args.catalog)
    if args.photoz:
        lc.attach_photoz(shapes, args.photoz)
    extra = read_spike_stars(args.spike_stars) if args.spike_stars else None
    cand, sel_counts = radial_candidates(model, shapes, args, extra)
    cx, cy = model.to_frame(cand["ra"], cand["dec"])
    gx, gy = radial_grid(args.max_radius, args.grid_arcsec)
    counts = line_counts(cx, cy, cand["pa_obs"], gx, gy, args.line_tol_arcsec, args.max_len_arcsec)
    peaks = convergence_peaks(counts, args.min_lines)
    cs = SkyCoord(shapes["ra"], shapes["dec"], unit="deg")
    rows = []
    for iy, ix in sorted(peaks, key=lambda p: -counts[p]):
        ra, dec = model.to_sky(gx[iy, ix], gy[iy, ix])
        sep = SkyCoord(ra, dec, unit="deg").separation(cs).arcsec
        k = int(np.argmin(sep))
        rows.append(
            {
                "x": float(gx[iy, ix]),
                "y": float(gy[iy, ix]),
                "ra": float(ra),
                "dec": float(dec),
                "n_lines": int(counts[iy, ix]),
                "nearest_label": int(shapes["label"][k]),
                "nearest_sep_arcsec": float(sep[k]),
                "nearest_snr": float(shapes["snr"][k]),
                "dark_centre": bool(sep[k] > args.dark_radius_arcsec),
            }
        )
    # null: random angles inside each arc's anti window, so lines that point at the mass centre
    # because they were selected as anti are in the null too
    rng = np.random.default_rng(args.seed)
    pa_t = np.asarray(cand["pa_pred_z2"], float)
    n_rand, max_rand = [], []
    for _ in range(args.n_random):
        pa_r = anti_window_draw(pa_t, rng)
        rc = line_counts(cx, cy, pa_r, gx, gy, args.line_tol_arcsec, args.max_len_arcsec)
        n_rand.append(len(convergence_peaks(rc, args.min_lines)))
        max_rand.append(int(rc.max()) if rc.size else 0)
    max_obs = int(counts.max()) if counts.size else 0
    for row in rows:  # how often random angles give a centre at least this strong
        row["p_random"] = (
            float(np.mean(np.asarray(max_rand) >= row["n_lines"])) if max_rand else None
        )
    peaks_t = Table(rows=rows) if rows else Table(names=("x", "y", "n_lines"))
    peaks_t.meta.update(model._meta())
    peaks_t.meta.update(provenance=schema.Provenance.DERIVED.value, source=str(args.catalog))
    dest = args.out / args.model
    lc._write(cand, dest / "radial_sources.ecsv")
    lc._write(peaks_t, dest / "radial_peaks.ecsv")
    summary = {
        "model": args.model,
        "catalog": str(args.catalog),
        "spike_stars": (
            {
                "path": str(args.spike_stars),
                "n": len(extra),
                "sha256": hashlib.sha256(Path(args.spike_stars).read_bytes()).hexdigest(),
            }
            if extra is not None
            else None
        ),
        **sel_counts,
        "n_peaks": len(rows),
        "n_dark_centre_peaks": int(sum(r["dark_centre"] for r in rows)),
        "random_peaks_mean": float(np.mean(n_rand)) if n_rand else None,
        "random_peaks_p95": float(np.percentile(n_rand, 95)) if n_rand else None,
        "max_lines": max_obs,
        "random_max_lines": {
            str(k): int(v) for k, v in zip(*np.unique(max_rand, return_counts=True), strict=True)
        },
        "p_random_max": float(np.mean(np.asarray(max_rand) >= max_obs)) if max_rand else None,
        "peaks": rows,
        "assumptions": {
            k: getattr(args, k)
            for k in (
                "min_ellipticity",
                "min_semimajor_px",
                "min_snr",
                "max_radius",
                "grid_arcsec",
                "line_tol_arcsec",
                "max_len_arcsec",
                "min_lines",
                "dark_radius_arcsec",
                "n_random",
            )
        }
        | {"max_radial_mu": MAX_RADIAL_MU},
    }
    (dest / "radial.json").write_text(json.dumps(summary, indent=1))
    return summary


# ------------------------------------------------------------------------------------ shear (D-050)

SHEAR_MAX_G = 0.5  # rows where the cluster's reduced shear exceeds this are masked (ASSUMPTION)
SCHIRMER_XC = 0.15  # Schirmer et al. (2007) filter core / aperture radius (ASSUMPTION)


def lensable_mask(shapes: Table, z_lens: float, z_margin: float = 0.1) -> np.ndarray:
    """Rows a lens at ``z_lens`` can lens: no photo-z, or one that puts the source behind
    (``orientation_table``'s z160 > z_lens + margin), minus bright point sources (stars)."""
    mag = np.asarray(shapes["mag"], float)
    star = ~np.asarray(shapes["is_extended"], bool) & np.isfinite(mag) & (mag < SPIKE_STAR_MAG)
    ok = ~star
    if "z_phot" in shapes.colnames:
        has_pz = np.isfinite(np.asarray(shapes["z_phot"], float))
        z160 = np.asarray(shapes["z160"], float)
        ok &= ~has_pz | (np.isfinite(z160) & (z160 > z_lens + z_margin))
    return ok


def aperture_filter(x, kind: str, x_min: float) -> np.ndarray:
    """Weight Q(x) of a source at x = θ/R from the centre (zero outside ``x_min`` <= x <= 1).

    ``schirmer``: Q_TANH of Schirmer et al. (2007, eqs. 15-16): E(x) tanh(x/x_c)/(x/x_c), with
    E(x) = 1/(1 + exp(6 - 150x) + exp(-47 + 50x)); x_c = 0.15 (ASSUMPTION; the paper scans 7);
    ``pointmass``: 1/x², the matched filter for γ_t ∝ θ⁻² (point lens, D-047 W1);
    ``tophat``: 1."""
    x = np.asarray(x, float)
    if kind == "schirmer":
        xs = np.maximum(x, 1e-6) / SCHIRMER_XC
        with np.errstate(over="ignore"):
            q = np.tanh(xs) / xs / (1.0 + np.exp(6.0 - 150.0 * x) + np.exp(-47.0 + 50.0 * x))
    elif kind == "pointmass":
        q = 1.0 / np.maximum(x, x_min) ** 2
    elif kind == "tophat":
        q = np.ones_like(x)
    else:
        raise ValueError(f"unknown filter {kind!r}")
    return np.where((x >= x_min) & (x <= 1.0), q, 0.0)


def intrinsic_ellipticity(a, b, pa_deg, psf_sigma: float) -> np.ndarray:
    """Complex ellipticity ε = |ε| exp(2i PA) (PA east of north) of PSF-deconvolved second moments.

    ``a``, ``b``: semimajor and semiminor sigmas (px) of the observed moments; an isotropic
    Gaussian PSF of ``psf_sigma`` px is subtracted (ASSUMPTION, as in ``inject_radial``). Rows whose
    deconvolved moments are not positive definite (unresolved) are NaN. |ε| = (1 - q)/(1 + q)
    estimates the reduced shear of a randomly oriented source."""
    pa = np.deg2rad(np.asarray(pa_deg, float))
    a2, b2 = np.asarray(a, float) ** 2, np.asarray(b, float) ** 2
    s, c = np.sin(pa), np.cos(pa)
    qee = a2 * s * s + b2 * c * c - psf_sigma**2
    qnn = a2 * c * c + b2 * s * s - psf_sigma**2
    qen = (a2 - b2) * s * c
    tr, det = qee + qnn, qee * qnn - qen * qen
    with np.errstate(invalid="ignore", divide="ignore"):
        chi = ((qnn - qee) + 2j * qen) / tr
        eps = chi / (1.0 + np.sqrt(np.maximum(1.0 - np.abs(chi) ** 2, 0.0)))
    return np.where((tr > 0) & (det > 0), eps, np.nan + 0j)


def remove_cluster_shear(eps: np.ndarray, g: np.ndarray) -> np.ndarray:
    """Source-plane ellipticity (ε - g) / (1 - g* ε) for |g| < 1 (Seitz & Schneider 1997)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        return (eps - g) / (1.0 - np.conj(g) * eps)


def cluster_reduced_shear(model, ra, dec, z_s) -> tuple[np.ndarray, np.ndarray]:
    """Complex reduced shear g = |g| exp(2i PA) of the cluster model at each row's ``z_s``, and κ.

    PA is the model's stretching axis (``tangential_pa``); |g| = γ / |1 - κ|."""
    pred = model.evaluate(ra, dec, z_s)
    g = np.asarray(pred["reduced_shear"], float)
    pa = np.deg2rad(np.asarray(pred["tangential_pa"], float))
    return np.where(np.isfinite(pa), g * np.exp(2j * pa), 0j), np.asarray(pred["kappa"], float)


def shear_responsivity(eps: np.ndarray, g: np.ndarray) -> tuple[float, float]:
    """Slope R (and its error) of the measured ellipticity's component along the cluster model's
    shear, Re(ε exp(-i arg g)), against |g|, a least-squares line through the origin.

    Isophotal catalogue moments respond to shear with R < 1 (0.4–0.6 in MACS0416 and Abell 2744,
    D-053); subtracting the full g would leave -(1 - R) g, a radial pattern of the W1 sign."""
    ag = np.abs(g)
    ok = np.isfinite(eps) & (ag > 0)
    if ok.sum() < 20:  # R = 1 would leave the W1-signed residual; refuse instead
        raise ValueError(f"{int(ok.sum())} rows cannot calibrate the shear responsivity (need 20)")
    y = np.real(eps[ok] * np.exp(-1j * np.angle(g[ok])))
    x = ag[ok]
    r = float(np.dot(x, y) / np.dot(x, x))
    err = float(np.std(y - r * x) / np.sqrt(np.dot(x, x)))
    return r, err


def shear_sources(
    model,
    shapes: Table,
    psf_sigma: float,
    min_snr: float,
    max_g: float = SHEAR_MAX_G,
    responsivity: float | None = None,
    extra: Table | None = None,
    spike_veto: bool = True,
) -> tuple[np.ndarray, dict]:
    """Rows usable by the shear screen and their cluster-corrected ellipticities.

    Returns (``e``, counts): ``e`` is complex, NaN for unused rows. Used: lensable (background or
    no photo-z, not a star), S/N >= ``min_snr``, resolved after PSF deconvolution, and the cluster
    model's |g| < ``max_g`` and κ < 1 at the row's photo-z (z = 2 without one). The removed
    shear is R g, with R from :func:`shear_responsivity` on the used rows unless
    ``responsivity`` is given. Diffraction-spike segments (:func:`spike_segments`, with
    ``extra`` stars; D-043) are dropped: spikes point radially at their star, the W1 sign.
    ``counts["e_all"]`` holds the corrected ε of every resolved row with κ < 1 and |g| < 1
    (cuts on S/N, ``max_g`` and lensability not applied), for injected sources."""
    snr = np.asarray(shapes["snr"], float)
    a = np.asarray(shapes["semimajor_px"], float)
    b = a * (1.0 - np.asarray(shapes["ellipticity"], float))
    eps = intrinsic_ellipticity(a, b, shapes["pa_obs"], psf_sigma)
    z = np.full(len(shapes), 2.0)
    if "z_phot" in shapes.colnames:
        zp = np.asarray(shapes["z_phot"], float)
        z = np.where(np.isfinite(zp), zp, 2.0)
    g, kappa = cluster_reduced_shear(model, shapes["ra"], shapes["dec"], z)
    lens = lensable_mask(shapes, model.z_lens)
    spike = spike_segments(shapes, shapes, extra) if spike_veto else np.zeros(len(shapes), bool)
    lens &= ~spike
    bright = lens & (snr >= min_snr)
    resolved = bright & np.isfinite(eps)
    weak = resolved & (np.abs(g) < max_g) & (kappa < 1.0)
    r_err = float("nan")
    if responsivity is None:
        responsivity, r_err = shear_responsivity(np.where(weak, eps, np.nan), g)
    e_all = remove_cluster_shear(eps, responsivity * g)
    e_all = np.where(np.isfinite(eps) & (np.abs(g) < 1.0) & (kappa < 1.0), e_all, np.nan + 0j)
    e = np.where(weak, e_all, np.nan + 0j)
    counts = {
        "n_rows": len(shapes),
        "n_spike_segments": int(spike.sum()),
        "n_lensable": int(lens.sum()),
        "n_snr": int(bright.sum()),
        "n_resolved": int(resolved.sum()),
        "n_weak_shear": int(weak.sum()),
        "responsivity": float(responsivity),
        "responsivity_err": r_err,
    }
    counts["e_all"] = e_all
    return e, counts


class ApertureMass:
    """Catalogue aperture-mass statistic (Schneider 1996) on a grid of centres (model frame).

    For centre c, S_c = -Σ Q e_t / sqrt(Σ Q² |e|² / 2): the S/N of *negative* tangential shear
    (radial alignment), which a W1 lens (ε < 0) produces; S_× uses the cross component. Sums run
    over sources with ``r_min`` <= θ <= ``radius`` and centres with fewer than ``min_n`` sources are
    NaN. Positions are model-frame arcsec (x West, y North); ellipticities are |e| exp(2i PA)."""

    def __init__(self, x, y, cx, cy, radius: float, kind: str, r_min: float, min_n: int = 5):
        from scipy.sparse import csr_matrix
        from scipy.spatial import cKDTree

        x, y = np.asarray(x, float), np.asarray(y, float)
        cx, cy = np.ravel(np.asarray(cx, float)), np.ravel(np.asarray(cy, float))
        self.n_centres, self.min_n = len(cx), min_n
        rows = cols = np.zeros(0, int)
        if len(x):
            hits = cKDTree(np.c_[x, y]).query_ball_point(np.c_[cx, cy], radius)
            lens_ = np.fromiter((len(h) for h in hits), int, len(hits))
            rows = np.repeat(np.arange(len(hits)), lens_)
            cols = np.fromiter((j for h in hits for j in h), int, int(lens_.sum()))
        de, dn = -(x[cols] - cx[rows]), y[cols] - cy[rows]  # East, North of the centre
        theta = np.hypot(de, dn)
        q = aperture_filter(theta / radius, kind, r_min / radius)
        keep = q > 0
        rows, cols, q, de, dn = rows[keep], cols[keep], q[keep], de[keep], dn[keep]
        phi = np.arctan2(de, dn)  # position angle of the source seen from the centre
        shape = (self.n_centres, len(x))
        self.w = csr_matrix((q * np.exp(-2j * phi), (rows, cols)), shape=shape)
        self.w2 = csr_matrix((q * q, (rows, cols)), shape=shape)
        self.n = np.bincount(rows, minlength=self.n_centres)

    def snr(self, e: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(S, S_×) per centre for ellipticities ``e`` (one column per realisation allowed)."""
        e = np.asarray(e, complex)
        proj = self.w @ e  # Σ Q e exp(-2iφ) = -(Σ Q e_t) - i(Σ Q e_×)
        noise = np.sqrt(self.w2 @ (np.abs(e) ** 2) / 2.0)
        with np.errstate(invalid="ignore", divide="ignore"):
            s, sx = proj.real / noise, -proj.imag / noise
        bad = (self.n < self.min_n) if e.ndim == 1 else (self.n < self.min_n)[:, None]
        return np.where(bad, np.nan, s), np.where(bad, np.nan, sx)

    def null_max(self, e: np.ndarray, n_draws: int, rng, chunk: int = 50) -> np.ndarray:
        """Field maximum of S for ``n_draws`` random rotations of every source's ellipticity
        (positions and |e| kept)."""
        out = []
        for k in range(0, n_draws, chunk):
            m = min(chunk, n_draws - k)
            rot = np.exp(1j * rng.uniform(0, 2 * np.pi, (len(e), m)))
            s, _ = self.snr(e[:, None] * rot)
            out.append(np.nanmax(s, axis=0))
        return np.concatenate(out) if out else np.array([])


def psf_sigma_px(shapes: Table) -> float:
    """1st percentile of the minor-axis sigma, semimajor x (1 - ellipticity), over S/N > 50 rows:
    the narrowest objects are PSF-limited in their minor axis (``derived``)."""
    snr = np.asarray(shapes["snr"], float)
    b = np.asarray(shapes["semimajor_px"], float) * (1.0 - np.asarray(shapes["ellipticity"], float))
    return float(np.nanpercentile(b[snr > 50], 1))


def add_shear_options(s: argparse.ArgumentParser) -> None:
    """The ``shear`` thresholds (ASSUMPTIONs, D-050) and their defaults."""
    s.add_argument("--aperture-arcsec", type=float, default=10.0, help="aperture radius R")
    s.add_argument("--filter", choices=("pointmass", "schirmer", "tophat"), default="schirmer")
    s.add_argument("--r-min-arcsec", type=float, default=1.0, help="inner radius")
    s.add_argument("--min-sources", type=int, default=5, help="per aperture")
    s.add_argument("--min-snr", type=float, default=10.0)
    s.add_argument("--max-g", type=float, default=SHEAR_MAX_G)
    s.add_argument(
        "--responsivity", type=float, default=None, help="shear responsivity R (default: measured)"
    )
    s.add_argument("--max-radius", type=float, default=60.0, help="grid half-width (arcsec)")
    s.add_argument("--grid-arcsec", type=float, default=1.0)
    s.add_argument("--n-random", type=int, default=200)
    s.add_argument("--seed", type=int, default=1)


def shear_defaults() -> dict:
    """The ``shear`` CLI defaults, e.g. for injection-recovery."""
    ap = argparse.ArgumentParser()
    add_shear_options(ap)
    return vars(ap.parse_args([]))


def shear_screen(model, shapes: Table, args, psf_sigma: float, extra: Table | None = None) -> dict:
    """Run the shear screen; returns the map, null and peaks (all ``derived``)."""
    e, counts = shear_sources(
        model,
        shapes,
        psf_sigma,
        args.min_snr,
        args.max_g,
        getattr(args, "responsivity", None),
        extra,
    )
    e_all = counts.pop("e_all")
    use = np.isfinite(e)
    x, y = model.to_frame(np.asarray(shapes["ra"])[use], np.asarray(shapes["dec"])[use])
    gx, gy = radial_grid(args.max_radius, args.grid_arcsec)
    ap = ApertureMass(
        x, y, gx, gy, args.aperture_arcsec, args.filter, args.r_min_arcsec, args.min_sources
    )
    s, sx = ap.snr(e[use])
    null = ap.null_max(e[use], args.n_random, np.random.default_rng(args.seed))
    return {
        "counts": counts,
        "gx": gx,
        "gy": gy,
        "s": s.reshape(gx.shape),
        "s_cross": sx.reshape(gx.shape),
        "null_max": null,
        "e": e,
        "e_all": e_all,
        "use": use,
        "ap": ap,
        "xy": (x, y),
    }


def cross_p_values(s_cross: np.ndarray, null_max: np.ndarray) -> dict:
    """B-mode check: under random rotations S_× and -S_× have the same field-maximum distribution
    as S, so ``null_max`` gives the p-value of each cross extreme. An E-mode peak no rarer than
    the B-mode extremes is not distinguishable from shape systematics."""
    out = {}
    for name, v in (("plus", s_cross), ("minus", -s_cross)):
        m = float(np.nanmax(v)) if np.isfinite(v).any() else float("nan")
        out[f"s_cross_{name}_max"] = m
        ok = len(null_max) and np.isfinite(m)
        out[f"p_random_cross_{name}"] = float(np.mean(null_max >= m)) if ok else None
    return out


def cmd_shear(args) -> dict:
    model, _, _ = lc.load_model(args.model)
    lc.apply_frame_offset(args.model, model)
    shapes = lc.load_shapes(args.catalog)
    if args.photoz:
        lc.attach_photoz(shapes, args.photoz)
    psf = psf_sigma_px(shapes)
    extra = read_spike_stars(args.spike_stars) if args.spike_stars else None
    res = shear_screen(model, shapes, args, psf, extra)
    s, null, gx, gy = res["s"], res["null_max"], res["gx"], res["gy"]
    smax = float(np.nanmax(s)) if np.isfinite(s).any() else float("nan")
    p_max = float(np.mean(null >= smax)) if len(null) and np.isfinite(smax) else None
    rows = []
    thr = np.nanpercentile(null, 50) if len(null) else np.inf
    for iy, ix in convergence_peaks(np.nan_to_num(s, nan=-np.inf), thr):
        ra, dec = model.to_sky(gx[iy, ix], gy[iy, ix])
        rows.append(
            {
                "x": float(gx[iy, ix]),
                "y": float(gy[iy, ix]),
                "ra": float(ra),
                "dec": float(dec),
                "s": float(s[iy, ix]),
                "s_cross": float(res["s_cross"][iy, ix]),
                "p_random": float(np.mean(null >= s[iy, ix])),
            }
        )
    rows.sort(key=lambda r: -r["s"])
    cols = ("x", "y", "ra", "dec", "s", "s_cross", "p_random")
    peaks_t = Table(rows=rows, names=cols) if rows else Table(names=cols)
    peaks_t.meta.update(model._meta())
    peaks_t.meta.update(provenance=schema.Provenance.DERIVED.value, source=str(args.catalog))
    dest = args.out / args.model
    lc._write(peaks_t, dest / "shear_peaks.ecsv")
    finite = np.isfinite(s)
    summary = {
        "model": args.model,
        "catalog": str(args.catalog),
        "photoz": str(args.photoz) if args.photoz else None,
        "psf_sigma_px": psf,
        **res["counts"],
        "n_centres_valid": int(finite.sum()),
        "area_arcsec2": float(finite.sum() * args.grid_arcsec**2),
        "s_max": smax,
        "p_random_max": p_max,
        **cross_p_values(res["s_cross"], null),
        "null_max_p50_p95": (
            [float(v) for v in np.percentile(null, [50, 95])] if len(null) else None
        ),
        "n_peaks_above_null_median": len(rows),
        "peaks": rows[:20],
        "assumptions": {
            k: getattr(args, k)
            for k in (
                "aperture_arcsec",
                "filter",
                "r_min_arcsec",
                "min_sources",
                "min_snr",
                "max_g",
                "responsivity",
                "max_radius",
                "grid_arcsec",
                "n_random",
            )
        },
    }
    (dest / "shear.json").write_text(json.dumps(summary, indent=1))
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--model", choices=sorted(lc.MODELS), default="smacs0723-iclv2")
    ap.add_argument("--out", type=Path, default=paths.outputs_dir() / "exotic_screens")
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fluxratio", help="under- or overluminous catalogued images (two bands)")
    f.add_argument("--band1", required=True, help="_i2d URI or path (S3 by byte range)")
    f.add_argument("--band2", required=True, help="_i2d URI or path (S3 by byte range)")
    r = sub.add_parser("radial", help="anti-tangential arcs converging on a dark centre")
    r.add_argument("--catalog", type=Path, required=True, help="JWST pipeline _cat.ecsv")
    r.add_argument("--photoz", type=Path, help="eazy zout FITS table (e.g. DJA)")
    r.add_argument(
        "--spike-stars",
        type=Path,
        help="bright-star table (ra, dec, mag/Gmag; e.g. Gaia DR3) for the spike veto (D-043)",
    )
    add_radial_options(r)
    s = sub.add_parser("shear", help="negative tangential shear around a dark centre (W1, D-050)")
    s.add_argument("--catalog", type=Path, required=True, help="JWST pipeline _cat.ecsv")
    s.add_argument("--photoz", type=Path, help="eazy zout FITS table (e.g. DJA)")
    s.add_argument("--spike-stars", type=Path, help="bright-star table for the spike veto (D-043)")
    add_shear_options(s)
    args = ap.parse_args(argv)
    (args.out / args.model).mkdir(parents=True, exist_ok=True)
    summary = {"fluxratio": cmd_fluxratio, "radial": cmd_radial, "shear": cmd_shear}[args.cmd](args)
    print(json.dumps(summary, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
