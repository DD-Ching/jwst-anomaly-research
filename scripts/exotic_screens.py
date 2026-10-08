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

Every threshold is an ASSUMPTION (defaults below); outputs are ``derived`` with ``model_prediction``
inputs. Results go to ``outputs/exotic_screens/<model>/``.
"""

from __future__ import annotations

import argparse
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


def spike_radius(mag) -> np.ndarray:
    """Spike length (arcsec) of a point source of magnitude ``mag``: 3" x 10^(0.2 (20 - m)),
    clipped to 3-20" (the D-027 bright-star mask form; ASSUMPTION)."""
    return np.clip(3.0 * 10 ** (0.2 * (SPIKE_STAR_MAG - np.asarray(mag, float))), 3.0, 20.0)


def spike_segments(src: Table, shapes: Table) -> np.ndarray:
    """Elongated sources that are diffraction-spike segments (El Gordo, D-034).

    A segment lies within ``spike_radius`` of a bright point source (not ``is_extended``,
    brighter than ``SPIKE_STAR_MAG``), its major axis is within ``SPIKE_ALIGN_DEG`` of the
    direction to that source, and that direction is one of the field's spike axes. JWST spikes
    sit at fixed position angles for one pointing: a hexagonal set (theta, theta + 60,
    theta + 120 deg) plus the weaker axis at theta + 90 deg. theta is estimated from all aligned
    pairs (mode of PA mod 60). Genuine radial arcs around a compact source survive unless they
    happen to lie on a spike axis."""
    mag = np.asarray(shapes["mag"], float)
    bright = ~np.asarray(shapes["is_extended"], bool) & np.isfinite(mag) & (mag < SPIKE_STAR_MAG)
    stars = shapes[bright]
    out = np.zeros(len(src), bool)
    if not len(stars) or not len(src):
        return out
    cs = SkyCoord(src["ra"], src["dec"], unit="deg")
    ct = SkyCoord(stars["ra"], stars["dec"], unit="deg")
    radius = spike_radius(stars["mag"])
    pa_src = np.asarray(src["pa_obs"], float)
    cand, cand_pa = [], []
    for k in range(len(stars)):
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


def cmd_radial(args) -> dict:
    model, _, par = lc.load_model(args.model)  # a Lenstool model or published deflection maps
    if par is not None:
        lc.apply_frame_offset(args.model, model)
    shapes = lc.load_shapes(args.catalog)
    if args.photoz:
        lc.attach_photoz(shapes, args.photoz)
    sel = (
        (np.asarray(shapes["ellipticity"], float) >= args.min_ellipticity)
        & (np.asarray(shapes["semimajor_px"], float) >= args.min_semimajor_px)
        & (np.asarray(shapes["snr"], float) >= args.min_snr)
    )
    src = shapes[sel]
    x, y = model.to_frame(src["ra"], src["dec"])
    src = src[np.hypot(x, y) <= args.max_radius]
    spike = spike_segments(src, shapes)
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
    cx, cy = model.to_frame(cand["ra"], cand["dec"])
    g = np.arange(-args.max_radius, args.max_radius + 1e-9, args.grid_arcsec)
    gx, gy = np.meshgrid(g, g)
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
    # null: each arc keeps its selection (anti, i.e. >= 60 deg from the predicted tangential
    # direction) but takes a random angle inside that window, so lines that point at the mass
    # centre because they were selected as anti are in the null too
    rng = np.random.default_rng(args.seed)
    pa_t = np.asarray(cand["pa_pred_z2"], float)
    n_rand, max_rand = [], []
    for _ in range(args.n_random):
        pa_r = np.mod(pa_t + 90.0 + rng.uniform(-30.0, 30.0, len(cx)), 180.0)
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
        "n_spike_segments_dropped": n_spike,
        "n_elongated": len(src),
        "n_not_background_dropped": n_not_background,
        "n_anti": len(anti),
        "n_anti_not_model_radial": len(cand),
        "n_mu_radial_nan": n_mu_nan,
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
    args = ap.parse_args(argv)
    (args.out / args.model).mkdir(parents=True, exist_ok=True)
    summary = {"fluxratio": cmd_fluxratio, "radial": cmd_radial}[args.cmd](args)
    print(json.dumps(summary, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
