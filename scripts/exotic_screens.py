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
    rate comes from the same search with the position angles randomised (``--n-random`` draws).
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
            others = lum[idx[idx != k]]
            ref = float(np.median(others))
            ratio[k] = lum[k] / ref
            rsnr[k] = abs(flux[k] - ref * abs(mu[k])) / err[k]
    out = Table({"system": systems, "lum_ratio": ratio, "ratio_snr": rsnr})
    out["usable"] = ok
    return out


def concentration(img, err, xx, yy) -> float:
    """Aperture flux at r = 0.2" over r = 0.4", both at the peak within 0.4" (local background)."""
    f_small, _ = lc.peak_flux(img, err, xx, yy)
    inner = np.where((np.hypot(xx, yy) <= 0.4) & np.isfinite(img), img, -np.inf)
    if not np.isfinite(inner).any() or not np.isfinite(f_small):
        return float("nan")
    k = int(np.argmax(inner))
    cx, cy = float(xx.flat[k]), float(yy.flat[k])
    r = np.hypot(xx - cx, yy - cy)
    ann = (r > 0.8) & (r < 1.2) & np.isfinite(img)
    if ann.sum() < 10:
        return float("nan")
    bkg = float(np.median(img[ann]))
    f_big = float(np.nansum(img[r <= 0.4] - bkg))
    return f_small / f_big if f_big > 0 else float("nan")


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
    files = lc.model_files(args.model)
    par = lensmodel.parse_lenstool_par(files["best.par"])
    model = lensmodel.LensModel.from_par(par)
    images = lensmodel.load_lenstool_images(files["arcs.dat"])
    bt = lensmodel.backtrace_images(model, images, par["z_m_limit"])
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
        fe = [lc.peak_flux(*st) for st in stamps]
        out[f"flux_{tag}"] = [f for f, _ in fe]
        out[f"err_{tag}"] = [e for _, e in fe]
        out[f"conc_{tag}"] = [concentration(*st) for st in stamps]
    # only unresolved images take part: a resolved arc's aperture flux does not scale with |μ|
    compact = np.asarray(out["conc_b1"], float) >= MIN_CONCENTRATION
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


def cmd_radial(args) -> dict:
    files = lc.model_files(args.model)
    par = lensmodel.parse_lenstool_par(files["best.par"])
    model = lensmodel.LensModel.from_par(par)
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
    ot = lc.orientation_table(model, src)
    anti = ot[ot["orientation_class"] == "anti"]
    # the model's own radial arcs are ordinary: drop sources it predicts to be radially stretched
    z = np.where(
        np.asarray(anti["z_basis"]) == "photo-z",
        np.asarray(anti["z_phot"], float) if "z_phot" in anti.colnames else 2.0,
        2.0,
    )
    mu_r = radial_magnification(model, anti["ra"], anti["dec"], z)
    anti["mu_radial"] = mu_r
    cand = anti[mu_r < MAX_RADIAL_MU]
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
    # false-alarm rate: the same search with the candidates' position angles randomised
    rng = np.random.default_rng(args.seed)
    n_rand, max_rand = [], []
    for _ in range(args.n_random):
        rc = line_counts(
            cx, cy, rng.uniform(0, 180, len(cx)), gx, gy, args.line_tol_arcsec, args.max_len_arcsec
        )
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
        "n_elongated": len(src),
        "n_anti": len(anti),
        "n_anti_not_model_radial": len(cand),
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
    r.add_argument("--n-random", type=int, default=50)
    r.add_argument("--seed", type=int, default=1)
    args = ap.parse_args(argv)
    (args.out / args.model).mkdir(parents=True, exist_ok=True)
    summary = {"fluxratio": cmd_fluxratio, "radial": cmd_radial}[args.cmd](args)
    print(json.dumps(summary, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
