"""Arc curvature against the published models: does each long background arc bend around the mass?

A tangential arc near a cluster critical curve is concave toward the mass, so its centre of
curvature lies on the side the model's deflection points away from (``-alpha``; for a point mass
alpha points outward). An arc that is significantly concave **away** from the mass is a lens-model
residual candidate; ordinary explanations are saddle points between clumps (the local
mass direction is not ``-alpha`` there), radial arcs, blends with neighbours and intrinsically
curved galaxies (spiral arms).

For every background, strong-shear source of a stored ``lens_consistency.py arcs`` table
(``results/arcs/<model>/arcs.ecsv``) with ``semimajor_px >= --min-semimajor-px`` and
``ellipticity >= --min-e``, a ``--size`` arcsec F200W cutout of the same observation's ``_i2d``
(S3 byte range, :func:`jwst_anomaly.cutouts.make_cutouts`) is thresholded at ``--nsigma`` above
the border background; the connected component at the centre is the arc. In the frame of its
flux-weighted major axis (u along, v across) a weighted parabola v = a + b u + c u^2 gives the
curvature side sign(c) v-hat. ``agree`` is True when that side points toward ``-alpha``
(z = 2; the direction does not depend on z for a single-plane model). A source counts only if
|c| / sigma_c >= ``--min-sig`` and the arc spans ``>= --min-length`` arcsec (ASSUMPTIONS).

Sensitivity: an arc bent around a centre at distance r (``r_arcsec``, the circular-lens
approximation; ASSUMPTION) has c_pred = 1 / (2 r), a sagitta of L^2 / (8 r) over length L. A source
is ``testable`` only when c_pred >= ``--min-sig`` sigma_c; elsewhere a "reverse" is shape noise
(intrinsic bends, knots, blends), not a lens-model test.

Known-case check: the aligned background arcs; an arc tool that works must give agree for most
of them. Outputs ``results/arcs_curvature/<model>/curvature.ecsv`` and ``summary.json`` (derived;
the pixels are observed).

    python scripts/arcs_curvature.py --model abell370-canucs
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS, FITSFixedWarning
from scipy import ndimage

from jwst_anomaly import paths, schema
from jwst_anomaly.cutouts import make_cutouts

_DIR = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("lens_consistency", _DIR / "lens_consistency.py")
lc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(lc)

ASSUMPTIONS = (
    "min_semimajor_px",
    "min_e",
    "size",
    "nsigma",
    "smooth_px",
    "radius_semimajor",
    "bkg_px",
    "min_sig",
    "min_length",
)
S3_URI = "s3://stpubdata/jwst/public/"

# F200W i2d of the observation whose catalogue each stored ``arcs`` run used.
IMAGES = {
    "abell2744-bergamini23": "jw02561-o001_t003_nircam_clear-f200w_i2d.fits",
    "elgordo-caminha23": "jw01176-o241_t012_nircam_clear-f200w_i2d.fits",
    "macs0416-canucs": "jw01208-o004_t002_nircam_clear-f200w_i2d.fits",
    "macs1149-cats": "jw01208-o008_t004_nircam_clear-f200w_i2d.fits",
    "abell370-canucs": "jw01208-o002_t001_nircam_clear-f200w_i2d.fits",
}


def image_uri(filename: str) -> str:
    prog, obs = filename[2:7], filename.split("-")[1].split("_")[0]
    return f"{S3_URI}jw{prog}/L3/t/{obs}/{filename}"


def arc_pixels(sci, nsigma: float, smooth_px: float, radius_px: float, bkg_px: int):
    """Boolean mask of the connected component above threshold at (or nearest) the centre.

    A ``bkg_px`` median filter removes cluster-galaxy and ICL light on larger scales; the noise is
    the 1.4826 MAD of the residual border pixels. Only pixels within ``radius_px`` of the centre
    count, so neighbours that touch the arc cannot dominate the fit. Returns ``(mask, img)``
    with ``img`` the background-subtracted cutout; the mask is empty when nothing near the centre
    passes.
    """
    img = np.where(np.isfinite(sci), sci, 0.0)
    img = img - ndimage.median_filter(img, size=bkg_px)
    sm = ndimage.gaussian_filter(img, smooth_px) if smooth_px > 0 else img
    border = np.concatenate([sm[0], sm[-1], sm[1:-1, 0], sm[1:-1, -1]])
    bkg = np.median(border)
    sigma = 1.4826 * np.median(np.abs(border - bkg))
    yy, xx = np.indices(img.shape)
    cy, cx = (np.array(img.shape) - 1) // 2
    inside = (yy - cy) ** 2 + (xx - cx) ** 2 <= radius_px**2
    lab, n = ndimage.label((sm - bkg > nsigma * sigma) & inside)
    if n == 0:
        return np.zeros_like(img, bool), img
    if lab[cy, cx] == 0:
        ys, xs = np.nonzero(lab)
        k = np.argmin((ys - cy) ** 2 + (xs - cx) ** 2)
        if (ys[k] - cy) ** 2 + (xs[k] - cx) ** 2 > 9:  # nothing within 3 px of the centre
            return np.zeros_like(img, bool), img
        return lab == lab[ys[k], xs[k]], img
    return lab == lab[cy, cx], img


def fit_curvature(east, north, w) -> dict:
    """Weighted parabola in the major-axis frame of points ``(east, north)`` (arcsec).

    ``pa`` is the major axis (deg E of N, mod 180); ``c`` the curvature coefficient (1/arcsec) with
    its formal error from the weighted residual scatter; ``side_e``, ``side_n`` the unit vector to
    the centre of curvature; ``length`` the weighted 5-95 % extent along the major axis.
    """
    w = np.clip(np.asarray(w, float), 0, None)
    if w.sum() <= 0 or len(w) < 6:
        return {}
    e0, n0 = np.average(east, weights=w), np.average(north, weights=w)
    de, dn = east - e0, north - n0
    cov = np.cov(np.vstack([de, dn]), aweights=w)
    vals, vecs = np.linalg.eigh(cov)
    ue, un = vecs[:, 1]  # major axis (east, north components)
    ve, vn = -un, ue  # 90 deg counter-clockwise on the (east, north) plane
    u_c, v_c = de * ue + dn * un, de * ve + dn * vn
    order = np.argsort(u_c)
    cw = np.cumsum(w[order]) / w.sum()
    length = float(u_c[order][np.searchsorted(cw, 0.95)] - u_c[order][np.searchsorted(cw, 0.05)])
    sw = np.sqrt(w)
    a_mat = np.vstack([np.ones_like(u_c), u_c, u_c**2]).T * sw[:, None]
    coef, *_ = np.linalg.lstsq(a_mat, v_c * sw, rcond=None)
    resid = v_c - np.vstack([np.ones_like(u_c), u_c, u_c**2]).T @ coef
    n_eff = w.sum() ** 2 / np.sum(w**2)
    s2 = np.sum(w * resid**2) / w.sum() * n_eff / max(n_eff - 3, 1)
    cov_c = np.linalg.pinv(a_mat.T @ a_mat) * s2 * w.sum() / n_eff
    c = float(coef[2])
    sgn = 1.0 if c >= 0 else -1.0
    return {
        "pa": float(np.mod(np.degrees(np.arctan2(ue, un)), 180.0)),
        "c": c,
        "c_err": float(np.sqrt(max(cov_c[2, 2], 0.0))),
        "side_e": sgn * ve,
        "side_n": sgn * vn,
        "length": length,
    }


def measure(path: str, ra0: float, dec0: float, semimajor_px: float, args) -> dict:
    with fits.open(path) as hdul, warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        sci = hdul["SCI"].data.astype(float)
        wcs = WCS(hdul["SCI"].header)
    radius = args.radius_semimajor * semimajor_px
    mask, img = arc_pixels(sci, args.nsigma, args.smooth_px, radius, args.bkg_px)
    if mask.sum() < 6:
        return {"n_pix": int(mask.sum())}
    ys, xs = np.nonzero(mask)
    ra, dec = wcs.pixel_to_world_values(xs, ys)
    east = (ra - ra0) * np.cos(np.radians(dec0)) * 3600.0
    north = (dec - dec0) * 3600.0
    return {"n_pix": int(mask.sum())} | fit_curvature(east, north, img[ys, xs])


def classify(fit: dict, mass_e: float, mass_n: float, args) -> str:
    if "c" not in fit:
        return "no_arc"
    if fit["length"] < args.min_length or abs(fit["c"]) < args.min_sig * fit["c_err"]:
        return "flat"
    return "agree" if fit["side_e"] * mass_e + fit["side_n"] * mass_n > 0 else "reverse"


def run(args) -> dict:
    stored = Table.read(args.arcs_dir / args.model / "arcs.ecsv")
    sel = (
        (stored["population"] == "background")
        & stored["strong_shear"]
        & (np.asarray(stored["semimajor_px"], float) >= args.min_semimajor_px)
        & (np.asarray(stored["ellipticity"], float) >= args.min_e)
    )
    t = stored[sel]
    model, _, _ = lc.load_model(args.model)
    lc.apply_frame_offset(args.model, model)
    defl = model.deflection(t["ra"], t["dec"], np.full(len(t), 2.0))
    # model frame: x West, y North; mass direction = -alpha, east component = +alpha_x
    ax, ay = np.asarray(defl["alpha_x"], float), np.asarray(defl["alpha_y"], float)
    norm = np.hypot(ax, ay)
    mass_e, mass_n = ax / norm, -ay / norm
    targets = Table({"source_uid": [f"{args.model}-{int(v)}" for v in t["label"]]})
    targets["ra"], targets["dec"] = t["ra"], t["dec"]
    cut = make_cutouts(
        image_uri(IMAGES[args.model]),
        targets,
        size_arcsec=args.size,
        out_dir=args.cutout_dir,
        extensions=[],
        band="F200W",
    )
    rows = []
    for r, c, me, mn in zip(t, cut, mass_e, mass_n, strict=True):
        fit = (
            measure(c["path"], r["ra"], r["dec"], float(r["semimajor_px"]), args)
            if c["path"]
            else {"n_pix": 0}
        )
        rows.append(
            {
                "label": int(r["label"]),
                "ra": float(r["ra"]),
                "dec": float(r["dec"]),
                "orientation_class": str(r["orientation_class"]),
                "g_min": float(r["g_min"]),
                "semimajor_px": float(r["semimajor_px"]),
                "ellipticity": float(r["ellipticity"]),
                "quality_flag": str(c["quality_flag"]),
                "n_pix": fit["n_pix"],
                "pa_arc": fit.get("pa", np.nan),
                "length": fit.get("length", np.nan),
                "c": fit.get("c", np.nan),
                "c_sig": abs(fit["c"]) / fit["c_err"] if fit.get("c_err") else np.nan,
                "c_pred": 0.5 / float(r["r_arcsec"]),
                "testable": bool(
                    fit.get("c_err") and 0.5 / float(r["r_arcsec"]) >= args.min_sig * fit["c_err"]
                ),
                "cos_mass": fit["side_e"] * me + fit["side_n"] * mn if "c" in fit else np.nan,
                "curvature_class": classify(fit, me, mn, args),
                "path": str(c["path"]),
            }
        )
    out_t = Table(rows=rows) if rows else Table()
    out_t.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source=f"{IMAGES[args.model]} pixels + {stored.meta.get('source', 'arcs.ecsv')}",
        assumptions={k: getattr(args, k) for k in ASSUMPTIONS},
    )
    classes = {
        k: int(np.sum(out_t["curvature_class"] == k)) if rows else 0
        for k in ("agree", "reverse", "flat", "no_arc")
    }
    rev = out_t[out_t["curvature_class"] == "reverse"] if rows else []
    al = out_t[out_t["orientation_class"] == "aligned"] if rows else []
    summary = {
        "model": args.model,
        "image": IMAGES[args.model],
        "assumptions": {k: getattr(args, k) for k in ASSUMPTIONS},
        "n_selected": len(t),
        "n_testable": int(np.sum(out_t["testable"])) if rows else 0,
        "median_abs_c_over_pred": round(
            float(np.nanmedian(np.abs(out_t["c"]) / out_t["c_pred"])), 1
        )
        if rows
        else None,
        "classes": classes,
        "control_aligned": {
            k: int(np.sum(al["curvature_class"] == k)) if len(al) else 0
            for k in ("agree", "reverse", "flat", "no_arc")
        },
        "reverse": [
            {
                "label": int(r["label"]),
                "orientation_class": r["orientation_class"],
                "g_min": round(float(r["g_min"]), 2),
                "c_sig": round(float(r["c_sig"]), 1),
                "length": round(float(r["length"]), 2),
                "cos_mass": round(float(r["cos_mass"]), 2),
                "c_over_pred": round(float(abs(r["c"]) / r["c_pred"]), 1),
                "testable": bool(r["testable"]),
            }
            for r in rev
        ],
    }
    out = args.out / args.model
    out.mkdir(parents=True, exist_ok=True)
    out_t[[c for c in out_t.colnames if c != "path"]].write(out / "curvature.ecsv", overwrite=True)
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--model", required=True, choices=sorted(IMAGES))
    p.add_argument("--arcs-dir", type=Path, default=paths.repo_root() / "results" / "arcs")
    p.add_argument("--out", type=Path, default=paths.repo_root() / "results" / "arcs_curvature")
    p.add_argument("--cutout-dir", type=Path, default=paths.outputs_dir() / "arcs_curvature")
    p.add_argument("--min-semimajor-px", type=float, default=10.0)
    p.add_argument("--min-e", type=float, default=0.5)
    p.add_argument("--size", type=float, default=6.0, help="cutout side (arcsec)")
    p.add_argument("--nsigma", type=float, default=3.0)
    p.add_argument("--smooth-px", type=float, default=1.0)
    p.add_argument("--radius-semimajor", type=float, default=3.0, help="fit window / semimajor")
    p.add_argument("--bkg-px", type=int, default=61, help="median-filter box for halo light (px)")
    p.add_argument("--min-sig", type=float, default=3.0, help="|c| / sigma_c for a curved arc")
    p.add_argument("--min-length", type=float, default=0.5, help="arcsec, 5-95 %% extent")
    args = p.parse_args(argv)
    print(json.dumps(run(args), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
