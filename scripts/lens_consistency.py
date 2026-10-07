"""Lens-model consistency checks against a published Lenstool model (M3, D-024).

Two subcommands, both writing to ``outputs/lens_consistency/<model>/`` (ECSV plus a JSON summary):

``validate``
    Checks the ported dPIE (``jwst_anomaly.lensmodel``) against the model's own published
    products before it is trusted on anything else:

    * the published convergence map (D_LS/D_S = 1): median, 95th and 99th percentile and maximum
      of |Δκ| on a regular sub-grid of its pixels;
    * the published multiple images (``arcs.dat``): each image traced to the source plane, the
      scatter around its system's mean source mapped back to the image plane, and the resulting
      χ² at the model's ``sigposArcsec`` (compared with the ``Chi2pos`` written in ``best.par``).

``arcs``
    Compares the major-axis orientation of elongated sources in a JWST pipeline catalog
    (``sky_orientation``, E of N) with the model's predicted stretching direction at
    z_s = 1, 2 and 4 and, when an eazy ``zout`` table is given, at each source's photo-z.

    * ``aligned``: within ``--aligned-deg`` of the prediction at every z_s;
    * ``anti``: at least ``--anti-deg`` away from it at every z_s;
    * ``mixed``: anything else.

    First, the known multiple images matched in the catalog check the orientation convention.
    Only background sources (photo-z 16th percentile above z_lens + ``--z-margin``) test the
    model, because cluster members and foreground galaxies are not lensed. ``anti`` background
    sources are candidates for inspection, not anomalies: intrinsic shapes, deblending
    fragments, photo-z failures and model error are the ordinary explanations to rule out.

Every threshold here is an ASSUMPTION. Inputs: the pinned model files (downloaded and verified by
sha256, ``SMACS0723_MAHLER22_ICLV2``), a level-3 ``_cat.ecsv`` and optionally a ``zout`` file.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import tarfile
import warnings
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS, FITSFixedWarning

from jwst_anomaly import lensmodel, paths, schema
from jwst_anomaly.photometry import fetch_catalog

MODELS = {
    "smacs0723-iclv2": {
        "files": lensmodel.SMACS0723_MAHLER22_ICLV2,
        "kappa_member": "tmp_k/0000_k.fits",
    },
}
Z_GRID = (1.0, 2.0, 4.0)


def model_files(name: str) -> dict[str, Path]:
    """Download (once) and verify the pinned files of model ``name``; extract its kappa map."""
    spec = MODELS[name]
    files = {key: fetch_catalog(url, sha) for key, (url, sha) in spec["files"].items()}
    archive = files.pop("kappa_0000")
    out_dir = archive.with_name(archive.name.split(".", 1)[0] + "_extracted")
    member = out_dir / spec["kappa_member"]
    if not member.exists():
        with tarfile.open(archive) as tar:
            tar.extract(spec["kappa_member"], out_dir, filter="data")
    files["kappa_map"] = member
    return files


def kappa_map_check(model: lensmodel.LensModel, map_path: Path, step: int = 5) -> dict:
    """|Δκ| between the analytic model and a published convergence map (D_LS/D_S = 1)."""
    with fits.open(map_path) as hdul:
        data = np.asarray(hdul[0].data, float)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FITSFixedWarning)
            wcs = WCS(hdul[0].header)
    jj, ii = np.mgrid[0 : data.shape[0] : step, 0 : data.shape[1] : step]
    ra, dec = wcs.pixel_to_world_values(ii, jj)
    x, y = model.to_frame(ra, dec)
    published = data[jj, ii]
    diff = np.abs(model.kappa_xy(x, y) - published)
    k = int(np.argmax(diff))
    return {
        "map": str(map_path),
        "grid_step_px": step,
        "n_points": int(diff.size),
        "extent_arcsec": [float(x.min()), float(x.max()), float(y.min()), float(y.max())],
        "median_abs_dkappa": float(np.median(diff)),
        "p95_abs_dkappa": float(np.percentile(diff, 95)),
        "p99_abs_dkappa": float(np.percentile(diff, 99)),
        "max_abs_dkappa": float(diff.flat[k]),
        "max_at_xy_arcsec": [float(x.flat[k]), float(y.flat[k])],
        "kappa_at_max": float(published.flat[k]),
        "median_ratio": float(np.median(model.kappa_xy(x, y) / published)),
    }


def chi2pos_from_par(path: Path) -> float | None:
    """``#Chi2pos:`` written by Lenstool at the top of ``best.par`` (None if absent)."""
    match = re.search(r"^#Chi2pos:\s*(\S+)", path.read_text(), re.M)
    return float(match.group(1)) if match else None


def backtrace_check(model, images, z_m_limit, sigpos: float, chi2_ref: float | None):
    """Back-trace table plus a summary (χ², rms, per-system rms, > 3σ images)."""
    bt = lensmodel.backtrace_images(model, images, z_m_limit)
    ok = np.isfinite(bt["dtheta_arcsec"])
    dtheta = np.asarray(bt["dtheta_arcsec"][ok])
    systems = {}
    for sys_id in np.unique(bt["system"][ok]):
        m = ok & (bt["system"] == sys_id)
        systems[str(sys_id)] = float(np.sqrt(np.mean(np.asarray(bt["dtheta_arcsec"][m]) ** 2)))
    outliers = bt[ok & (bt["dtheta_arcsec"] > 3 * sigpos)]
    summary = {
        "n_images": len(bt),
        "n_traced": int(ok.sum()),
        "rms_dtheta_arcsec": float(np.sqrt(np.mean(dtheta**2))),
        "median_dtheta_arcsec": float(np.median(dtheta)),
        "sigpos_arcsec": sigpos,
        "chi2_pos": float(np.sum((dtheta / sigpos) ** 2)),
        "chi2_pos_lenstool": chi2_ref,
        "rms_by_system": systems,
        "images_over_3sigma": [str(i) for i in outliers["image_id"]],
    }
    return bt, summary


def load_shapes(cat_path: Path) -> Table:
    """Shape columns of a JWST pipeline ``_cat.ecsv`` (observed)."""
    cat = Table.read(cat_path)
    sky = cat["sky_centroid"]
    with np.errstate(divide="ignore", invalid="ignore"):
        snr = np.asarray(cat["isophotal_flux"], float) / np.asarray(
            cat["isophotal_flux_err"], float
        )
    out = Table(
        {
            "label": np.asarray(cat["label"]),
            "ra": sky.ra.deg,
            "dec": sky.dec.deg,
            "pa_obs": np.mod(np.asarray(cat["sky_orientation"], float), 180.0),
            "ellipticity": np.asarray(cat["ellipticity"], float),
            "semimajor_px": np.asarray(cat["semimajor_sigma"], float),
            "area_px": np.asarray(cat["isophotal_area"], float),
            "snr": snr,
            "is_extended": np.asarray(cat["is_extended"], bool),
        }
    )
    out.meta.update(provenance=schema.Provenance.OBSERVED.value, source=str(cat_path))
    return out


def attach_photoz(shapes: Table, zout_path: Path, radius_arcsec: float = 0.3) -> None:
    """Add ``z_phot``, ``z160``, ``z840`` from the nearest eazy ``zout`` entry (NaN if none)."""
    zout = Table.read(zout_path)
    c1 = SkyCoord(shapes["ra"], shapes["dec"], unit="deg")
    c2 = SkyCoord(np.asarray(zout["ra"], float), np.asarray(zout["dec"], float), unit="deg")
    idx, sep, _ = c1.match_to_catalog_sky(c2)
    ok = sep.arcsec <= radius_arcsec
    for col in ("z_phot", "z160", "z840"):
        vals = np.asarray(zout[col], float)[idx]
        shapes[col] = np.where(ok & (vals > 0), vals, np.nan)
    shapes.meta["photoz_source"] = str(zout_path)


def orientation_table(
    model: lensmodel.LensModel,
    sources: Table,
    z_grid=Z_GRID,
    aligned_deg: float = 30.0,
    anti_deg: float = 60.0,
) -> Table:
    """Predicted stretching direction at each z_s and the offset from the observed axis."""
    out = sources.copy(copy_data=True)
    offsets, pas = [], []
    for z in z_grid:
        pred = model.evaluate(out["ra"], out["dec"], z)
        tag = f"z{z:g}"
        out[f"pa_pred_{tag}"] = pred["tangential_pa"]
        out[f"g_{tag}"] = pred["reduced_shear"]
        out[f"mu_{tag}"] = pred["magnification"]
        off = lensmodel.axis_offset_deg(out["pa_obs"], pred["tangential_pa"])
        out[f"offset_{tag}"] = off
        offsets.append(off)
        pas.append(np.asarray(pred["tangential_pa"], float))
    if "z_phot" in out.colnames:
        zp = np.asarray(out["z_phot"], float)
        good = np.isfinite(zp) & (zp > model.z_lens)
        off = np.full(len(out), np.nan)
        if good.any():
            pred = model.evaluate(out["ra"][good], out["dec"][good], zp[good])
            off[good] = lensmodel.axis_offset_deg(out["pa_obs"][good], pred["tangential_pa"])
        out["offset_zphot"] = off
    stack = np.vstack(offsets)
    out["offset_min"] = stack.min(axis=0)
    out["offset_max"] = stack.max(axis=0)
    spread = np.zeros(len(out))
    for a in range(len(pas)):
        for b in range(a + 1, len(pas)):
            spread = np.maximum(spread, lensmodel.axis_offset_deg(pas[a], pas[b]))
    out["pa_pred_spread"] = spread  # > 45 deg: the prediction flips near a critical curve
    cls = np.where(out["offset_max"] <= aligned_deg, "aligned", "mixed")
    cls = np.where(out["offset_min"] >= anti_deg, "anti", cls)
    out["orientation_class"] = cls
    out.meta.update(model._meta())
    out.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source=f"{sources.meta.get('source', 'catalog')} against {model.source}",
        z_grid=list(z_grid),
        assumptions={"aligned_deg": aligned_deg, "anti_deg": anti_deg},
    )
    return out


def _binom_p(k: int, n: int, p: float) -> float:
    """One-sided binomial P(X >= k)."""
    from scipy.stats import binomtest

    return float(binomtest(k, n, p, alternative="greater").pvalue) if n else float("nan")


def _write(table: Table, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    table.write(path, overwrite=True)


def cmd_validate(args) -> dict:
    files = model_files(args.model)
    par = lensmodel.parse_lenstool_par(files["best.par"])
    model = lensmodel.LensModel.from_par(par)
    sigpos = lensmodel.parse_lenstool_par(files["input.par"])["sigpos_arcsec"]
    images = lensmodel.load_lenstool_images(files["arcs.dat"])
    bt, bsum = backtrace_check(
        model, images, par["z_m_limit"], sigpos, chi2pos_from_par(files["best.par"])
    )
    summary = {
        "model": args.model,
        "model_sha256": model.sha256,
        "n_potentials": len(model.components),
        "kappa_map": kappa_map_check(model, files["kappa_map"], step=args.step),
        "backtrace": bsum,
    }
    out = args.out / args.model
    _write(bt, out / "backtrace.ecsv")
    (out / "validate.json").write_text(json.dumps(summary, indent=1))
    return summary


def cmd_arcs(args) -> dict:
    files = model_files(args.model)
    par = lensmodel.parse_lenstool_par(files["best.par"])
    model = lensmodel.LensModel.from_par(par)
    shapes = load_shapes(args.catalog)
    if args.photoz:
        attach_photoz(shapes, args.photoz)
    x, y = model.to_frame(shapes["ra"], shapes["dec"])
    shapes["r_arcsec"] = np.hypot(x, y)
    sel = (
        shapes["is_extended"]
        & (shapes["ellipticity"] >= args.min_ellipticity)
        & (shapes["semimajor_px"] >= args.min_semimajor_px)
        & (shapes["snr"] >= args.min_snr)
        & (shapes["r_arcsec"] <= args.max_radius)
    )
    table = orientation_table(
        model, shapes[sel], aligned_deg=args.aligned_deg, anti_deg=args.anti_deg
    )
    strong = np.asarray(table[f"g_z{Z_GRID[0]:g}"], float) >= args.min_shear
    table["strong_shear"] = strong
    if "z160" in table.colnames:
        zl = model.z_lens + args.z_margin
        pop = np.where(table["z160"] > zl, "background", "unknown")
        pop = np.where(table["z840"] < zl, "member_or_foreground", pop)
    else:
        pop = np.full(len(table), "unknown")
    table["population"] = pop

    # Convention check: catalogued multiple images, at their model redshift.
    images = lensmodel.load_lenstool_images(files["arcs.dat"])
    z_img = np.array(
        [
            z if z > 0 else par["z_m_limit"].get(str(s), np.nan)
            for z, s in zip(images["z"], images["system"], strict=True)
        ]
    )
    ci = SkyCoord(images["ra"], images["dec"], unit="deg")
    cs = SkyCoord(shapes["ra"], shapes["dec"], unit="deg")
    idx, sep, _ = ci.match_to_catalog_sky(cs)
    m = (sep.arcsec <= args.image_match_arcsec) & np.isfinite(z_img)
    m &= np.asarray(shapes["ellipticity"][idx] >= 0.3)
    pred = model.evaluate(images["ra"][m], images["dec"][m], z_img[m])
    img_off = lensmodel.axis_offset_deg(np.asarray(shapes["pa_obs"][idx[m]]), pred["tangential_pa"])

    def stats(rows):
        n = len(rows)
        a = int(np.sum(rows["orientation_class"] == "aligned"))
        b = int(np.sum(rows["orientation_class"] == "anti"))
        # For random orientations each class has probability 1/3 (30 deg out of 90).
        return {"n": n, "aligned": a, "anti": b, "p_aligned_excess": _binom_p(a, n, 1 / 3)}

    test = table[strong]
    summary = {
        "model": args.model,
        "catalog": str(args.catalog),
        "photoz": str(args.photoz) if args.photoz else None,
        "assumptions": {
            k: getattr(args, k)
            for k in (
                "min_ellipticity",
                "min_semimajor_px",
                "min_snr",
                "max_radius",
                "min_shear",
                "aligned_deg",
                "anti_deg",
                "z_margin",
            )
        },
        "multiple_images": {
            "n_matched_elongated": int(m.sum()),
            "median_offset_deg": float(np.median(img_off)) if m.any() else None,
            "n_within_30deg": int(np.sum(img_off <= 30)),
            "n_over_60deg": int(np.sum(img_off >= 60)),
        },
        "selected": len(table),
        "strong_shear": stats(test),
        "strong_shear_by_population": {
            p: stats(test[test["population"] == p])
            for p in ("background", "member_or_foreground", "unknown")
        },
    }
    out = args.out / args.model
    _write(table, out / "arcs.ecsv")
    anti = test[
        (test["orientation_class"] == "anti") & (test["population"] != "member_or_foreground")
    ]
    anti.sort("ellipticity", reverse=True)
    _write(anti, out / "arcs_anti.ecsv")
    summary["anti_candidates"] = [int(v) for v in anti["label"]]
    (out / "arcs.json").write_text(json.dumps(summary, indent=1))
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--model", choices=sorted(MODELS), default="smacs0723-iclv2")
    ap.add_argument("--out", type=Path, default=paths.outputs_dir() / "lens_consistency")
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate", help="compare with the published kappa map and multiple images")
    v.add_argument("--step", type=int, default=5, help="kappa-map sub-grid step in pixels")
    a = sub.add_parser("arcs", help="observed source orientation against the predicted shear")
    a.add_argument("--catalog", type=Path, required=True, help="JWST pipeline _cat.ecsv")
    a.add_argument("--photoz", type=Path, help="eazy zout FITS table (e.g. DJA)")
    a.add_argument("--min-ellipticity", type=float, default=0.5)
    a.add_argument("--min-semimajor-px", type=float, default=2.0)
    a.add_argument("--min-snr", type=float, default=10.0)
    a.add_argument("--max-radius", type=float, default=50.0, help="arcsec from the model reference")
    a.add_argument("--min-shear", type=float, default=0.2, help="reduced shear at z_s = 1")
    a.add_argument("--aligned-deg", type=float, default=30.0)
    a.add_argument("--anti-deg", type=float, default=60.0)
    a.add_argument("--z-margin", type=float, default=0.1)
    a.add_argument("--image-match-arcsec", type=float, default=0.5)
    args = ap.parse_args(argv)
    summary = cmd_validate(args) if args.cmd == "validate" else cmd_arcs(args)
    json.dump(summary, sys.stdout, indent=1)
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
