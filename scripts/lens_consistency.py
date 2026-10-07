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
    (``sky_orientation``, E of N) with the model's predicted stretching direction. Each source is
    tested at its own redshift range when an eazy ``zout`` table puts it behind the lens (z16,
    z_phot and z84); otherwise at z_s = 1, 2 and 4.

    * ``aligned``: within ``--aligned-deg`` of the prediction at every tested z_s;
    * ``anti``: at least ``--anti-deg`` away from it at every tested z_s;
    * ``mixed``: anything else.

    First, the known multiple images matched in the catalog check the orientation convention.
    Only background sources (photo-z 16th percentile above z_lens + ``--z-margin``) test the
    model, because cluster members and foreground galaxies are not lensed; the headline statistic
    uses them alone. ``anti`` background sources are candidates for inspection, not anomalies:
    intrinsic shapes, deblending fragments, photo-z failures and model error are the ordinary
    explanations to rule out.

Every threshold here is an ASSUMPTION. Inputs: the pinned model files (downloaded and verified by
sha256, ``SMACS0723_MAHLER22_ICLV2``), a level-3 ``_cat.ecsv`` and optionally a ``zout`` file.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import tarfile
import tempfile
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
    """Download (once) and verify the pinned files of model ``name``; extract its kappa map.

    The map member is copied out of the (sha256-verified) archive into a temporary file, checked
    against the member's size, and renamed into place, so an interrupted run leaves no partial map.
    """
    spec = MODELS[name]
    files = {key: fetch_catalog(url, sha) for key, (url, sha) in spec["files"].items()}
    archive = files.pop("kappa_0000")
    target = archive.with_name(
        archive.name.split(".", 1)[0] + "_" + Path(spec["kappa_member"]).name
    )
    with tarfile.open(archive) as tar:
        member = tar.getmember(spec["kappa_member"])
        if not target.exists() or target.stat().st_size != member.size:
            fd, tmp = tempfile.mkstemp(dir=target.parent, suffix=".part")
            try:
                with os.fdopen(fd, "wb") as out, tar.extractfile(member) as src:
                    while chunk := src.read(1 << 20):
                        out.write(chunk)
                if os.path.getsize(tmp) != member.size:
                    raise ValueError(f"{archive}: truncated member {spec['kappa_member']}")
                os.replace(tmp, target)
            finally:
                if os.path.exists(tmp):
                    os.unlink(tmp)
    files["kappa_map"] = target
    return files


def read_sigpos(path: Path) -> float:
    """``sigposArcsec`` from a Lenstool input file (only that keyword is read)."""
    match = re.search(r"^\s*sigposArcsec\s+(\S+)", path.read_text(encoding="latin-1"), re.M | re.I)
    if not match:
        raise SystemExit(f"error: {path} has no sigposArcsec; pass a model with one")
    return float(match.group(1))


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
    analytic = model.kappa_xy(x, y)
    diff = np.abs(analytic - published)
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
        "median_ratio": float(np.median(analytic / published)),
    }


def chi2pos_from_par(path: Path) -> float | None:
    """``#Chi2pos:`` written by Lenstool at the top of ``best.par`` (None if absent)."""
    match = re.search(r"^#Chi2pos:\s*(\S+)", path.read_text(encoding="latin-1"), re.M)
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
    flux = np.asarray(cat["isophotal_flux"], float)
    with np.errstate(divide="ignore", invalid="ignore"):
        snr = flux / np.asarray(cat["isophotal_flux_err"], float)
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
    z_margin: float = 0.1,
) -> Table:
    """Predicted stretching direction and its offset from the observed axis for each source.

    Columns per grid redshift (``pa_pred_z1``, ``offset_z1``, ``g_z1``, ``mu_z1`` ...) are always
    written for reference. The class uses each source's own redshift range (z160, z_phot, z840)
    when its photo-z puts it behind the lens (``z_basis = "photo-z"``), else the grid
    (``z_basis = "grid"``). ``g_min`` is the smallest reduced shear over the redshifts used.
    """
    out = sources.copy(copy_data=True)
    x, y = model.to_frame(out["ra"], out["dec"])
    fields = model.fields_xy(x, y)  # independent of z_s: computed once, scaled per redshift

    def predict(z):
        pred = model.evaluate(out["ra"], out["dec"], z, fields=fields)
        pa = np.asarray(pred["tangential_pa"], float)
        off = lensmodel.axis_offset_deg(out["pa_obs"], pa)
        return pa, off, np.asarray(pred["reduced_shear"], float), np.asarray(pred["magnification"])

    grid = []
    for z in z_grid:
        pa, off, g, mu = predict(z)
        tag = f"z{z:g}"
        out[f"pa_pred_{tag}"], out[f"offset_{tag}"], out[f"g_{tag}"], out[f"mu_{tag}"] = (
            pa,
            off,
            g,
            mu,
        )
        grid.append((pa, off, g))
    n = len(out)
    use_pz = np.zeros(n, bool)
    if {"z_phot", "z160", "z840"} <= set(out.colnames):
        zs = [np.asarray(out[c], float) for c in ("z160", "z_phot", "z840")]
        use_pz = np.isfinite(zs[0]) & (zs[0] > model.z_lens + z_margin)
        own = []
        for z in zs:
            pa, off, g, _ = predict(np.where(use_pz, z, np.nan))
            own.append((pa, off, g))
        out["offset_zphot"] = np.where(use_pz, own[1][1], np.nan)
    else:
        own = grid
    sets = [
        (
            np.where(use_pz, o[0], gr[0]),
            np.where(use_pz, o[1], gr[1]),
            np.where(use_pz, o[2], gr[2]),
        )
        for o, gr in zip(own, grid, strict=True)
    ]
    offsets = np.vstack([s[1] for s in sets])
    out["z_basis"] = np.where(use_pz, "photo-z", "grid")
    out["offset_min"] = offsets.min(axis=0)
    out["offset_max"] = offsets.max(axis=0)
    out["g_min"] = np.vstack([s[2] for s in sets]).min(axis=0)
    spread = np.zeros(n)
    for a in range(len(sets)):
        for b in range(a + 1, len(sets)):
            spread = np.maximum(spread, lensmodel.axis_offset_deg(sets[a][0], sets[b][0]))
    out["pa_pred_spread"] = spread  # > 45 deg: the prediction flips near a critical curve
    cls = np.where(out["offset_max"] <= aligned_deg, "aligned", "mixed")
    out["orientation_class"] = np.where(out["offset_min"] >= anti_deg, "anti", cls)
    out.meta.update(model._meta())
    out.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source=f"{sources.meta.get('source', 'catalog')} against {model.source}",
        z_grid=list(z_grid),
        assumptions={"aligned_deg": aligned_deg, "anti_deg": anti_deg, "z_margin": z_margin},
    )
    return out


def _binom_p(k: int, n: int, p: float) -> float:
    """One-sided binomial P(X >= k)."""
    from scipy.stats import binomtest

    return float(binomtest(k, n, p, alternative="greater").pvalue) if n else float("nan")


def class_stats(rows: Table, aligned_deg: float) -> dict:
    """Counts per class and the chance of that many ``aligned`` under random orientations.

    A random axis falls within ``aligned_deg`` of a given direction with probability
    ``aligned_deg / 90``; requiring it at every tested redshift makes this an upper bound.
    """
    n = len(rows)
    a = int(np.sum(rows["orientation_class"] == "aligned"))
    b = int(np.sum(rows["orientation_class"] == "anti"))
    p0 = aligned_deg / 90.0
    return {"n": n, "aligned": a, "anti": b, "p_random": p0, "p_aligned_excess": _binom_p(a, n, p0)}


def _write(table: Table, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    table.write(path, overwrite=True)


def cmd_validate(args) -> dict:
    files = model_files(args.model)
    par = lensmodel.parse_lenstool_par(files["best.par"])
    model = lensmodel.LensModel.from_par(par)
    sigpos = read_sigpos(files["input.par"])
    images = lensmodel.load_lenstool_images(files["arcs.dat"])
    chi2_ref = chi2pos_from_par(files["best.par"])
    bt, bsum = backtrace_check(model, images, par["z_m_limit"], sigpos, chi2_ref)
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


def convention_check(model, par, images, shapes, match_arcsec: float) -> dict:
    """Offsets of catalogued multiple images (e >= 0.3) from the prediction at their model z."""
    z_img = lensmodel.image_redshifts(images, par["z_m_limit"])
    ci = SkyCoord(images["ra"], images["dec"], unit="deg")
    cs = SkyCoord(shapes["ra"], shapes["dec"], unit="deg")
    idx, sep, _ = ci.match_to_catalog_sky(cs)
    m = (sep.arcsec <= match_arcsec) & np.isfinite(z_img) & (shapes["ellipticity"][idx] >= 0.3)
    if not m.any():
        return {"n_matched_elongated": 0}
    pred = model.evaluate(images["ra"][m], images["dec"][m], z_img[m])
    off = lensmodel.axis_offset_deg(np.asarray(shapes["pa_obs"][idx[m]]), pred["tangential_pa"])
    strong = np.asarray(shapes["ellipticity"][idx[m]]) >= 0.6
    return {
        "n_matched_elongated": int(m.sum()),
        "median_offset_deg": float(np.median(off)),
        "n_within_30deg": int(np.sum(off <= 30)),
        "n_over_60deg": int(np.sum(off >= 60)),
        "e_ge_0.6": {
            "n": int(strong.sum()),
            "max_offset_deg": float(off[strong].max()) if strong.any() else None,
        },
        "over_60deg": [str(i) for i in images["image_id"][m][off >= 60]],
    }


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
        model,
        shapes[sel],
        aligned_deg=args.aligned_deg,
        anti_deg=args.anti_deg,
        z_margin=args.z_margin,
    )
    table["strong_shear"] = np.asarray(table["g_min"], float) >= args.min_shear
    if "z160" in table.colnames:
        zl = model.z_lens + args.z_margin
        pop = np.where(table["z160"] > zl, "background", "unknown")
        pop = np.where(table["z840"] < zl, "member_or_foreground", pop)
    else:
        pop = np.full(len(table), "unknown")
    table["population"] = pop

    images = lensmodel.load_lenstool_images(files["arcs.dat"])
    test = table[table["strong_shear"]]
    by_pop = {
        p: class_stats(test[test["population"] == p], args.aligned_deg)
        for p in ("background", "member_or_foreground", "unknown")
    }
    keys = (
        "min_ellipticity",
        "min_semimajor_px",
        "min_snr",
        "max_radius",
        "min_shear",
        "aligned_deg",
        "anti_deg",
        "z_margin",
    )
    summary = {
        "model": args.model,
        "catalog": str(args.catalog),
        "photoz": str(args.photoz) if args.photoz else None,
        "assumptions": {k: getattr(args, k) for k in keys},
        "multiple_images": convention_check(model, par, images, shapes, args.image_match_arcsec),
        "selected": len(table),
        "strong_shear_background": by_pop["background"],
        "strong_shear_by_population": by_pop,
        "strong_shear_all": class_stats(test, args.aligned_deg),
    }
    out = args.out / args.model
    _write(table, out / "arcs.ecsv")
    keep = (test["orientation_class"] == "anti") & (test["population"] != "member_or_foreground")
    anti = test[keep]
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
    a.add_argument("--min-shear", type=float, default=0.2, help="min reduced shear over tested z")
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
