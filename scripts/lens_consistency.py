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

``images``
    Forward-predicts every image of each catalogued multiple-image system: the system's mean
    back-traced source position is solved for all images (``lensmodel.find_images``) and each
    predicted image is compared with ``arcs.dat`` and a pipeline catalog. A predicted image's
    magnitude is the reference image's catalog magnitude scaled by the magnification ratio (the
    reference is the catalogued image with the smallest predicted |μ|, i.e. the least sensitive to
    the critical curves). Classes: ``observed`` (an ``arcs.dat`` image within ``--match-arcsec``),
    ``candidate`` (an uncatalogued source there whose photo-z allows the system redshift;
    without ``--photoz`` every nearby source qualifies), ``other_source`` (a source there whose
    photo-z excludes it),
    ``demagnified`` (|μ| < 0.5, e.g. central images), ``missing`` (predicted brighter than the
    depth, nothing there), ``faint`` (predicted below the depth), ``no_flux_ref`` (no catalogued
    image of the system matched in the catalog, so no predicted magnitude) and ``outside`` (no
    catalog source within ``--footprint-arcsec``: off the image or inside a bright galaxy's
    segment).
    ``arcs.dat`` images that no predicted image reproduces are reported as ``unpredicted``.

    Pipeline segments of arcs near cluster galaxies are unreliable flux references (SMACS: the
    catalogued images sit 0.7–2.6″ from their nearest segment). With ``--forced-image`` (an
    ``_i2d`` URI or path; S3 is read by byte range), every non-observed predicted image is also
    checked by forced aperture photometry: the reference flux is measured at the system's
    catalogued images (recentred on the peak within 0.4″), scaled by |μ| to a predicted flux and
    S/N, and compared with the best aperture within ``--forced-search-arcsec``. Forced classes:
    ``undetectable`` (predicted < 5σ), ``recovered`` (a ≥ 5σ source there at 1/3–3× the predicted
    flux), ``confused`` (a ≥ 5σ source more than 3× brighter dominates, e.g. a cluster galaxy),
    ``absent`` (predicted ≥ 10σ; best < 3σ or < 1/3 of the predicted flux), ``ambiguous``,
    ``no_reference`` and ``off_image`` (< 80 % valid pixels). Offsets are West/North arcsec.
    ERR is scaled by 1.5 (D-027 amendment).

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
            "mag": (
                np.asarray(cat["isophotal_abmag"], float)
                if "isophotal_abmag" in cat.colnames
                else np.full(len(cat), np.nan)
            ),
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


def depth_mag(shapes: Table, snr_range=(4.0, 6.0)) -> float:
    """Median magnitude of catalog sources at S/N 4-6: a 5-sigma depth proxy (ASSUMPTION)."""
    m = (
        (shapes["snr"] >= snr_range[0])
        & (shapes["snr"] <= snr_range[1])
        & np.isfinite(shapes["mag"])
    )
    return float(np.median(shapes["mag"][m])) if m.any() else float("nan")


def predict_counter_images(
    model,
    grid,
    backtrace: Table,
    shapes: Table,
    match_arcsec: float = 1.5,
    footprint_arcsec: float = 5.0,
    ref_match_arcsec: float = 0.5,
    depth: float = float("nan"),
    z_margin: float = 0.1,
    min_abs_mu: float = 0.5,
) -> tuple[Table, list[str]]:
    """Predicted images of every catalogued system, classified against the observations."""
    cs = SkyCoord(shapes["ra"], shapes["dec"], unit="deg")
    rows, unpredicted = [], []
    systems = np.asarray(backtrace["system"]).astype(str)
    for sys_id in dict.fromkeys(systems):
        sel = (systems == sys_id) & np.isfinite(backtrace["beta_x"])
        if sel.sum() < 2:
            continue
        obs = backtrace[sel]
        z = float(obs["z_used"][0])
        pred = lensmodel.find_images(
            model, grid, float(np.mean(obs["beta_x"])), float(np.mean(obs["beta_y"])), z
        )
        co = SkyCoord(obs["ra"], obs["dec"], unit="deg")
        cp = SkyCoord(pred["ra"], pred["dec"], unit="deg")
        # reference image: catalogued, matched in the catalog, smallest predicted |mu|
        idx, sep, _ = co.match_to_catalog_sky(cs)
        ok = (
            (sep.arcsec <= ref_match_arcsec)
            & np.isfinite(np.asarray(shapes["mag"])[idx])
            & np.isfinite(np.asarray(obs["magnification"], float))
        )
        m_ref = mu_ref = np.nan
        ref_id = ""
        if ok.any():
            k = np.where(ok)[0][np.argmin(np.abs(np.asarray(obs["magnification"])[ok]))]
            m_ref, mu_ref, ref_id = (
                float(shapes["mag"][idx[k]]),
                float(obs["magnification"][k]),
                str(obs["image_id"][k]),
            )
        for img_id, c in zip(obs["image_id"], co, strict=True):
            if not len(pred) or c.separation(cp).arcsec.min() > match_arcsec:
                unpredicted.append(str(img_id))
        for p, c in zip(pred, cp, strict=True):
            d_obs = c.separation(co).arcsec
            j = int(np.argmin(d_obs))
            d_cat = c.separation(cs).arcsec
            q = int(np.argmin(d_cat))
            mu = float(p["magnification"])
            m_pred = m_ref - 2.5 * np.log10(abs(mu) / abs(mu_ref)) if np.isfinite(m_ref) else np.nan
            near = d_cat[q] <= match_arcsec
            zlo = float(shapes["z160"][q]) if "z160" in shapes.colnames else np.nan
            zhi = float(shapes["z840"][q]) if "z840" in shapes.colnames else np.nan
            z_ok = not (np.isfinite(zlo) and np.isfinite(zhi)) or (
                zlo - z_margin <= z <= zhi + z_margin
            )
            if d_obs[j] <= match_arcsec:
                cls = "observed"
            elif abs(mu) < min_abs_mu:
                cls = "demagnified"  # e.g. a central image: expected undetectable
            elif near and z_ok:
                cls = "candidate"
            elif near:
                cls = "other_source"  # a catalog source there, but its photo-z excludes z_sys
            elif d_cat[q] > footprint_arcsec:
                cls = "outside"
            elif not np.isfinite(m_pred):
                cls = "no_flux_ref"  # no catalogued image of the system matched in the catalog
            elif np.isfinite(depth) and m_pred < depth:
                cls = "missing"
            else:
                cls = "faint"
            rows.append(
                {
                    "system": sys_id,
                    "z_sys": z,
                    "x": float(p["x"]),
                    "y": float(p["y"]),
                    "ra": float(p["ra"]),
                    "dec": float(p["dec"]),
                    "magnification": mu,
                    "image_class": cls,
                    "matched_image": str(obs["image_id"][j]) if d_obs[j] <= match_arcsec else "",
                    "sep_image_arcsec": float(d_obs[j]),
                    "ref_image": ref_id,
                    "mag_pred": float(m_pred),
                    "nearest_label": int(shapes["label"][q]),
                    "sep_catalog_arcsec": float(d_cat[q]),
                    "nearest_mag": float(shapes["mag"][q]),
                    "nearest_z160": zlo,
                    "nearest_z840": zhi,
                }
            )
    out = Table(rows=rows) if rows else Table()
    out.meta.update(model._meta())
    out.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        depth_mag=depth,
        assumptions={
            "match_arcsec": match_arcsec,
            "footprint_arcsec": footprint_arcsec,
            "z_margin": z_margin,
        },
    )
    return out, unpredicted


FORCED_R_AP = 0.2  # aperture radius, arcsec (ASSUMPTION)
FORCED_ANNULUS = (0.6, 1.0)  # background annulus, arcsec (ASSUMPTION)
ERR_SCALE = 1.5  # ERR underestimates the noise by 1.2-1.5x (D-027 amendment)
MIN_VALID = 0.8  # minimum finite fraction of aperture and annulus pixels (ASSUMPTION)


def aperture_snr(img, err, xx, yy, cx: float, cy: float) -> tuple[float, float]:
    """Background-subtracted flux and error in a ``FORCED_R_AP`` aperture at offset (cx, cy).

    ``xx``/``yy`` are each pixel's offset from the target in arcsec (West, North). NaN when fewer
    than ``MIN_VALID`` of the aperture or annulus pixels are finite (gaps, edges, masks)."""
    r = np.hypot(xx - cx, yy - cy)
    ap = r <= FORCED_R_AP
    ann = (r > FORCED_ANNULUS[0]) & (r < FORCED_ANNULUS[1])
    good = np.isfinite(img) & np.isfinite(err)
    if not ap.any() or good[ap].mean() < MIN_VALID or good[ann].mean() < MIN_VALID:
        return float("nan"), float("nan")
    bkg = float(np.nanmedian(img[ann & good]))
    flux = float(np.nansum(img[ap] - bkg))
    ferr = float(np.sqrt(np.nansum(err[ap] ** 2))) * ERR_SCALE
    return flux, ferr


def peak_flux(img, err, xx, yy, radius: float = 0.4) -> tuple[float, float]:
    """Aperture flux at the brightest pixel within ``radius`` of the stamp centre."""
    inner = np.where((np.hypot(xx, yy) <= radius) & np.isfinite(img), img, -np.inf)
    if not np.isfinite(inner).any():
        return float("nan"), float("nan")
    k = int(np.argmax(inner))
    return aperture_snr(img, err, xx, yy, float(xx.flat[k]), float(yy.flat[k]))


def best_within(
    img, err, xx, yy, search: float, step: float = 0.1
) -> tuple[float, float, float, float]:
    """Highest aperture S/N on a grid of centres within ``search``: (snr, flux, dx, dy)."""
    best = (-np.inf, np.nan, 0.0, 0.0)
    for cx in np.arange(-search, search + 1e-9, step):
        for cy in np.arange(-search, search + 1e-9, step):
            if np.hypot(cx, cy) > search:
                continue
            f, e = aperture_snr(img, err, xx, yy, cx, cy)
            if np.isfinite(f) and e > 0 and f / e > best[0]:
                best = (f / e, f, float(cx), float(cy))
    return best


MAX_FLUX_RATIO = 3.0  # best/predicted flux above this: a brighter source dominates (ASSUMPTION)
MIN_FLUX_RATIO = 1 / 3  # below this the best source is too faint to be the image (ASSUMPTION)


def forced_class(pred_snr: float, best_snr: float, flux_ratio: float = 1.0) -> str:
    """Forced-photometry verdict for one predicted image (thresholds are ASSUMPTIONs).

    ``flux_ratio`` is the best aperture's flux over the predicted flux."""
    if not np.isfinite(pred_snr):
        return "no_reference"
    if pred_snr < 5:
        return "undetectable"
    if best_snr >= 5 and flux_ratio > MAX_FLUX_RATIO:
        return "confused"
    if best_snr >= 5 and flux_ratio >= MIN_FLUX_RATIO:
        return "recovered"
    if pred_snr >= 10 and (best_snr < 3 or flux_ratio < MIN_FLUX_RATIO):
        return "absent"  # nothing there, or only a source far fainter than predicted
    return "ambiguous"


def forced_check(
    table: Table,
    backtrace: Table,
    stamp,
    search_arcsec: float = 1.0,
    max_ref_mu: float = 50.0,
) -> None:
    """Add forced-photometry columns to the predicted-image ``table`` in place.

    ``stamp(ra, dec)`` returns ``(sci, err, xx, yy)`` around a position. The reference is the
    system's catalogued image at S/N > 5 with the smallest |μ| below ``max_ref_mu`` (μ near a
    critical curve is too uncertain to scale from; ASSUMPTION)."""
    n = len(table)
    keys = ("pred_snr", "best_snr", "flux_ratio", "best_dx", "best_dy")
    cols = {k: np.full(n, np.nan) for k in keys}
    fclass = np.full(n, "", dtype="U12")
    ref_cache: dict[str, tuple[float, float]] = {}
    systems = np.asarray(backtrace["system"]).astype(str)
    for i, row in enumerate(table):
        if row["image_class"] in ("observed", "demagnified", "outside"):
            continue
        sys_id = str(row["system"])
        if sys_id not in ref_cache:
            ref = (np.nan, np.nan)
            for b in backtrace[systems == sys_id]:
                mu = abs(float(b["magnification"]))
                if not np.isfinite(mu) or mu > max_ref_mu:
                    continue
                f, e = peak_flux(*stamp(float(b["ra"]), float(b["dec"])))
                if (
                    np.isfinite(f)
                    and e > 0
                    and f / e > 5
                    and (not np.isfinite(ref[1]) or mu < ref[1])
                ):
                    ref = (f, mu)
            ref_cache[sys_id] = ref
        f_ref, mu_ref = ref_cache[sys_id]
        if not np.isfinite(f_ref):
            fclass[i] = "no_reference"
            continue
        sci, err, xx, yy = stamp(float(row["ra"]), float(row["dec"]))
        _, e0 = aperture_snr(sci, err, xx, yy, 0.0, 0.0)
        if not (np.isfinite(e0) and e0 > 0):
            fclass[i] = "off_image"  # off the footprint, or in a gap or masked region
            continue
        pred = f_ref * abs(float(row["magnification"])) / mu_ref
        cols["pred_snr"][i] = pred / e0 if e0 > 0 else np.nan
        snr, flux, dx, dy = best_within(sci, err, xx, yy, search_arcsec)
        cols["best_snr"][i], cols["best_dx"][i], cols["best_dy"][i] = snr, dx, dy
        cols["flux_ratio"][i] = flux / pred if pred > 0 else np.nan
        fclass[i] = forced_class(cols["pred_snr"][i], snr, cols["flux_ratio"][i])
    for k, v in cols.items():
        table[k] = v
    table["forced_class"] = fclass
    table.meta["forced"] = {
        "r_ap_arcsec": FORCED_R_AP,
        "annulus_arcsec": list(FORCED_ANNULUS),
        "err_scale": ERR_SCALE,
        "search_arcsec": search_arcsec,
        "max_ref_mu": max_ref_mu,
        "flux_ratio_range": [MIN_FLUX_RATIO, MAX_FLUX_RATIO],
        "min_valid_fraction": MIN_VALID,
    }


def image_stamper(uri: str, half_arcsec: float = 1.5):
    """``stamp(ra, dec)`` reading SCI/ERR sections of an ``_i2d`` (local path or S3 by byte range).

    Returns ``(stamper, close)``."""
    import fsspec

    opts = {"anon": True} if uri.startswith("s3://") else {}
    fs, path = fsspec.core.url_to_fs(uri, **opts)
    fo = fs.open(path, "rb", block_size=2**20, cache_type="readahead")
    hdul = fits.open(fo, lazy_load_hdus=True)
    sci, err = hdul["SCI"], hdul["ERR"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        wcs = WCS(sci.header)
    pix = float(abs(wcs.proj_plane_pixel_scales()[0].to_value("arcsec")))
    hp = int(np.ceil(half_arcsec / pix))
    ny, nx = sci.header["NAXIS2"], sci.header["NAXIS1"]
    jj, ii = np.mgrid[-hp : hp + 1, -hp : hp + 1].astype(float)
    # pixel offsets -> (dRA cos dec, dDec) in arcsec; West = -dRA cos dec, as in the model frame
    cd = wcs.pixel_scale_matrix * 3600.0

    def stamp(ra: float, dec: float):
        xf, yf = (float(v) for v in wcs.world_to_pixel_values(ra, dec))
        x, y = int(round(xf)), int(round(yf))
        di, dj = ii + (x - xf), jj + (y - yf)  # offsets from the exact target position
        xx = -(cd[0, 0] * di + cd[0, 1] * dj)
        yy = cd[1, 0] * di + cd[1, 1] * dj
        if not (hp <= x < nx - hp and hp <= y < ny - hp):
            nan = np.full(xx.shape, np.nan)
            return nan, nan, xx, yy
        sl = (slice(y - hp, y + hp + 1), slice(x - hp, x + hp + 1))
        return np.asarray(sci.section[sl], float), np.asarray(err.section[sl], float), xx, yy

    def close() -> None:
        hdul.close()
        fo.close()

    return stamp, close


def cmd_images(args) -> dict:
    files = model_files(args.model)
    par = lensmodel.parse_lenstool_par(files["best.par"])
    model = lensmodel.LensModel.from_par(par)
    cache = (
        paths.cache_dir()
        / "external"
        / f"{model.sha256[:12]}_alpha_{args.half_width:g}_{args.step:g}.npz"
    )
    grid = lensmodel.DeflectionGrid.cached(model, cache, args.half_width, args.step)
    images = lensmodel.load_lenstool_images(files["arcs.dat"])
    bt = lensmodel.backtrace_images(model, images, par["z_m_limit"])
    shapes = load_shapes(args.catalog)
    if args.photoz:
        attach_photoz(shapes, args.photoz)
    depth = depth_mag(shapes)
    table, unpredicted = predict_counter_images(
        model, grid, bt, shapes, args.match_arcsec, args.footprint_arcsec, depth=depth
    )
    if not len(table):
        raise SystemExit("no system has two or more back-traced images; nothing to predict")
    if args.forced_image:
        stamp, close = image_stamper(
            args.forced_image, args.forced_search_arcsec + FORCED_ANNULUS[1] + 0.1
        )
        try:
            forced_check(table, bt, stamp, args.forced_search_arcsec)
        finally:
            close()
        table.meta["forced"]["image"] = args.forced_image
    out = args.out / args.model
    _write(table, out / "images_predicted.ecsv")
    classes = (
        "observed",
        "demagnified",
        "candidate",
        "other_source",
        "missing",
        "faint",
        "no_flux_ref",
        "outside",
    )
    counts = {c: int(np.sum(table["image_class"] == c)) for c in classes}
    flagged = table[np.isin(table["image_class"], ["candidate", "missing"])]
    summary = {
        "model": args.model,
        "catalog": str(args.catalog),
        "depth_mag_5sigma_proxy": depth,
        "n_systems": len(set(table["system"])),
        "n_predicted": len(table),
        "classes": counts,
        "unpredicted_catalogued_images": unpredicted,
        "forced": (
            {
                **table.meta["forced"],
                "classes": {
                    c: int(np.sum(table["forced_class"] == c))
                    for c in (
                        "recovered",
                        "confused",
                        "absent",
                        "undetectable",
                        "ambiguous",
                        "no_reference",
                        "off_image",
                    )
                },
                "rows": [
                    {
                        "system": str(r["system"]),
                        "catalog_class": str(r["image_class"]),
                        "forced_class": str(r["forced_class"]),
                        "xy": [round(float(r["x"]), 2), round(float(r["y"]), 2)],
                        "mu": round(float(r["magnification"]), 2),
                        "pred_snr": round(float(r["pred_snr"]), 1),
                        "best_snr": round(float(r["best_snr"]), 1),
                        "flux_ratio": round(float(r["flux_ratio"]), 2),
                        "best_offset_arcsec": [
                            round(float(r["best_dx"]), 1),
                            round(float(r["best_dy"]), 1),
                        ],
                    }
                    for r in table
                    if r["forced_class"]
                ],
            }
            if "forced_class" in table.colnames
            else None
        ),
        "flagged": [
            {
                "system": str(r["system"]),
                "class": str(r["image_class"]),
                "xy": [round(float(r["x"]), 2), round(float(r["y"]), 2)],
                "mu": round(float(r["magnification"]), 2),
                "mag_pred": round(float(r["mag_pred"]), 2),
                "nearest_label": int(r["nearest_label"]),
                "sep_catalog_arcsec": round(float(r["sep_catalog_arcsec"]), 2),
            }
            for r in flagged
        ],
    }
    (out / "images.json").write_text(json.dumps(summary, indent=1))
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
    i = sub.add_parser("images", help="forward-predict every image of the catalogued systems")
    i.add_argument("--catalog", type=Path, required=True, help="JWST pipeline _cat.ecsv")
    i.add_argument("--photoz", type=Path, help="eazy zout FITS table (e.g. DJA)")
    i.add_argument("--match-arcsec", type=float, default=1.5)
    i.add_argument("--footprint-arcsec", type=float, default=5.0)
    i.add_argument("--half-width", type=float, default=60.0, help="solver grid half-width, arcsec")
    i.add_argument("--step", type=float, default=0.1, help="solver grid step, arcsec")
    i.add_argument(
        "--forced-image", help="_i2d URI or path for forced photometry (S3 by byte range)"
    )
    i.add_argument("--forced-search-arcsec", type=float, default=1.0)
    args = ap.parse_args(argv)
    summary = {"validate": cmd_validate, "arcs": cmd_arcs, "images": cmd_images}[args.cmd](args)
    json.dump(summary, sys.stdout, indent=1)
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
