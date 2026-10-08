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

``fluxratios``
    Compares each catalogued system's image fluxes and colours with the model (lensing is
    achromatic and scales flux by |μ|). Each ``arcs.dat`` image is matched to the nearest DJA
    ``fix_phot`` source within ``--match-arcsec``. The flux is the detection-image Kron total
    ``mag_auto``: DJA aperture fluxes are not totals for extended arcs (``<band>_tot_corr`` is 1
    and ``tot_corr`` a point-source correction). The implied source magnitude is
    ``s = mag_auto + 2.5 log10 |μ|``, with μ at the catalogued position and system redshift;
    an image's residual is ``s`` minus the median ``s`` of the *other* usable images of its
    system (leave-one-out; two or more usable images needed), and ``flux_outlier`` marks the
    system's largest |residual| (both images of a pair) when it exceeds ``--outlier-mag``.
    The colour test does the same with the ``--colour`` aperture colour (``colour_outlier``
    above ``--colour-tol``; it does not use μ). Images are unusable when unmatched, at S/N <
    ``--min-snr``, when one DJA source matches several catalogued images, or when the DJA
    segment exceeds ``--max-npix`` pixels (it swallows host or ICL light; SMACS 1.1), or, with
    ``--photoz``, when the counterpart's 95 % photo-z interval, widened by ``--z-margin`` ×
    (1 + z) because eazy intervals are too narrow (SMACS 3.3: 1.83–1.86 at z = 1.99), excludes the
    system redshift (another object at the position; El Gordo 9a). The flux
    test also drops |μ| > ``--max-abs-mu`` (near the critical curves the model's μ is
    unreliable). Ordinary explanations of an outlier: blending, segmentation, microlensing,
    substructure the model omits, time delay with variability. The model parity (sign of μ) is
    reported but not tested.

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

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS, FITSFixedWarning

from jwst_anomaly import lensmodel, paths, schema
from jwst_anomaly.features import snr_from_mag_err
from jwst_anomaly.photometry import fetch_catalog, load_dja_catalog

# sigpos: "input.par" (its sigposArcsec), "arcs" (each image's error column, as Lenstool's
# ``forme -10``) or a number. El Gordo's CDS best_fit.par has an empty image section: a uniform
# 0.621" (the file's smallest error) reproduces Lenstool's Chi2pos best (82.5 vs 80.22; the
# per-image column gives 52.0), so it is an ASSUMPTION backed by that match (D-030).
MODELS = {
    "smacs0723-iclv2": {
        "files": lensmodel.SMACS0723_MAHLER22_ICLV2,
        "kappa_member": "tmp_k/0000_k.fits",
        "sigpos": "input.par",
    },
    # frame_offset_arcsec: (dRA cos dec, dDec) from the model's image frame to the JWST frame;
    # El Gordo's image list is on the HST/RELICS frame (median offset to DJA v7.0, issue #41).
    "elgordo-caminha23": {
        "files": lensmodel.ELGORDO_CAMINHA23,
        "sigpos": 0.621,
        "frame_offset_arcsec": (0.224, -0.016),
    },
    "abell2744-bergamini23": {"files": lensmodel.ABELL2744_BERGAMINI23, "sigpos": "arcs"},
    # CANUCS JWST-era Lenstool models (D-044): the independent second model for vetting. The
    # Abell 370 fit is source-plane (image-plane rms 2.3"), so its image list is gated off.
    "macs0416-canucs": {"files": lensmodel.MACS0416_CANUCS, "sigpos": "input.par"},
    "abell370-canucs": {
        "files": lensmodel.ABELL370_CANUCS,
        "sigpos": "input.par",
        "frame_offset_arcsec": (-0.148, 0.002),
        "image_list_ok": False,
    },
    # map models: published deflection maps (D_LS/D_S = 1), no Lenstool par or image list
    "whl0137-relics-lenstool": {
        "files": lensmodel.WHL0137_RELICS_LENSTOOL,
        "kind": "maps",
        "z_lens": 0.566,
        "cosmology": (70.0, 0.3),
        "mag_maps": {6.2: "mag_z6.2"},
    },
}

# HFF CATS map models (D-035): "<cluster>-cats". The image list is used only where our
# image-plane rms is within IMAGE_LIST_RMS_FACTOR x the release's quoted rms (ASSUMPTION);
# ``validate`` re-measures it and reports whether its measurement agrees with this setting.
# Abell 370, Abell S1063 and Abell 2744 have no params.txt (no fitted image redshifts).
IMAGE_LIST_RMS_FACTOR = 1.5
_HFF_QUOTED_RMS = {"macs0416": 0.72, "macs1149": 0.63, "macs0717": 2.41, "abells1063": 0.48}
# macs0416 passes once find_images refines cells on folds (D-040)
_HFF_IMAGE_LIST_OK = {"macs1149", "macs0717", "macs0416"}
# (dRA cos dec, dDec) from arcs.txt to the JWST frame where it exceeds 0.1" (D-034, D-038)
_HFF_FRAME_OFFSET = {"macs0416": (0.208, -0.025), "abell370": (-0.121, -0.015)}
for _c, (_v, _zl, _files) in lensmodel.HFF_CATS.items():
    MODELS[f"{_c}-cats"] = {
        "files": _files,
        "kind": "maps",
        "z_lens": _zl,
        "cosmology": (70.0, 0.3),
        "mag_maps": {2.0: "mag_z2"},
        "version": _v,
        "image_list_ok": _c in _HFF_IMAGE_LIST_OK,
        "quoted_rms_arcsec": _HFF_QUOTED_RMS.get(_c),
    }
    if _c in _HFF_FRAME_OFFSET:
        MODELS[f"{_c}-cats"]["frame_offset_arcsec"] = _HFF_FRAME_OFFSET[_c]


def is_map_model(name: str) -> bool:
    return MODELS[name].get("kind") == "maps"


def has_image_list(name: str) -> bool:
    """Whether ``images`` and ``fluxratio`` may use the model's multiple-image list."""
    spec = MODELS[name]
    return "arcs.dat" in spec["files"] and spec.get("image_list_ok", True)


def image_list(name: str, files: dict[str, Path], par) -> tuple[Table, dict[str, float]]:
    """The model's catalogued multiple images and its fixed system redshifts (``z_m_limit``;
    none for a map model, whose image list carries the redshifts it used, 0 = unknown)."""
    if "arcs.dat" not in files:
        files.update(model_files(name, ("arcs.dat",)))
    images = lensmodel.load_lenstool_images(files["arcs.dat"])
    if par is not None:
        return images, par["z_m_limit"]
    if "params.txt" in MODELS[name]["files"]:
        if "params.txt" not in files:
            files.update(model_files(name, ("params.txt",)))
        return images, lensmodel.read_z_m_limit(files["params.txt"])
    return images, {}


def apply_frame_offset(name: str, model, images: Table | None = None) -> tuple[float, float]:
    """Shift a model and its image list from the image list's frame to the JWST frame by
    ``MODELS[name]["frame_offset_arcsec"]`` (dRA cos dec, dDec; D-034), in place. The model's
    reference point moves with the images, so every model position lands in the JWST frame. A map
    model's WCS moves too: its maps are looked up by sky position, so moving only the reference
    point would leave them in place (D-040)."""
    dra, ddec = MODELS[name].get("frame_offset_arcsec", (0.0, 0.0))
    if dra or ddec:
        model.shift_frame(dra, ddec)
        if images is not None:
            shift_images(name, images)
    return float(dra), float(ddec)


def shift_images(name: str, images: Table) -> Table:
    """Shift an image list into the JWST frame by the model's ``frame_offset_arcsec``, in place."""
    dra, ddec = MODELS[name].get("frame_offset_arcsec", (0.0, 0.0))
    if dra or ddec:
        cos = np.cos(np.deg2rad(np.asarray(images["dec"], float)))
        images["ra"] = np.asarray(images["ra"], float) + dra / 3600.0 / cos
        images["dec"] = np.asarray(images["dec"], float) + ddec / 3600.0
    return images


def load_model(name: str):
    """``(model, files, par)`` for ``MODELS[name]``; ``par`` is None for a map model."""
    from astropy.cosmology import FlatLambdaCDM

    if is_map_model(name):
        spec = MODELS[name]
        files = model_files(name, ("alpha_x", "alpha_y"))  # the check maps only for validate
        h0, om = spec["cosmology"]
        model = lensmodel.MapLensModel.from_fits(
            files["alpha_x"],
            files["alpha_y"],
            spec["z_lens"],
            FlatLambdaCDM(H0=h0, Om0=om),
            source=name,
            centre=spec.get("centre"),
        )
        return model, files, None
    files = model_files(name)
    par = lensmodel.parse_lenstool_par(files["best.par"])
    return lensmodel.LensModel.from_par(par), files, par


Z_GRID = (1.0, 2.0, 4.0)


def model_files(name: str, keys=None) -> dict[str, Path]:
    """Download (once) and verify the pinned files of model ``name``; extract its kappa map.

    The map member is copied out of the (sha256-verified) archive into a temporary file, checked
    against the member's size, and renamed into place, so an interrupted run leaves no partial map.
    """
    spec = MODELS[name]
    files = {
        key: fetch_catalog(url, sha)
        for key, (url, sha) in spec["files"].items()
        if keys is None or key in keys
    }
    if "kappa_member" not in spec:
        return files
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


def model_sigpos(name: str, files: dict[str, Path], images: Table) -> np.ndarray:
    """Per-image position error (arcsec) under the model's ``sigpos`` rule (see ``MODELS``)."""
    rule = MODELS[name]["sigpos"]
    if rule == "input.par":
        return np.full(len(images), read_sigpos(files["input.par"]))
    sigma = np.asarray(images["a"], float) if rule == "arcs" else np.full(len(images), float(rule))
    if not np.all(np.isfinite(sigma) & (sigma > 0)):
        raise SystemExit(f"error: {name}: position errors must be finite and positive")
    return sigma


def grid_cache_path(model: lensmodel.LensModel, half_width: float, step: float) -> Path:
    """Cache file of a model's deflection grid (shared by ``validate`` and ``images``)."""
    return paths.cache_dir() / "external" / f"{model.sha256[:12]}_alpha_{half_width:g}_{step:g}.npz"


def grid_half_width(model: lensmodel.LensModel, images: Table, margin: float = 20.0) -> float:
    """Solver grid half-width (arcsec) that covers every catalogued image plus ``margin``."""
    x, y = model.to_frame(np.asarray(images["ra"], float), np.asarray(images["dec"], float))
    return float(np.ceil(max(np.abs(x).max(), np.abs(y).max()) + margin))


def imageplane_check(model, grid, bt: Table, sigma: np.ndarray, chi2_ref: float | None):
    """Exact image-plane residuals plus a summary (χ², rms, > 3σ images)."""
    ip = lensmodel.imageplane_residuals(model, grid, bt)
    ip["sigma_arcsec"] = sigma
    d = np.asarray(ip["dtheta_arcsec"], float)
    ok = np.isfinite(d)
    chi = d / sigma
    any_ok = bool(ok.any())
    summary = {
        "n_images": len(ip),
        "n_solved": int(ok.sum()),
        "unsolved_images": [str(i) for i in ip["image_id"][~ok]],
        "shared_match_images": [str(i) for i in ip["image_id"][ip["shared_match"]]],
        "rms_dtheta_arcsec": float(np.sqrt(np.mean(d[ok] ** 2))) if any_ok else float("nan"),
        "max_dtheta_arcsec": float(d[ok].max()) if any_ok else float("nan"),
        "chi2_pos": float(np.sum(chi[ok] ** 2)),
        "chi2_pos_lenstool": chi2_ref,
        "images_over_3sigma": [str(i) for i in ip["image_id"][ok & (chi > 3)]],
        "grid_half_width_arcsec": float(grid.x[-1]),
        "grid_step_arcsec": float(grid.x[1] - grid.x[0]),
    }
    return ip, summary


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


def backtrace_check(model, images, z_m_limit, sigpos, chi2_ref: float | None):
    """Back-trace table plus a summary (χ², rms, per-system rms, > 3σ images).

    ``sigpos`` is one error for all images or one per image (arcsec).
    """
    bt = lensmodel.backtrace_images(model, images, z_m_limit)
    sig = np.broadcast_to(np.asarray(sigpos, float), (len(bt),))
    ok = np.isfinite(bt["dtheta_arcsec"])
    dtheta = np.asarray(bt["dtheta_arcsec"][ok])
    systems = {}
    for sys_id in np.unique(bt["system"][ok]):
        m = ok & (bt["system"] == sys_id)
        systems[str(sys_id)] = float(np.sqrt(np.mean(np.asarray(bt["dtheta_arcsec"][m]) ** 2)))
    outliers = bt[ok & (bt["dtheta_arcsec"] > 3 * sig)]
    summary = {
        "n_images": len(bt),
        "n_traced": int(ok.sum()),
        "rms_dtheta_arcsec": float(np.sqrt(np.mean(dtheta**2))),
        "median_dtheta_arcsec": float(np.median(dtheta)),
        "sigpos_arcsec": float(np.median(sig)),
        "chi2_pos": float(np.sum((dtheta / sig[ok]) ** 2)),
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


def map_check(model, map_path: Path, quantity: str, z_s: float | None, step: int = 37) -> dict:
    """Model against a published κ (D_LS/D_S = 1) or |μ| map (at ``z_s``) on a pixel sub-grid.

    Only map pixels with 0.05 < κ < 2, or |μ| < 10, are compared (critical curves excluded)."""
    with fits.open(map_path, memmap=True) as hdul:
        data = hdul[0].data
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FITSFixedWarning)
            wcs = WCS(hdul[0].header)
        jj, ii = np.mgrid[step // 2 : data.shape[0] : step, step // 2 : data.shape[1] : step]
        pub = np.asarray(data[jj, ii], float)
    ra, dec = wcs.pixel_to_world_values(ii, jj)
    if quantity == "kappa":
        val = model.kappa_xy(*model.to_frame(ra, dec))
        sel = (pub > 0.05) & (pub < 2)
    else:
        val = np.abs(np.asarray(model.evaluate(ra.ravel(), dec.ravel(), z_s)["magnification"]))
        val, pub = val.reshape(ra.shape), np.abs(pub)
        sel = (pub > 0) & (pub < 10)
    ok = sel & np.isfinite(val) & np.isfinite(pub)
    if not ok.any():
        raise SystemExit(f"error: {map_path}: no comparable pixel")
    rel = np.abs(val[ok] - pub[ok]) / pub[ok]
    return {
        "map": str(map_path),
        "quantity": quantity,
        "z_s": z_s,
        "n_points": int(ok.sum()),
        "median_rel_diff": float(np.median(rel)),
        "p95_rel_diff": float(np.percentile(rel, 95)),
        "median_ratio": float(np.median(val[ok] / pub[ok])),
    }


def cmd_validate_maps(args) -> dict:
    """Map model: κ and magnification against the published maps (no image list)."""
    model, _, _ = load_model(args.model)
    files = model_files(args.model)
    step = args.map_step
    summary = {
        "model": args.model,
        "model_sha256": model.sha256,
        "kappa_map": map_check(model, files["kappa_map"], "kappa", None, step),
    }
    for z_s, key in MODELS[args.model].get("mag_maps", {}).items():
        summary[f"magnification_z{z_s:g}"] = map_check(model, files[key], "mu", z_s, step)
    out = args.out / args.model
    out.mkdir(parents=True, exist_ok=True)
    if "arcs.dat" in MODELS[args.model]["files"]:  # validate measures it even when not trusted
        images, zml = image_list(args.model, files, None)
        bt = lensmodel.backtrace_images(model, images, zml)
        traced = images[np.isfinite(bt["beta_x"])]
        if not len(traced):
            raise SystemExit(f"error: {args.model}: no catalogued image has a usable redshift")
        half = grid_half_width(model, traced)
        grid = lensmodel.DeflectionGrid.cached(
            model, grid_cache_path(model, half, args.grid_step), half, args.grid_step
        )
        # no published position error for map models: unit sigma, so chi2_pos = sum dtheta^2;
        # the image-plane rms is what compares with the model's quoted rms
        ip, isum = imageplane_check(model, grid, bt, np.ones(len(bt)), None)
        for k in ("chi2_pos", "chi2_pos_lenstool", "images_over_3sigma"):
            isum[k] = None  # no published position error: only the rms is meaningful
        quoted = MODELS[args.model].get("quoted_rms_arcsec")
        if quoted:
            ok = isum["rms_dtheta_arcsec"] <= IMAGE_LIST_RMS_FACTOR * quoted
            isum["image_list_gate"] = {
                "quoted_rms_arcsec": quoted,
                "factor": IMAGE_LIST_RMS_FACTOR,
                "passes": bool(ok),
                "agrees_with_models_setting": bool(ok) == bool(has_image_list(args.model)),
            }
        summary["image_plane"] = isum
        _write(ip, out / "imageplane.ecsv")
    (out / "validate.json").write_text(json.dumps(summary, indent=1))
    return summary


def cmd_validate(args) -> dict:
    if is_map_model(args.model):
        return cmd_validate_maps(args)
    files = model_files(args.model)
    par = lensmodel.parse_lenstool_par(files["best.par"])
    model = lensmodel.LensModel.from_par(par)
    images = lensmodel.load_lenstool_images(files["arcs.dat"])
    sigma = model_sigpos(args.model, files, images)
    chi2_ref = chi2pos_from_par(files["best.par"])
    image_plane_opt = (
        "image plane optimization" in files["best.par"].read_text(encoding="latin-1").lower()
    )
    bt, bsum = backtrace_check(model, images, par["z_m_limit"], sigma, chi2_ref)
    half = grid_half_width(model, images)
    grid = lensmodel.DeflectionGrid.cached(
        model, grid_cache_path(model, half, args.grid_step), half, args.grid_step
    )
    # a source-plane fit's Chi2pos is not an image-plane chi2: no reference for that comparison
    ip, isum = imageplane_check(model, grid, bt, sigma, chi2_ref if image_plane_opt else None)
    summary = {
        "model": args.model,
        "model_sha256": model.sha256,
        "n_potentials": len(model.components),
        "image_plane_optimised": image_plane_opt,
        "kappa_map": (
            kappa_map_check(model, files["kappa_map"], step=args.step)
            if "kappa_map" in files
            else None
        ),
        "backtrace": bsum,
        "image_plane": isum,
    }
    out = args.out / args.model
    _write(bt, out / "backtrace.ecsv")
    _write(ip, out / "imageplane.ecsv")
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
    model, files, par = load_model(args.model)
    apply_frame_offset(args.model, model)  # every model, maps included (D-040)
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

    images = None
    if par is not None:
        images = shift_images(args.model, lensmodel.load_lenstool_images(files["arcs.dat"]))
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
        "multiple_images": (
            convention_check(model, par, images, shapes, args.image_match_arcsec)
            if par is not None
            else None  # map model: no multiple-image list
        ),
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


def aperture_snr(
    img, err, xx, yy, cx: float, cy: float, r_ap: float = FORCED_R_AP
) -> tuple[float, float]:
    """Background-subtracted flux and error in an ``r_ap`` aperture at offset (cx, cy).

    ``xx``/``yy`` are each pixel's offset from the target in arcsec (West, North). NaN when fewer
    than ``MIN_VALID`` of the aperture or annulus pixels are finite (gaps, edges, masks)."""
    r = np.hypot(xx - cx, yy - cy)
    ap = r <= r_ap
    ann = (r > FORCED_ANNULUS[0]) & (r < FORCED_ANNULUS[1])
    good = np.isfinite(img) & np.isfinite(err)
    if not ap.any() or good[ap].mean() < MIN_VALID or good[ann].mean() < MIN_VALID:
        return float("nan"), float("nan")
    bkg = float(np.nanmedian(img[ann & good]))
    flux = float(np.nansum(img[ap] - bkg))
    ferr = float(np.sqrt(np.nansum(err[ap] ** 2))) * ERR_SCALE
    return flux, ferr


def peak_position(img, xx, yy, radius: float = 0.4) -> tuple[float, float]:
    """Offset (arcsec) of the brightest finite pixel within ``radius`` of the stamp centre."""
    inner = np.where((np.hypot(xx, yy) <= radius) & np.isfinite(img), img, -np.inf)
    if not np.isfinite(inner).any():
        return float("nan"), float("nan")
    k = int(np.argmax(inner))
    return float(xx.flat[k]), float(yy.flat[k])


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


MAX_SEARCH_ARCSEC = 3.0  # cap of the residual-scaled search radius (ASSUMPTION)
MAX_REF_SPREAD = 3.0  # sibling references disagreeing by more than this are inconsistent
MAX_PEAK_OVER_POSITION = 2.0  # the 0.4" peak may not be much brighter than the catalogued spot
MIN_REF_CONCENTRATION = 0.6  # f(0.2")/f(0.4") of a reference image (as exotic_screens, D-031)


def system_search_radius(table: Table, sys_id: str, base: float, has_unpredicted: bool) -> float:
    """Search radius for a system's predicted images: ``base``, or 1.5x the largest offset of
    its catalogued images from their predictions when that is larger (model position error),
    and 1.5x the match radius when one of its catalogued images has no prediction at all. The
    widening (not ``base`` itself) is capped at ``MAX_SEARCH_ARCSEC`` (D-034)."""
    obs = table[(np.asarray(table["system"]).astype(str) == sys_id)]
    obs = obs[obs["image_class"] == "observed"]
    widen = 0.0
    if len(obs) and "sep_image_arcsec" in obs.colnames:
        widen = 1.5 * float(np.nanmax(obs["sep_image_arcsec"]))
    if has_unpredicted:  # the model is off here
        match = float(table.meta.get("assumptions", {}).get("match_arcsec", 1.5))
        widen = max(widen, 1.5 * match)
    return max(base, min(widen, MAX_SEARCH_ARCSEC))


def forced_check(
    table: Table,
    backtrace: Table,
    stamp,
    search_arcsec: float = 1.0,
    max_ref_mu: float = 50.0,
    unpredicted: list[str] | None = None,
) -> None:
    """Add forced-photometry columns to the predicted-image ``table`` in place.

    ``stamp(ra, dec)`` returns ``(sci, err, xx, yy)`` around a position (it must cover
    ``MAX_SEARCH_ARCSEC`` plus the background annulus).

    Reference flux (rules from El Gordo, D-034; thresholds are ASSUMPTIONs): a catalogued image
    with |μ| <= ``max_ref_mu``, S/N > 5 at its catalogued position, a 0.4"-recentred peak at
    most ``MAX_PEAK_OVER_POSITION`` x brighter (else the peak is a neighbour), and compact
    (f(0.2")/f(0.4") >= ``MIN_REF_CONCENTRATION``). Of those, the
    least magnified is the reference. When two or more qualify and their f/|μ| differ by more
    than ``MAX_REF_SPREAD``, the system is ``inconsistent_reference`` (a flux-ratio question for
    ``exotic_screens.py fluxratio``, not a counter-image test). The search radius grows with
    the system's own model residuals (``system_search_radius``)."""
    n = len(table)
    keys = ("pred_snr", "best_snr", "flux_ratio", "best_dx", "best_dy", "search_arcsec")
    cols = {k: np.full(n, np.nan) for k in keys}
    fclass = np.full(n, "", dtype="U22")
    ref_cache: dict[str, tuple[float, float, str]] = {}
    systems = np.asarray(backtrace["system"]).astype(str)
    for i, row in enumerate(table):
        if row["image_class"] in ("observed", "demagnified", "outside"):
            continue
        sys_id = str(row["system"])
        if sys_id not in ref_cache:
            usable = []
            for b in backtrace[systems == sys_id]:
                mu = abs(float(b["magnification"]))
                if not np.isfinite(mu) or mu > max_ref_mu:
                    continue
                img, err, xx, yy = stamp(float(b["ra"]), float(b["dec"]))
                f, e = peak_flux(img, err, xx, yy)
                f0, e0 = aperture_snr(img, err, xx, yy, 0.0, 0.0)
                if not (np.isfinite(f) and np.isfinite(f0) and e > 0 and e0 > 0):
                    continue
                if not (f / e > 5 and f0 / e0 > 5 and f <= MAX_PEAK_OVER_POSITION * max(f0, 0)):
                    continue
                # compact only: a resolved arc's (or a neighbour's wing) aperture flux does not
                # scale with |mu| (surface brightness is conserved; D-031)
                px, py = peak_position(img, xx, yy)
                f_big, _ = aperture_snr(img, err, xx, yy, px, py, r_ap=0.4)
                if np.isfinite(f_big) and f_big > 0 and f / f_big >= MIN_REF_CONCENTRATION:
                    usable.append((f, mu))
            if not usable:
                ref_cache[sys_id] = (np.nan, np.nan, "no_reference")
            else:
                lum = [f / mu for f, mu in usable]
                f_ref, mu_ref = min(usable, key=lambda t: t[1])
                state = "inconsistent_reference" if max(lum) > MAX_REF_SPREAD * min(lum) else ""
                ref_cache[sys_id] = (f_ref, mu_ref, state)
        f_ref, mu_ref, state = ref_cache[sys_id]
        if state:
            fclass[i] = state
            continue
        sci, err, xx, yy = stamp(float(row["ra"]), float(row["dec"]))
        _, e0 = aperture_snr(sci, err, xx, yy, 0.0, 0.0)
        if not (np.isfinite(e0) and e0 > 0):
            fclass[i] = "off_image"  # off the footprint, or in a gap or masked region
            continue
        has_unpred = any(lensmodel.image_family(u) == sys_id for u in (unpredicted or []))
        radius = system_search_radius(table, sys_id, search_arcsec, has_unpred)
        pred = f_ref * abs(float(row["magnification"])) / mu_ref
        cols["pred_snr"][i] = pred / e0
        snr, flux, dx, dy = best_within(sci, err, xx, yy, radius)
        cols["best_snr"][i], cols["best_dx"][i], cols["best_dy"][i] = snr, dx, dy
        cols["search_arcsec"][i] = radius
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
        "max_search_arcsec": MAX_SEARCH_ARCSEC,
        "max_ref_mu": max_ref_mu,
        "max_ref_spread": MAX_REF_SPREAD,
        "max_peak_over_position": MAX_PEAK_OVER_POSITION,
        "min_ref_concentration": MIN_REF_CONCENTRATION,
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
    if not has_image_list(args.model):
        raise SystemExit(
            f"error: {args.model}: no usable multiple-image list (none published, or excluded by "
            "the image-plane rms gate of D-035)"
        )
    model, files, par = load_model(args.model)
    images, zml = image_list(args.model, files, par)
    dra, ddec = apply_frame_offset(args.model, model, images)  # into the JWST frame
    bt = lensmodel.backtrace_images(model, images, zml)
    traced = images[np.isfinite(bt["beta_x"])]
    if not len(traced):
        raise SystemExit(f"error: {args.model}: no catalogued image has a usable redshift")
    half = args.half_width or grid_half_width(model, traced)  # as validate: images + 20"
    grid = lensmodel.DeflectionGrid.cached(
        model, grid_cache_path(model, half, args.step), half, args.step
    )
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
            args.forced_image,
            max(args.forced_search_arcsec, MAX_SEARCH_ARCSEC) + FORCED_ANNULUS[1] + 0.1,
        )
        try:
            forced_check(table, bt, stamp, args.forced_search_arcsec, unpredicted=unpredicted)
        finally:
            close()
        table.meta["forced"]["image"] = args.forced_image
        table.meta["forced"]["frame_offset_arcsec"] = [dra, ddec]
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
                        "inconsistent_reference",
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
                        "search_arcsec": round(float(r["search_arcsec"]), 2),
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


def load_dja_photometry(
    path: Path, bands: list[str], aperture: int = 1, zout: Path | None = None
) -> Table:
    """DJA ``fix_phot`` for ``fluxratios``: ``photometry.load_dja_catalog`` plus ``magerr_auto``,
    ``npix`` and, with ``zout`` (eazy, row-aligned by ``id``), ``z_phot, z025, z975``.

    ``mag_auto`` (detection-image Kron) is the only total: ``<band>_tot_corr`` is 1 and
    ``tot_corr`` a point-source correction capped near 1.2 in DJA v7, so aperture magnitudes
    are not totals for extended images.
    """
    out = load_dja_catalog(path, bands, aperture)
    raw = Table.read(path)
    missing = [c for c in ("mag_auto", "magerr_auto", "npix") if c not in raw.colnames]
    if missing:
        raise SystemExit(f"error: {path} has no {missing}")
    out["magerr_auto"] = np.asarray(np.ma.filled(raw["magerr_auto"], np.nan), float)
    # Unknown segment size counts as too large (unusable), not as a clean compact source.
    out["npix"] = np.asarray(np.ma.filled(raw["npix"], np.iinfo(np.int32).max), int)
    out.rename_column("id", "dja_id")
    if zout is not None:
        z = Table.read(zout)
        if len(z) != len(out) or not np.array_equal(np.asarray(z["id"]), out["dja_id"]):
            raise SystemExit(f"error: {zout} rows do not match {path} ids")
        for c in ("z_phot", "z025", "z975"):
            v = np.asarray(np.ma.filled(z[c], np.nan), float)
            out[c] = np.where(v >= 0, v, np.nan)  # eazy writes -1 for no fit
        out.meta["photoz_source"] = str(zout)
    out.meta.update(aperture=aperture)
    return out


def _loo_residual(
    values: np.ndarray, errors: np.ndarray, use: np.ndarray, systems: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Leave-one-out residuals: ``values`` minus the median of the *other* usable images of its
    system (NaN if < 2 usable), so a discrepant image does not pull its own reference; in a pair
    each image carries the full pair difference.

    Returns ``(resid, resid_err, worst)``. ``resid_err`` adds the reference's error (the others'
    median error / sqrt(n), an approximation) in quadrature. ``worst`` marks each system's
    largest |resid| (both images of a pair), so one bad image in a system of three or more
    does not also flag the good ones, whose references it shifts."""
    resid = np.full(len(values), np.nan)
    rerr = np.full(len(values), np.nan)
    worst = np.zeros(len(values), bool)
    for sys_id in dict.fromkeys(systems):
        ii = np.flatnonzero((systems == sys_id) & use)
        if ii.size < 2:
            continue
        for i in ii:
            others = ii[ii != i]
            resid[i] = values[i] - np.median(values[others])
            rerr[i] = np.hypot(errors[i], np.median(errors[others]) / np.sqrt(others.size))
        a = np.abs(resid[ii])
        worst[ii[np.isclose(a, a.max())]] = True
    return resid, rerr, worst


def _r(v, nd: int = 2):
    """Round for JSON; NaN becomes None (JSON has no NaN)."""
    v = float(v)
    return round(v, nd) if np.isfinite(v) else None


def flux_ratio_table(
    model,
    images: Table,
    z_m_limit: dict[str, float],
    phot: Table,
    colour: tuple[str, str] = ("f150w", "f444w"),
    *,
    match_arcsec: float = 0.3,
    offset_arcsec: tuple[float, float] = (0.0, 0.0),
    min_snr: float = 10.0,
    max_abs_mu: float = 20.0,
    outlier_mag: float = 0.75,
    colour_tol: float = 0.3,
    max_npix: int = 20000,
    z_margin: float = 0.15,
) -> Table:
    """Per-image source-magnitude and colour residuals (see ``fluxratios``)."""
    z = lensmodel.image_redshifts(images, z_m_limit)
    ok_z = np.isfinite(z)
    mu = np.full(len(images), np.nan)
    if ok_z.any():
        pred = model.evaluate(images["ra"][ok_z], images["dec"][ok_z], z[ok_z])
        mu[ok_z] = np.asarray(pred["magnification"], float)
    ci = SkyCoord(images["ra"], images["dec"], unit="deg")
    if any(offset_arcsec):  # catalogue frame minus image-list frame, (dRA cos dec, dDec)
        ci = ci.spherical_offsets_by(offset_arcsec[0] * u.arcsec, offset_arcsec[1] * u.arcsec)
    cp = SkyCoord(phot["ra"], phot["dec"], unit="deg")
    idx, sep, _ = ci.match_to_catalog_sky(cp)
    matched = sep.arcsec <= match_arcsec
    ids = np.where(matched, np.asarray(phot["dja_id"])[idx], -1)
    uniq, counts = np.unique(ids[matched], return_counts=True)
    blended = matched & np.isin(ids, uniq[counts > 1])
    npix = np.where(matched, np.asarray(phot["npix"])[idx], -1)
    large = npix > max_npix  # segment swallowing host or ICL light (SMACS 1.1)

    def col(name: str) -> np.ndarray:
        return np.where(matched, np.asarray(phot[name], float)[idx], np.nan)

    # A counterpart whose 95 % photo-z interval excludes the system redshift is another object
    # (El Gordo 9a: z_phot 0.89 for a z = 4.32 system). Missing photo-z excludes nothing.
    pz_excl = np.zeros(len(images), bool)
    if "z025" in phot.colnames:
        dz = z_margin * (1.0 + z)
        with np.errstate(invalid="ignore"):
            pz_excl = (z < col("z025") - dz) | (z > col("z975") + dz)
    bad = blended | large | pz_excl
    mag, magerr = col("mag_auto"), col("magerr_auto")
    m1, m2 = col(f"{colour[0]}_mag"), col(f"{colour[1]}_mag")
    e1, e2 = col(f"{colour[0]}_mag_err"), col(f"{colour[1]}_mag_err")
    with np.errstate(divide="ignore", invalid="ignore"):
        src = mag + 2.5 * np.log10(np.abs(mu))
        snr_auto = np.where(magerr > 0, 1.0857 / magerr, np.nan)  # SEP: Pogson magerr_auto
    snr1, snr2 = (np.where(e > 0, snr_from_mag_err(e), np.nan) for e in (e1, e2))
    c_obs = m1 - m2
    c_err = np.hypot(e1, e2)
    systems = np.asarray(images["system"]).astype(str)
    base = matched & ~bad & np.isfinite(mu) & (np.abs(mu) <= max_abs_mu)
    use_m = base & np.isfinite(src) & (snr_auto >= min_snr)
    use_c = matched & ~bad & np.isfinite(c_obs) & (snr1 >= min_snr) & (snr2 >= min_snr)
    resid, resid_err, worst = _loo_residual(src, magerr, use_m, systems)
    c_resid, c_resid_err, c_worst = _loo_residual(c_obs, c_err, use_c, systems)
    cname = f"{colour[0]}_{colour[1]}"
    out = Table(
        {
            "image_id": images["image_id"],
            "system": images["system"],
            "ra": images["ra"],
            "dec": images["dec"],
            "z_used": z,
            "magnification": mu,
            "parity": np.sign(mu),
            "dja_id": ids,
            "sep_arcsec": sep.arcsec,
            "npix": npix,
            "blended": blended,
            "large_segment": large,
            "photoz_excluded": pz_excl,
            "mag_auto": mag,
            "source_mag": src,
            "resid_mag": resid,
            "resid_err": resid_err,
            f"colour_{cname}": c_obs,
            "colour_resid": c_resid,
            "colour_err": c_resid_err,
        }
    )
    cls = np.full(len(images), "untested", dtype="U16")
    cls[np.isfinite(resid)] = "consistent"
    cls[worst & (np.abs(resid) > outlier_mag)] = "flux_outlier"
    out["flux_class"] = cls
    ccls = np.full(len(images), "untested", dtype="U16")
    ccls[np.isfinite(c_resid)] = "consistent"
    ccls[c_worst & (np.abs(c_resid) > colour_tol)] = "colour_outlier"
    out["colour_class"] = ccls
    out.meta.update(model._meta())
    out.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source=f"{images.meta.get('source', 'arcs.dat')} and {phot.meta.get('source')}"
        f" against {model.source}",
        thresholds={
            "match_arcsec": match_arcsec,
            "offset_arcsec": list(offset_arcsec),
            "min_snr": min_snr,
            "max_abs_mu": max_abs_mu,
            "outlier_mag": outlier_mag,
            "colour_tol": colour_tol,
            "colour": list(colour),
            "max_npix": max_npix,
            "z_margin": z_margin,
            "photoz": "z025-z975" if "z025" in phot.colnames else None,
        },
    )
    return out


def _resid_stats(r: np.ndarray) -> dict:
    r = r[np.isfinite(r)]
    if not r.size:
        return {"n_images": 0}
    return {
        "n_images": int(r.size),
        "median_abs": round(float(np.median(np.abs(r))), 3),
        "rms": round(float(np.sqrt(np.mean(r**2))), 3),
    }


def cmd_fluxratios(args) -> dict:
    files = model_files(args.model)
    par = lensmodel.parse_lenstool_par(files["best.par"])
    model = lensmodel.LensModel.from_par(par)
    images = lensmodel.load_lenstool_images(files["arcs.dat"])
    colour = tuple(b.strip().lower() for b in args.colour.split(","))
    if len(colour) != 2:
        raise SystemExit("--colour takes two bands, e.g. f150w,f444w")
    phot = load_dja_photometry(args.photometry, list(colour), args.aperture, args.photoz)
    offset = tuple(args.offset_arcsec)
    table = flux_ratio_table(
        model,
        images,
        par["z_m_limit"],
        phot,
        colour,
        match_arcsec=args.match_arcsec,
        offset_arcsec=offset,
        min_snr=args.min_snr,
        max_abs_mu=args.max_abs_mu,
        outlier_mag=args.outlier_mag,
        colour_tol=args.colour_tol,
        max_npix=args.max_npix,
        z_margin=args.z_margin,
    )
    out = args.out / args.model
    _write(table, out / "flux_ratios.ecsv")
    summary = {
        "model": args.model,
        "photometry": str(args.photometry),
        "thresholds": table.meta["thresholds"],
        "n_images": len(table),
        "n_matched": int(np.sum(table["dja_id"] >= 0)),
        "n_blended": int(np.sum(table["blended"])),
        "n_large_segment": int(np.sum(table["large_segment"])),
        "n_photoz_excluded": int(np.sum(table["photoz_excluded"])),
        "flux_classes": {
            c: int(np.sum(table["flux_class"] == c))
            for c in ("consistent", "flux_outlier", "untested")
        },
        "colour_classes": {
            c: int(np.sum(table["colour_class"] == c))
            for c in ("consistent", "colour_outlier", "untested")
        },
        "flux_resid_mag": _resid_stats(np.asarray(table["resid_mag"], float)),
        "colour_resid_mag": _resid_stats(np.asarray(table["colour_resid"], float)),
        "flagged": [
            {
                "image_id": str(r["image_id"]),
                "flux_class": str(r["flux_class"]),
                "colour_class": str(r["colour_class"]),
                "mu": _r(r["magnification"]),
                "resid_mag": _r(r["resid_mag"]),
                "colour_resid": _r(r["colour_resid"]),
                "photoz_excluded": bool(r["photoz_excluded"]),
                "npix": int(r["npix"]),
            }
            for r in table
            if r["flux_class"] == "flux_outlier" or r["colour_class"] == "colour_outlier"
        ],
    }
    (out / "fluxratios.json").write_text(json.dumps(summary, indent=1))
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--model", choices=sorted(MODELS), default="smacs0723-iclv2")
    ap.add_argument("--out", type=Path, default=paths.outputs_dir() / "lens_consistency")
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate", help="compare with the published kappa map and multiple images")
    v.add_argument("--step", type=int, default=5, help="kappa-map sub-grid step in pixels")
    v.add_argument("--map-step", type=int, default=37, help="map-model check sub-grid step, px")
    v.add_argument("--grid-step", type=float, default=0.25, help="solver grid step, arcsec")
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
    i.add_argument(
        "--half-width", type=float, help='solver grid half-width, arcsec (default: images + 20")'
    )
    i.add_argument("--step", type=float, default=0.1, help="solver grid step, arcsec")
    i.add_argument(
        "--forced-image", help="_i2d URI or path for forced photometry (S3 by byte range)"
    )
    i.add_argument("--forced-search-arcsec", type=float, default=1.0)
    f = sub.add_parser("fluxratios", help="image flux ratios against model magnification ratios")
    f.add_argument("--photometry", type=Path, required=True, help="DJA fix_phot FITS table")
    f.add_argument("--colour", default="f150w,f444w", help="two bands for the colour test")
    f.add_argument("--aperture", type=int, default=1, help="DJA aperture index (1 = 0.5 arcsec)")
    f.add_argument("--match-arcsec", type=float, default=0.3)
    f.add_argument(
        "--offset-arcsec",
        nargs=2,
        type=float,
        default=(0.0, 0.0),
        metavar=("DRA", "DDEC"),
        help="photometry frame minus image list, arcsec (El Gordo: 0.221 -0.018)",
    )
    f.add_argument("--min-snr", type=float, default=10.0)
    f.add_argument("--max-abs-mu", type=float, default=20.0)
    f.add_argument("--outlier-mag", type=float, default=0.75)
    f.add_argument("--colour-tol", type=float, default=0.3)
    f.add_argument("--max-npix", type=int, default=20000, help="larger DJA segments are blends")
    f.add_argument("--photoz", type=Path, help="DJA eazy zout matching --photometry row by row")
    f.add_argument("--z-margin", type=float, default=0.15, help="photo-z margin, x (1 + z)")
    args = ap.parse_args(argv)
    summary = {
        "validate": cmd_validate,
        "arcs": cmd_arcs,
        "images": cmd_images,
        "fluxratios": cmd_fluxratios,
    }[args.cmd](args)
    json.dump(summary, sys.stdout, indent=1)
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
