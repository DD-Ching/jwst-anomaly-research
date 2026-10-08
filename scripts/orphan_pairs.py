"""Blind search for "orphan image pairs" behind a lensing cluster or in a deep field
(docs/orphan_pairs.md).

An orphan pair is two background sources 0.3-3" apart whose SEDs and photo-z are statistically
identical, with no catalogued galaxy between them and no published multiple image nearby. It is a
*possible* signature of an unseen compact deflector (a "dark lens") - a hypothesis, never a
conclusion. The search is lens-model independent; a lens model is used afterwards only to test
the ordinary explanation (a cluster-scale image pair near a critical curve). Deep fields have no
cluster lens model: there the ordinary-lensing test is "no lens model needed" (catalogue mu, or 1).

Catalogues are read through an adapter (``as_standard``) into one column layout, so CANUCS DR1
(``FLUX_COLOR03_TOTAL_*``, EAzY columns in the catalogue) and DJA grizli (``<band>_flux_aper_N``
plus a separate eazy-py ``zout``) run through the same code (D-051).

Steps (every threshold is an ASSUMPTION; defaults below):

1. Sources: ``USE_PHOT_APER03``, not ``FLAG_BCG``, S/N >= ``MIN_SNR`` in the summed F277W+F356W+
   F444W colour-aperture flux (the bands of those three the catalogue has), S/N >= ``HI_SNR`` in
   at least ``MIN_HI_BANDS`` bands (without it the chi^2 test has no power: 25 % of random faint
   pairs "match"), and behind the cluster (``Z025 > z_cluster + Z_MARGIN``, or
   ``Z_SPEC > z_cluster + Z_MARGIN`` when a spectroscopic redshift exists). Deep fields use the
   notional deflector redshift ``Z_LENS_REF`` in place of z_cluster.
2. Pairs at ``SEP_MIN``-``SEP_MAX`` arcsec. SED match: chi^2 of one SED against a free multiple of
   the other over every band valid in both (``FLUX_COLOR03_TOTAL_*``, error floor ``ERR_FLOOR``),
   at least ``MIN_BANDS`` bands, chi^2 survival probability >= ``P_MIN``, and overlapping
   16-84 % photo-z intervals.
3. Classes (first that applies): ``published`` (a member within ``IMAGE_RADIUS`` of a published
   multiple image), ``same_galaxy`` (separation < ``SAME_GALAXY_FACTOR`` x the sum of the Kron
   aperture radii, or both deblended with overlapping segment boxes), ``visible_lens`` (another
   catalogued source of any redshift within ``LENS_RADIUS`` of the midpoint, or between the
   members - within ``LENS_RADIUS`` of the joining segment or inside the circle with the pair as
   diameter - while clear of both members by ``LENS_RADIUS``), else ``orphan``.
4. Null: the SED-match fraction of (a) cross pairs between the catalogue and copies shifted by
   10-60" and (b) real pairs at ``FAR_MIN``-``FAR_MAX`` arcsec, times the number of close pairs,
   is the expected number of chance matches; per class, the fraction times the close pairs that
   the geometric class rules put in that class. (c) The match rate among photo-z-overlapping far
   pairs times the photo-z-overlapping close pairs tests whether an excess is just redshift
   clustering (physical neighbours share SEDs). (d) is (c) at 3-6" and (e) is (c) conditioned on
   the pair's S/N, size and colour (``pair_cells``; D-051).
5. ``--cutouts``: F150W/F277W/F444W cutouts of the top orphans (S3 byte ranges of MAST level-3
   ``_i2d`` files) on one contact sheet, plus a lens-model check (CATS map magnification and
   parity at both members; deep fields: the catalogue magnification).

Outputs are ``derived`` (from observed catalogues) with ``model_prediction`` columns from the lens
model; results go to ``outputs/orphan_pairs/<field>/``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import sys
import tarfile
import tempfile
import urllib.request
from pathlib import Path

import numpy as np
from astropy.table import Table
from scipy.spatial import cKDTree
from scipy.stats import chi2 as chi2_dist

sys.path.insert(0, str(Path(__file__).resolve().parent))

from jwst_anomaly import paths, schema  # noqa: E402
from jwst_anomaly.photometry import BAD_FLAGS, fetch_catalog  # noqa: E402

SEP_MIN, SEP_MAX = 0.3, 3.0  # arcsec, the orphan-pair annulus (task definition)
FAR_MIN, FAR_MAX = 10.0, 30.0  # arcsec, null annulus (b) (ASSUMPTION)
NEAR_MIN, NEAR_MAX = 3.0, 6.0  # arcsec, null annulus (d): small-scale clustering (ASSUMPTION)
SHIFT_RANGE = (10.0, 60.0)  # arcsec, null (a) shift lengths (ASSUMPTION)
MIN_SNR = 10.0  # summed F277W+F356W+F444W colour-aperture S/N (ASSUMPTION)
SNR_BANDS = ("F277W", "F356W", "F444W")
Z_MARGIN = 0.1  # behind the cluster: Z025 (or Z_SPEC) > z_cluster + this (ASSUMPTION)
ERR_FLOOR = 0.03  # fractional flux-error floor added in quadrature (ASSUMPTION)
MIN_BANDS = 8  # bands valid in both members (ASSUMPTION)
MIN_HI_BANDS = 8  # bands with S/N >= HI_SNR per source: an informative SED (ASSUMPTION)
HI_SNR = 10.0
P_MIN = 0.01  # chi^2 survival probability for "identical SED" (ASSUMPTION)
AGREE_SIGMA = 2.0  # a band "agrees" when its scaled residual is below this (ASSUMPTION)
IMAGE_RADIUS = 1.0  # arcsec to a published multiple image (task definition)
LENS_RADIUS = 0.3  # arcsec from the midpoint to a catalogued (visible) lens (task definition)
SAME_GALAXY_FACTOR = 0.5  # x sum of Kron aperture radii (task definition)
KRON_SCALE = 2.5  # Kron aperture = 2.5 x the first-moment Kron radius (SExtractor convention)
PIXSCALE = 0.04  # arcsec/pixel of the CANUCS detection image (fitted from X, Y vs RA, Dec)
MAX_ORDINARY_MU = 10.0  # |mu| above this: a cluster-scale image pair is ordinary (task definition)
CUTOUT_BANDS = ("f150w", "f277w", "f444w")
CUTOUT_ARCSEC = 4.0
# Deep fields (D-051): a notional deflector redshift. Sources must lie behind it by Z_MARGIN, so
# the redshift floor is z_low > 0.5, the same as behind a z = 0.4 cluster (ASSUMPTION).
Z_LENS_REF = 0.4
DJA_APERTURE = 0  # DJA aperture index: 0.36" diameter, the closest to CANUCS's 0.3" (ASSUMPTION)
DJA_DEBLENDED = 0x0001  # SEP object flag OBJ_MERGED: the object was deblended (sep.h)
# DJA aperture radius for the same_galaxy rule = this x flux_radius (half-light radius, px): the
# median ratio of the D-048 Kron aperture radius to FLUX_RADIUS among selected CANUCS sources
# (3.2-3.5 in MACS0416, its NCF and the MACS1149 NCF; ASSUMPTION).
DJA_KRON_PER_R50 = 3.3
# MIRI bands are left out of the DJA SED: 0.36" apertures are smaller than the MIRI PSF core.
MIRI_BANDS = {"F560W", "F770W", "F1000W", "F1130W", "F1280W", "F1500W", "F1800W", "F2100W"}


_HLSP = "https://archive.stsci.edu/hlsps/canucs/dr1/{f}/{d}/hlsp_canucs_jwst-hst_multi_{f}-{n}"
# Inputs are downloaded once and verified by sha256 (SOURCES.md, D-048). ``image_lists`` are
# extra (url, sha256) lists already in the JWST frame; ``image_models`` are lens_consistency
# MODELS whose pinned image list (and frame offset) is used.
FIELDS = {
    "macs0416": {
        "catalogue": (
            _HLSP.format(f="macs0416", d="clu", n="clu_multi_v1_photometry-cat.fits.gz"),
            "339107a5c5041621d7bed4ecc8b4a51b5148a6913d4e513c3a5e783363f11bff",
        ),
        "z_cluster": 0.396,
        "cats": "macs0416-cats",
        "image_lists": [
            (  # CANUCS frame offset under 0.1", not applied (D-044)
                _HLSP.format(f="macs0416", d="model", n="allmultim-cat_multi_v1_model.txt"),
                "c8978003d8dd617cb980ed7ba5acde1485cd742db43846c25ed251f110dc2417",
            )
        ],
        "i2d": "jw01208-o004_t002",
    },
    "macs1149": {
        "catalogue": (
            _HLSP.format(f="macs1149", d="clu", n="clu_multi_v1_photometry-cat.fits.gz"),
            "08ab67347f2c3dfe4f743cc1b2a9d4b77a1e40eb66ddeb611648bb296adc0739",
        ),
        "z_cluster": 0.543,
        "cats": "macs1149-cats",
        "image_lists": [],  # CANUCS DR1 publishes no MACS1149 image list (readme: with v2)
        "i2d": "jw01208-o008_t004",
    },
    "abell370": {
        "catalogue": (
            _HLSP.format(f="a370", d="clu", n="clu_multi_v1_photometry-cat.fits.gz"),
            "f5622f2867aa3094df6861981ee67f7f1879b8381b348224748c7381436ae17b",
        ),
        "z_cluster": 0.375,
        "cats": "abell370-cats",
        "image_lists": [],
        # the pinned CANUCS Lenstool list, moved by its (-0.148, 0.002)" frame offset (D-044).
        # Its image-plane gate (image_list_ok False) concerns model constraints; here only the
        # positions are used, to mark pairs within 1" of a published image.
        "image_models": ["abell370-canucs"],
        "i2d": "jw01208-o002_t001",
    },
}

# Deep (non-cluster) fields, D-051. No published multiple images and no cluster lens model.
# ``mast``: the program and level-3 observation prefixes whose ``_i2d`` files serve the cutouts
# (the observation covering each pair is chosen from its MAST footprint).
_NCF = "https://archive.stsci.edu/hlsps/canucs/dr1/{c}/ncf/hlsp_canucs_jwst-hst_multi_{c}-ncf_multi_v1_photometry-cat.fits.gz"
_DJA = "https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/"
DEEP_FIELDS = {
    # CANUCS NIRCam flanking fields: same catalogue format as the cluster fields, ~3-6' from
    # the cluster cores; the catalogue MU (CANUCS lens model) is 1.00-1.4 there.
    "macs0416-ncf": {
        "catalogue": (
            _NCF.format(c="macs0416"),
            "72a2c609015f830af4ca05bcd49c7754723c0a55434176f780508bfe62a5afe3",
        ),
        "mast": ("1208", ["jw01208-o023_t002"]),
    },
    "macs1149-ncf": {
        "catalogue": (
            _NCF.format(c="macs1149"),
            "f83510e39332657bf149a08d8060b8dd63a9bb9b839ca5d157714c88ea273ceb",
        ),
        "mast": ("1208", ["jw01208-o029_t004"]),
    },
    "abell370-ncf": {
        "catalogue": (
            _NCF.format(c="a370"),
            "b35c6eed37ae964d7f996a563dca5794ccfa8ccb0e8b364800fc04b54aa3039f",
        ),
        "mast": ("1208", ["jw01208-o020_t001"]),
    },
    "macs0417-ncf": {
        "catalogue": (
            _NCF.format(c="macs0417"),
            "f82940179ebf63ddf2ec872e7a0297f2fce1785faa0de83ee43b401e96acc49d",
        ),
        "mast": ("1208", ["jw01208-o026_t003"]),
    },
    "macs1423-ncf": {
        "catalogue": (
            _NCF.format(c="macs1423"),
            "b4881e0ea7fc4dd11fbad76afa73a3ca56f0848d0f85384f59646070c9b7b63b",
        ),
        "mast": ("1208", ["jw01208-o032_t005"]),
    },
    # DJA v7.3 GOODS-North (JADES/FRESCO area). The catalogue is 223.8 MB and the photo-z
    # tarball 371.1 MB: both over 200 MB, stated reason in D-051. Only the tarball's zout member
    # is kept (streamed; the tarball itself is never written to disk). The standalone
    # ``gdn-grizli-v7.3-fix.eazypy.zout.fits`` on S3 belongs to an older catalogue (63,069 rows,
    # ids do not match): do not use it.
    "goodsn-dja": {
        "format": "dja",
        "catalogue": (
            _DJA + "gdn-grizli-v7.3-fix_phot.fits",
            "9b18b41731c3a86085cb9c4fdb7a4c9f15c5477431f4eea904d410a1f173b6c1",
        ),
        "max_bytes": 230_000_000,
        "photoz": {
            "url": _DJA + "gdn-grizli-v7.3-fix.photoz.tar.gz",
            "sha256": "1d89eef3f613eeb7d592ecf03d94869dfd76110c8208eb728291592ddca8db82",
            "member": "gdn-grizli-v7.3-fix.eazypy.zout.fits",
            "member_sha256": "363176053431708410d4957a77e56826f5ced3d6ed072722aeccc9a2edde48fc",
            "max_bytes": 380_000_000,
        },
        "mast": ("1181", ["jw01181-"]),
    },
}
for _spec in DEEP_FIELDS.values():
    _spec.setdefault("format", "canucs")
    _spec.update(deep=True, z_cluster=None, cats=None, image_lists=[])
FIELDS.update(DEEP_FIELDS)


def z_lens_of(spec: dict) -> float:
    """The deflector redshift a field's background cut is measured from (deep: Z_LENS_REF)."""
    return spec["z_cluster"] if spec.get("z_cluster") is not None else Z_LENS_REF


# ------------------------------------------------------------------------------- catalogue


# Every catalogue is read into this layout (``as_standard``); the search uses only these columns.
# ``flux``/``err`` are (N, B) arrays in nJy with invalid entries NaN, bands in ``meta["bands"]``.
STANDARD_COLUMNS = (
    "src_id",  # catalogue id
    "ra",
    "dec",
    "x_pix",  # detection-image pixel position and segment box (box overlap = same segment)
    "y_pix",
    "xmin",
    "xmax",
    "ymin",
    "ymax",
    "pa_deg",  # major-axis angle in the pixel frame, degrees
    "ap_radius",  # Kron aperture semi-major axis, arcsec (isophotal radius as fallback)
    "deblend",  # deblended in detection
    "use",  # catalogue says the colour-aperture photometry is usable
    "bcg",
    "pointsrc",
    "z_low",  # Z_SPEC where > 0, else the 2.5 % photo-z percentile
    "z16",
    "z84",
    "z_best",
    "z_spec",
    "mu",  # catalogue magnification (CANUCS lens model); 1 where the catalogue has none
    "flux",
    "err",
)


def _masked_float(col) -> np.ndarray:
    return np.asarray(np.ma.filled(np.ma.asarray(col).astype(float), np.nan), float)


def _clean(f: np.ndarray, e: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    bad = ~np.isfinite(f) | ~np.isfinite(e) | (e <= 0)
    f[bad], e[bad] = np.nan, np.nan
    return f, e


def canucs_bands(cat: Table) -> list[str]:
    """CANUCS bands with colour-aperture fluxes, in catalogue order."""
    pre = "FLUX_COLOR03_TOTAL_"
    return [c[len(pre) :] for c in cat.colnames if c.startswith(pre)]


def canucs_standard(cat: Table) -> Table:
    """CANUCS DR1 photometry catalogue (``FLUX_COLOR03_TOTAL_*``, nJy) -> standard layout."""
    bands = canucs_bands(cat)
    f = np.column_stack([np.asarray(cat[f"FLUX_COLOR03_TOTAL_{b}"], float) for b in bands])
    e = np.column_stack([np.asarray(cat[f"FLUXERR_COLOR03_TOTAL_{b}"], float) for b in bands])
    f, e = _clean(f, e)
    zspec = np.asarray(cat["Z_SPEC"], float)
    has_spec = np.isfinite(zspec) & (zspec > 0)
    kron = KRON_SCALE * np.asarray(cat["KRON_RADIUS"], float) * np.asarray(cat["A"], float)
    iso = np.sqrt(np.clip(np.asarray(cat["AREA_ISO"], float), 0, None) / np.pi)
    out = Table(
        {
            "src_id": np.asarray(cat["SOURCE"]),
            "ra": np.asarray(cat["RA"], float),
            "dec": np.asarray(cat["DEC"], float),
            "x_pix": np.asarray(cat["X"], float),
            "y_pix": np.asarray(cat["Y"], float),
            "xmin": np.asarray(cat["X_MIN"], float),
            "xmax": np.asarray(cat["X_MAX"], float),
            "ymin": np.asarray(cat["Y_MIN"], float),
            "ymax": np.asarray(cat["Y_MAX"], float),
            "pa_deg": np.asarray(cat["PA"], float),
            "ap_radius": np.where(np.isfinite(kron) & (kron > 0), kron, iso) * PIXSCALE,
            "deblend": np.asarray(cat["FLAG_DEBLEND"], bool),
            "use": np.asarray(cat["USE_PHOT_APER03"], bool),
            "bcg": np.asarray(cat["FLAG_BCG"], bool),
            "pointsrc": np.asarray(cat["FLAG_POINTSRC"], bool),
            "z_low": np.where(has_spec, zspec, np.asarray(cat["Z025"], float)),
            "z16": np.asarray(cat["Z160"], float),
            "z84": np.asarray(cat["Z840"], float),
            "z_best": np.asarray(cat["Z_ML"], float),
            "z_spec": zspec,
            "mu": np.asarray(cat["MU"], float),
            "flux": f,
            "err": e,
        }
    )
    out.meta.update(
        bands=bands,
        format="canucs",
        pixscale=PIXSCALE,
        source="CANUCS DR1 photometry catalogue (colour apertures, EAzY photo-z)",
    )
    return out


def dja_bands(cat: Table, aperture: int = DJA_APERTURE) -> list[str]:
    """DJA HST+NIRCam bands (upper case) with aperture fluxes; MIRI and ``<band>u`` duplicates
    of a band that is also present are left out (ASSUMPTION: the same filter, re-imaged)."""
    suffix = f"_flux_aper_{aperture}"
    names = [c[: -len(suffix)] for c in cat.colnames if c.endswith(suffix)]
    names = [n for n in names if re.fullmatch(r"f\d{3,4}(w|m|n|lp)u?", n)]
    keep = [n for n in names if not (n.endswith("u") and n[:-1] in names)]
    return [n.upper() for n in keep if n.upper() not in MIRI_BANDS]


def dja_standard(cat: Table, zout: Table, aperture: int = DJA_APERTURE) -> Table:
    """DJA grizli ``*_phot.fits`` (µJy apertures) plus its eazy-py ``zout`` -> standard layout.

    The zout must be row-aligned with the catalogue (same ids in the same order; checked).
    Fluxes are converted to nJy; SEP aperture flags in ``photometry.BAD_FLAGS`` invalidate a
    band. The DJA apertures are not PSF-matched (D-051): colours of a resolved source drift with
    band, but both members of a lensed pair of compact images drift alike. ``ap_radius`` is
    ``DJA_KRON_PER_R50 x flux_radius`` (the DJA Kron apertures follow another convention)."""
    ids, zids = np.asarray(cat["id"]), np.asarray(zout["id"])
    if len(ids) != len(zids) or not np.array_equal(ids, zids):
        raise ValueError("DJA catalogue and zout are not row-aligned (different versions?)")
    bands = dja_bands(cat, aperture)
    fl, el = [], []
    for b in bands:
        lo = b.lower()
        col = cat[f"{lo}_flux_aper_{aperture}"]
        unit = str(col.unit or "").lower()
        if unit not in ("ujy", "µjy"):
            raise ValueError(f"{lo}_flux_aper_{aperture} has unit {unit!r}, expected uJy")
        f = _masked_float(col) * 1000.0
        e = _masked_float(cat[f"{lo}_fluxerr_aper_{aperture}"]) * 1000.0
        flag_col = f"{lo}_flag_aper_{aperture}"
        if flag_col in cat.colnames:
            flags = np.asarray(np.ma.filled(cat[flag_col], 0), int)
            f[(flags & BAD_FLAGS) != 0] = np.nan
        fl.append(f)
        el.append(e)
    f, e = _clean(np.column_stack(fl), np.column_stack(el))
    pixscale = float(cat.meta["ASEC_0"]) / float(cat.meta["APER_0"])
    # DJA Kron apertures (SEP, kron_radius 2.4-3.8 x a) are ~2.6x the CANUCS ones in units of
    # the half-light radius, so the same_galaxy rule would swallow every pair under ~1.5". Use
    # the half-light radius times the CANUCS ratio instead (D-051).
    kron = DJA_KRON_PER_R50 * _masked_float(cat["flux_radius"])
    iso = np.sqrt(np.clip(_masked_float(cat["area_iso"]), 0, None) / np.pi)
    zspec = _masked_float(zout["z_spec"])
    has_spec = np.isfinite(zspec) & (zspec > 0)
    det_flag = np.asarray(np.ma.filled(cat[f"flag_aper_{aperture}"], 0), int)
    n = len(cat)
    out = Table(
        {
            "src_id": ids,
            "ra": _masked_float(cat["ra"]),
            "dec": _masked_float(cat["dec"]),
            "x_pix": _masked_float(cat["x"]),
            "y_pix": _masked_float(cat["y"]),
            "xmin": _masked_float(cat["xmin"]),
            "xmax": _masked_float(cat["xmax"]),
            "ymin": _masked_float(cat["ymin"]),
            "ymax": _masked_float(cat["ymax"]),
            "pa_deg": np.rad2deg(_masked_float(cat["theta_image"])),
            "ap_radius": np.where(np.isfinite(kron) & (kron > 0), kron, iso) * pixscale,
            "deblend": (np.asarray(np.ma.filled(cat["flag"], 0), int) & DJA_DEBLENDED) != 0,
            "use": (det_flag & BAD_FLAGS) == 0,
            "bcg": np.zeros(n, bool),
            "pointsrc": np.zeros(n, bool),
            "z_low": np.where(has_spec, zspec, _masked_float(zout["z025"])),
            "z16": _masked_float(zout["z160"]),
            "z84": _masked_float(zout["z840"]),
            "z_best": _masked_float(zout["z_ml"]),
            "z_spec": zspec,
            "mu": np.ones(n),
            "flux": f,
            "err": e,
        }
    )
    out.meta.update(
        bands=bands,
        format="dja",
        pixscale=pixscale,
        source=f"DJA grizli catalogue, aperture {aperture} "
        f'({float(cat.meta[f"ASEC_{aperture}"]):.2f}" diameter), eazy-py zout',
    )
    return out


def as_standard(cat: Table) -> Table:
    """The standard layout of ``cat``: returned as is when it already is one, else CANUCS."""
    if "flux" in cat.colnames and "bands" in cat.meta:
        return cat
    if any(c.startswith("FLUX_COLOR03_TOTAL_") for c in cat.colnames):
        return canucs_standard(cat)
    raise ValueError("unknown catalogue layout: pass dja_standard(cat, zout) for DJA catalogues")


def bands_of(cat: Table) -> list[str]:
    """Bands of a catalogue (standard layout, or CANUCS)."""
    return list(as_standard(cat).meta["bands"])


def flux_matrix(cat: Table, bands: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """(N, B) fluxes and errors (nJy) of ``bands``; invalid entries are NaN."""
    std = as_standard(cat)
    idx = [list(std.meta["bands"]).index(b) for b in bands]
    f, e = np.asarray(std["flux"], float), np.asarray(std["err"], float)
    return f[:, idx].copy(), e[:, idx].copy()


def snr_bands_of(cat: Table) -> list[str]:
    """The bands of ``SNR_BANDS`` the catalogue has (CANUCS macs0417/macs1423 NCF lack F356W)."""
    have = bands_of(cat)
    out = [b for b in SNR_BANDS if b in have]
    if not out:
        raise ValueError(f"no band of {SNR_BANDS} in the catalogue")
    return out


def summed_snr(cat: Table) -> np.ndarray:
    """S/N of the summed SNR-band flux; -inf where any of those bands is invalid."""
    fs, es = flux_matrix(cat, snr_bands_of(cat))
    with np.errstate(invalid="ignore"):
        snr = fs.sum(1) / np.sqrt((es**2).sum(1))
    return np.nan_to_num(snr, nan=-np.inf)


def select_sources(cat: Table, z_cluster: float) -> np.ndarray:
    """Boolean mask of step 1 (see module docstring)."""
    std = as_standard(cat)
    # an invalid band (NaN, or error <= 0) makes the summed S/N invalid
    snr = summed_snr(std)
    behind = np.asarray(std["z_low"], float) > z_cluster + Z_MARGIN
    keep = np.asarray(std["use"], bool) & ~np.asarray(std["bcg"], bool)
    f, e = np.asarray(std["flux"], float), np.asarray(std["err"], float)
    with np.errstate(invalid="ignore"):
        n_hi = (f / e >= HI_SNR).sum(1)
    return keep & (snr >= MIN_SNR) & (n_hi >= MIN_HI_BANDS) & behind


def tangent_xy(ra, dec, ra0: float, dec0: float) -> tuple[np.ndarray, np.ndarray]:
    """Flat-sky offsets in arcsec (x towards East, y North); fine over a cluster field."""
    ra, dec = np.asarray(ra, float), np.asarray(dec, float)
    return (ra - ra0) * np.cos(np.deg2rad(dec0)) * 3600.0, (dec - dec0) * 3600.0


def aperture_radius(cat: Table) -> np.ndarray:
    """Kron aperture semi-major axis (arcsec); the isophotal radius where KRON_RADIUS is 0."""
    return np.asarray(as_standard(cat)["ap_radius"], float)


# ------------------------------------------------------------------------------- matching


def sed_chi2(f1, e1, f2, e2, floor: float = ERR_FLOOR, n_iter: int = 5):
    """Chi^2 of SED 1 against ``scale`` x SED 2 over bands valid in both, per row.

    Inputs are (P, B) arrays (NaN = invalid). Returns ``chi2, nband, scale, n_agree`` where
    ``n_agree`` counts bands whose scaled residual is below ``AGREE_SIGMA``."""
    f1, e1, f2, e2 = (np.atleast_2d(np.asarray(a, float)) for a in (f1, e1, f2, e2))
    ok = np.isfinite(f1) & np.isfinite(f2) & np.isfinite(e1) & np.isfinite(e2)
    f1z, f2z = np.where(ok, f1, 0.0), np.where(ok, f2, 0.0)
    v1 = np.where(ok, e1**2 + (floor * f1z) ** 2, 1.0)
    v2 = np.where(ok, e2**2 + (floor * f2z) ** 2, 1.0)
    w = ok / v1
    scale = (w * f1z * f2z).sum(1) / np.maximum((w * f2z**2).sum(1), 1e-300)
    for _ in range(n_iter):  # the error of a*f2 depends on a: iterate the weights
        w = ok / (v1 + scale[:, None] ** 2 * v2)
        scale = (w * f1z * f2z).sum(1) / np.maximum((w * f2z**2).sum(1), 1e-300)
    resid2 = w * (f1z - scale[:, None] * f2z) ** 2
    return resid2.sum(1), ok.sum(1), scale, ((resid2 < AGREE_SIGMA**2) & ok).sum(1)


def match_table(i, j, f, e, z16, z84) -> Table:
    """SED and photo-z comparison of index pairs ``(i, j)``; ``match`` applies step 2's cuts."""
    chi2, nb, scale, n_agree = sed_chi2(f[i], e[i], f[j], e[j])
    dof = np.maximum(nb - 1, 1)
    p = chi2_dist.sf(chi2, dof)
    zover = (z16[i] <= z84[j]) & (z16[j] <= z84[i])
    t = Table({"i": i, "j": j, "chi2": chi2, "nband": nb, "chi2_nu": chi2 / dof, "p_sed": p})
    t["scale"], t["n_agree"], t["z_overlap"] = scale, n_agree, zover
    t["match"] = (nb >= MIN_BANDS) & (p >= P_MIN) & zover
    return t


def close_pairs(x, y, rmin: float, rmax: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Unordered index pairs with rmin <= separation <= rmax (arcsec), and the separations."""
    pts = np.column_stack([x, y])
    ij = cKDTree(pts).query_pairs(rmax, output_type="ndarray")
    if len(ij) == 0:
        return np.zeros(0, int), np.zeros(0, int), np.zeros(0)
    d = np.hypot(*(pts[ij[:, 0]] - pts[ij[:, 1]]).T)
    keep = d >= rmin
    return ij[keep, 0], ij[keep, 1], d[keep]


def shifted_null(x, y, f, e, z16, z84, n_shift: int, rng) -> dict:
    """Null (a): SED-match fraction of cross pairs between the sources and shifted copies."""
    tree = cKDTree(np.column_stack([x, y]))
    n_pairs = n_match = 0
    for _ in range(n_shift):
        r, phi = rng.uniform(*SHIFT_RANGE), rng.uniform(0, 2 * np.pi)
        other = cKDTree(np.column_stack([x + r * np.cos(phi), y + r * np.sin(phi)]))
        sdm = tree.sparse_distance_matrix(other, SEP_MAX, output_type="ndarray")
        sdm = sdm[sdm["v"] >= SEP_MIN]
        if len(sdm):
            m = match_table(sdm["i"], sdm["j"], f, e, z16, z84)
            n_pairs += len(m)
            n_match += int(m["match"].sum())
    return {"n_pairs": n_pairs, "n_match": n_match, "fraction": n_match / max(n_pairs, 1)}


# ------------------------------------------------------------------------------- exclusions


def read_image_list(path: Path) -> Table:
    """``id ra dec ...`` image lists (CANUCS allmultim / Lenstool multim; '#' comments)."""
    ids, ra, dec = [], [], []
    for line in Path(path).read_text(encoding="latin-1").splitlines():
        parts = line.split()
        if not parts or parts[0].startswith("#") or len(parts) < 3:
            continue
        ids.append(parts[0])
        ra.append(float(parts[1]))
        dec.append(float(parts[2]))
    return Table({"image_id": ids, "ra": ra, "dec": dec})


def classify_pairs(pairs: Table, cat: Table, images: Table, ra0: float, dec0: float) -> Table:
    """Add ``near_image``, ``same_galaxy``, ``lens_dist`` and ``pair_class`` (step 3).

    ``pairs`` holds catalogue row indices ``ci``, ``cj``."""
    cat = as_standard(cat)
    x, y = tangent_xy(cat["ra"], cat["dec"], ra0, dec0)
    ci, cj = np.asarray(pairs["ci"], int), np.asarray(pairs["cj"], int)
    if len(images):
        ix, iy = tangent_xy(images["ra"], images["dec"], ra0, dec0)
        itree = cKDTree(np.column_stack([ix, iy]))
        di = itree.query(np.column_stack([x[ci], y[ci]]))[0]
        dj = itree.query(np.column_stack([x[cj], y[cj]]))[0]
        near = np.minimum(di, dj) <= IMAGE_RADIUS
    else:
        near = np.zeros(len(pairs), bool)
    rad = aperture_radius(cat)
    sep = np.asarray(pairs["sep"], float)
    deb = np.asarray(cat["deblend"], bool)

    def overlap(lo, hi):
        a_lo, a_hi = np.asarray(cat[lo], float), np.asarray(cat[hi], float)
        return (a_lo[ci] <= a_hi[cj]) & (a_lo[cj] <= a_hi[ci])

    boxes = overlap("xmin", "xmax") & overlap("ymin", "ymax")
    same = (sep < SAME_GALAXY_FACTOR * (rad[ci] + rad[cj])) | (deb[ci] & deb[cj] & boxes)
    # Lens tests use only catalogued sources (any redshift, BCGs included) at least LENS_RADIUS
    # from both members, so a fragment of a member is not a lens. lens_dist is the distance from
    # the midpoint to the nearest such source (inf: none within 0.5 sep + LENS_RADIUS); line_dist
    # its distance to the segment joining the members (a lens need not sit at the midpoint).
    good = np.isfinite(x) & np.isfinite(y)
    idx = np.flatnonzero(good)
    tree = cKDTree(np.column_stack([x[idx], y[idx]]))
    mx, my = 0.5 * (x[ci] + x[cj]), 0.5 * (y[ci] + y[cj])
    lens_dist = np.full(len(pairs), np.inf)
    line_dist = np.full(len(pairs), np.inf)
    lens_id = np.zeros(len(pairs), int)
    n_between = np.zeros(len(pairs), int)  # sources inside the circle on the pair's diameter
    source = np.asarray(cat["src_id"])
    near_mid = tree.query_ball_point(np.column_stack([mx, my]), 0.5 * sep + LENS_RADIUS)
    for p, cand in enumerate(near_mid):
        cand = idx[np.asarray(cand, int)]
        cand = cand[(cand != ci[p]) & (cand != cj[p])]
        if len(cand) == 0:
            continue
        # every lens test needs the source clear of both members (not a fragment of one)
        clear = np.minimum(
            np.hypot(x[cand] - x[ci[p]], y[cand] - y[ci[p]]),
            np.hypot(x[cand] - x[cj[p]], y[cand] - y[cj[p]]),
        )
        cand = cand[clear >= LENS_RADIUS]
        if len(cand) == 0:
            continue
        dm = np.hypot(x[cand] - mx[p], y[cand] - my[p])
        k = int(np.argmin(dm))
        lens_dist[p] = dm[k]
        # distance to the segment, for sources projecting inside it
        ux, uy = x[cj[p]] - x[ci[p]], y[cj[p]] - y[ci[p]]
        t = ((x[cand] - x[ci[p]]) * ux + (y[cand] - y[ci[p]]) * uy) / (ux * ux + uy * uy)
        dl = np.abs((x[cand] - x[ci[p]]) * uy - (y[cand] - y[ci[p]]) * ux) / np.hypot(ux, uy)
        inside = (t > 0) & (t < 1)
        if inside.any():
            line_dist[p] = dl[inside].min()
        between = dm <= 0.5 * sep[p]
        n_between[p] = int(between.sum())
        # the catalogued lens: nearest qualifying source to the midpoint (0 when none qualifies)
        hit = (dm <= LENS_RADIUS) | between | (inside & (dl <= LENS_RADIUS))
        if hit.any():
            lens_id[p] = source[cand[hit][int(np.argmin(dm[hit]))]]
    pairs["near_image"], pairs["same_galaxy"] = near, same
    pairs["lens_dist"], pairs["lens_source"], pairs["line_dist"] = lens_dist, lens_id, line_dist
    pairs["n_between"] = n_between
    cls = np.where(near, "published", np.where(same, "same_galaxy", "orphan"))
    lens = (lens_dist <= LENS_RADIUS) | (line_dist <= LENS_RADIUS) | (n_between > 0)
    cls = np.where((cls == "orphan") & lens, "visible_lens", cls)
    pairs["pair_class"] = cls
    return pairs


def mirror_angle(cat: Table, ci, cj) -> np.ndarray:
    """|dphi_i + dphi_j| (deg, 0-90), dphi = major-axis angle minus the joining-line angle.

    A fold image pair is mirror-symmetric across the line's perpendicular bisector, so the two
    members' angles to the joining line are opposite (value near 0). Pixel-frame angles."""
    cat = as_standard(cat)
    X, Y, pa = (np.asarray(cat[c], float) for c in ("x_pix", "y_pix", "pa_deg"))
    line = np.rad2deg(np.arctan2(Y[cj] - Y[ci], X[cj] - X[ci]))
    s = np.mod(pa[ci] - line + pa[cj] - line, 180.0)
    return np.minimum(s, 180.0 - s)


# ------------------------------------------------------------------------------- dark-lens numbers


def dark_lens_numbers(sep_arcsec, z_l: float, z_s) -> Table:
    """HYPOTHESIS numbers for a point/SIS deflector at ``z_l`` making the pair.

    theta_E = sep / 2; M_E = theta_E^2 c^2 D_l D_s / (4 G D_ls) (mass inside theta_E); SIS
    sigma = c sqrt(theta_E D_s / (4 pi D_ls)). Cosmology: flat LCDM H0=70, Om=0.3 (ASSUMPTION)."""
    import astropy.units as u
    from astropy.constants import G, c
    from astropy.cosmology import FlatLambdaCDM

    cosmo = FlatLambdaCDM(H0=70.0, Om0=0.3)
    theta = (np.asarray(sep_arcsec, float) / 2.0 * u.arcsec).to_value(u.rad)
    z_s = np.asarray(z_s, float)
    d_l = cosmo.angular_diameter_distance(z_l)
    d_s = cosmo.angular_diameter_distance(z_s)
    d_ls = cosmo.angular_diameter_distance(z_l, z_s)
    m_e = (theta**2 * c**2 * d_l * d_s / (4 * G * d_ls)).to_value(u.Msun)
    sigma = (c * np.sqrt(theta * d_s / (4 * np.pi * d_ls))).to_value(u.km / u.s)
    return Table({"theta_e": np.asarray(sep_arcsec) / 2.0, "mass_e_msun": m_e, "sis_sigma": sigma})


# ------------------------------------------------------------------------------- per field


def search(cat: Table, z_cluster: float, images: Table, n_shift: int = 20, seed: int = 1):
    """Steps 1-4 on one catalogue: ``(pairs, summary)``; ``pairs`` holds SED-matched close pairs.

    ``cat`` is a standard-layout table (``as_standard``; a CANUCS table is converted);
    ``z_cluster`` is the deflector redshift of the background cut (deep fields: Z_LENS_REF)."""
    cat = as_standard(cat)
    sel = select_sources(cat, z_cluster)
    rows = np.flatnonzero(sel)
    sub = cat[rows]
    ra0, dec0 = float(np.nanmedian(cat["ra"])), float(np.nanmedian(cat["dec"]))
    x, y = tangent_xy(sub["ra"], sub["dec"], ra0, dec0)
    f, e = np.asarray(sub["flux"], float), np.asarray(sub["err"], float)
    z16, z84 = np.asarray(sub["z16"], float), np.asarray(sub["z84"], float)

    i, j, d = close_pairs(x, y, SEP_MIN, SEP_MAX)
    m = match_table(i, j, f, e, z16, z84)
    m["sep"] = d
    fi, fj, fd = close_pairs(x, y, FAR_MIN, FAR_MAX)
    far = match_table(fi, fj, f, e, z16, z84)
    null_a = shifted_null(x, y, f, e, z16, z84, n_shift, np.random.default_rng(seed))
    frac_b = float(far["match"].mean()) if len(far) else float("nan")
    zf = far[far["z_overlap"]]
    frac_c = float(zf["match"].mean()) if len(zf) else float("nan")

    m["ci"], m["cj"] = rows[m["i"]], rows[m["j"]]
    m = classify_pairs(m, cat, images, ra0, dec0)  # geometry only: also for unmatched pairs
    pairs = m[m["match"]]
    for k in ("ci", "cj"):
        tag = "a" if k == "ci" else "b"
        cols = (("src_id", "id"), ("ra", "ra"), ("dec", "dec"), ("z_best", "z"))
        cols += (("z_spec", "zspec"), ("mu", "mu_cat"), ("pointsrc", "pointsrc"))
        for col, out in cols:
            pairs[f"{out}_{tag}"] = np.asarray(cat[col])[pairs[k]]
    pairs["mirror_deg"] = mirror_angle(cat, pairs["ci"], pairs["cj"])
    pairs.meta.update(provenance=schema.Provenance.DERIVED.value, source=cat.meta["source"])
    n_close = len(m)
    counts = {c: int((pairs["pair_class"] == c).sum()) for c in CLASSES}
    summary = {
        "n_catalogue": len(cat),
        "n_selected": int(sel.sum()),
        "n_close_pairs": n_close,
        "n_close_testable": int((m["nband"] >= MIN_BANDS).sum()),
        "n_sed_matched": len(pairs),
        "null_a_shift": null_a,
        "null_a_expected": null_a["fraction"] * n_close,
        "null_b_far": {"n_pairs": len(far), "n_match": int(far["match"].sum()), "fraction": frac_b},
        "null_b_expected": frac_b * n_close,
        # (c) redshift clustering: SED-match rate among photo-z-overlapping far pairs
        "null_c_zmatched": {
            "far_fraction": frac_c,
            "n_close_zoverlap": int(m["z_overlap"].sum()),
            "expected": frac_c * int(m["z_overlap"].sum()),
        },
        "classes": counts,
        "expected_by_class": {
            c: float(np.sum(m["pair_class"] == c) * null_a["fraction"]) for c in CLASSES
        },
        "expected_by_class_zclustered": {
            c: float(np.sum((m["pair_class"] == c) & m["z_overlap"]) * frac_c) for c in CLASSES
        },
        "poisson_p_orphan_excess": poisson_excess(
            counts.get("orphan", 0), float(np.sum(m["pair_class"] == "orphan") * null_a["fraction"])
        ),
        "poisson_p_orphan_excess_zclustered": poisson_excess(
            counts.get("orphan", 0),
            float(np.sum((m["pair_class"] == "orphan") & m["z_overlap"]) * frac_c),
        ),
        "matched_fraction_by_sep": sep_bins(m),
        "bands": list(cat.meta["bands"]),
        "snr_bands": snr_bands_of(cat),
        "z_lens": z_cluster,
    }
    # (d) small-scale clustering (D-051): the match rate of photo-z-overlapping pairs at
    # NEAR_MIN-NEAR_MAX, where physical companions are commoner than at 10-30" but blending and
    # lensing by galaxies are not an issue. Its z-overlap match rate times the z-overlapping
    # close pairs per class; and the z-overlap match rate per separation bin out to FAR_MIN.
    ni, nj, nd = close_pairs(x, y, NEAR_MIN, NEAR_MAX)
    near = match_table(ni, nj, f, e, z16, z84)
    near["sep"] = nd
    zn = near[near["z_overlap"]]
    frac_d = float(zn["match"].mean()) if len(zn) else float("nan")
    exp_d = {c: float(np.sum((m["pair_class"] == c) & m["z_overlap"]) * frac_d) for c in CLASSES}
    summary["null_d_near_zmatched"] = {
        "sep": [NEAR_MIN, NEAR_MAX],
        "n_pairs_zoverlap": len(zn),
        "fraction": frac_d,
        "expected_by_class": exp_d,
        "poisson_p_orphan_excess": poisson_excess(counts.get("orphan", 0), exp_d["orphan"]),
    }
    # (e) the z-clustered null (c) conditioned on what makes two SEDs easy to match (D-051):
    # the far z-overlap match rate in cells of the pair's fainter-member S/N, larger member's
    # aperture radius and LW/SW colour, summed over each class's z-overlapping close pairs.
    key_c, key_f = pair_cells(sub, m), pair_cells(sub, far)
    zfar = np.asarray(far["z_overlap"], bool)
    rate = {k: float(np.mean(far["match"][zfar & (key_f == k)])) for k in np.unique(key_f[zfar])}
    exp_cell = np.array([rate.get(k, frac_c) for k in key_c])  # empty cell: the global rate
    zc = np.asarray(m["z_overlap"], bool)
    exp_e = {c: float(exp_cell[zc & (m["pair_class"] == c)].sum()) for c in CLASSES}
    summary["null_e_conditioned"] = {
        "cells": "fainter-member summed S/N x larger aperture radius x LW/SW colour",
        "expected_by_class": exp_e,
        "poisson_p_orphan_excess": poisson_excess(counts.get("orphan", 0), exp_e["orphan"]),
    }
    mid = close_pairs(x, y, NEAR_MAX, FAR_MIN)
    mid_t = match_table(*mid[:2], f, e, z16, z84)
    mid_t["sep"] = mid[2]
    zall = [t[t["z_overlap"]] for t in (m, near, mid_t)]
    summary["zoverlap_match_fraction_by_sep"] = sep_bins(
        vstack_rows(zall), edges=(0.3, 0.6, 1.0, 1.5, 2.0, 3.0, 4.5, 6.0, 8.0, 10.0)
    )
    return pairs, summary


CELL_SNR = (20.0, 40.0, 100.0)  # fainter member's summed S/N bin edges (ASSUMPTION)
CELL_RADIUS = (0.4, 0.6, 0.9)  # larger member's aperture radius, arcsec (ASSUMPTION)
CELL_COLOUR = (0.0, 0.3, 0.6)  # log10(LW / SW flux) of member i (ASSUMPTION)
SW_BANDS = ("F090W", "F115W", "F150W")


def pair_cells(sub: Table, pairs: Table) -> np.ndarray:
    """Integer cell of each pair for null (e): S/N x size x colour bins (see CELL_*)."""
    snr = summed_snr(sub)
    rad = np.asarray(sub["ap_radius"], float)
    lw, _ = flux_matrix(sub, snr_bands_of(sub))
    sw_have = [b for b in SW_BANDS if b in sub.meta["bands"]]
    sw = flux_matrix(sub, sw_have)[0] if sw_have else np.full((len(sub), 1), np.nan)
    import warnings

    with np.errstate(invalid="ignore", divide="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN rows: colour NaN -> 0.0
        colour = np.log10(np.nanmean(lw, 1) / np.nanmean(sw, 1))
    i, j = np.asarray(pairs["i"], int), np.asarray(pairs["j"], int)
    s = np.digitize(np.minimum(snr[i], snr[j]), CELL_SNR)
    r = np.digitize(np.fmax(rad[i], rad[j]), CELL_RADIUS)
    c = np.digitize(np.nan_to_num(colour[i], nan=0.0), CELL_COLOUR)
    return s * 100 + r * 10 + c


def vstack_rows(tables: list[Table]) -> Table:
    """``sep`` and ``match`` columns of several pair tables, stacked."""
    return Table(
        {
            "sep": np.concatenate([np.asarray(t["sep"], float) for t in tables]),
            "match": np.concatenate([np.asarray(t["match"], bool) for t in tables]),
        }
    )


CLASSES = ("published", "same_galaxy", "visible_lens", "orphan")


def poisson_excess(observed: int, expected: float) -> float:
    """P(N >= observed | Poisson(expected)): chance of at least this many matched pairs."""
    from scipy.stats import poisson

    return float(poisson.sf(observed - 1, expected)) if expected > 0 else float("nan")


def sep_bins(m: Table, edges=(0.3, 0.6, 1.0, 1.5, 2.0, 3.0)) -> list[dict]:
    """SED-match fraction of close pairs per separation bin (blending raises it at small sep)."""
    out = []
    for lo, hi in zip(edges[:-1], edges[1:], strict=True):
        s = (m["sep"] >= lo) & (m["sep"] < hi)
        frac = float(m["match"][s].mean()) if s.any() else float("nan")
        out.append({"sep": [lo, hi], "n": int(s.sum()), "frac": frac})
    return out


def rank_orphans(pairs: Table) -> Table:
    """Orphans ordered by agreeing bands (desc), then reduced chi^2 (asc)."""
    o = pairs[pairs["pair_class"] == "orphan"]
    if len(o) == 0:
        return o
    o = o[np.lexsort((np.asarray(o["chi2_nu"]), -np.asarray(o["n_agree"])))]
    o["rank"] = np.arange(1, len(o) + 1)
    return o


# ------------------------------------------------------------------------------- follow-up


def lens_check(top: Table, model_name: str | None, z_cluster: float) -> Table:
    """CATS map magnification (signed) at both members at their photo-z (model_prediction).

    ``model_name`` None (deep fields): no lens model is needed; the catalogue magnification
    (CANUCS ``MU``, ~1 in the flanking fields; 1 for DJA) is the only lensing at the pair, so
    ``cluster_explains`` is ``|mu_cat| > MAX_ORDINARY_MU``."""
    if model_name is None:
        return _no_model_check(top, z_cluster)
    import lens_consistency as lc

    model, _, _ = lc.load_model(model_name)
    offset = lc.apply_frame_offset(model_name, model)
    zs = np.maximum(
        0.5 * (np.asarray(top["z_a"], float) + np.asarray(top["z_b"], float)), z_cluster + 0.05
    )
    for tag in ("a", "b"):
        mu = model.magnification(top[f"ra_{tag}"], top[f"dec_{tag}"], zs)
        top[f"mu_cats_{tag}"] = np.asarray(mu["magnification"], float)
    top["z_pair"] = zs
    mu_a, mu_b = np.asarray(top["mu_cats_a"]), np.asarray(top["mu_cats_b"])
    known = np.isfinite(mu_a) & np.isfinite(mu_b)  # NaN: outside the model's maps
    top["parity_flip"] = known & (np.sign(mu_a) != np.sign(mu_b))
    big = np.fmax(np.abs(mu_a), np.abs(mu_b))
    # outside the CATS maps: fall back to the CANUCS catalogue |mu| (its own lens model)
    canucs = np.fmax(
        np.abs(np.asarray(top["mu_cat_a"], float)), np.abs(np.asarray(top["mu_cat_b"], float))
    )
    big = np.where(np.isfinite(big), big, canucs)
    top["model_covered"] = known
    top["cluster_explains"] = top["parity_flip"] | (big > MAX_ORDINARY_MU)
    dl = dark_lens_numbers(top["sep"], z_cluster, zs)
    for c in dl.colnames:
        top[c] = dl[c]
    top.meta["source"] = (
        f"{top.meta.get('source', 'CANUCS DR1 photometry catalogue')}; magnifications from "
        f"{model.source} (frame offset {offset[0]:+.3f}, {offset[1]:+.3f} arcsec) and the "
        "CANUCS catalogue MU"
    )
    top.meta["model_prediction_columns"] = [
        "mu_cats_a",
        "mu_cats_b",
        "z_pair",
        "parity_flip",
        "model_covered",
        "cluster_explains",
    ]
    top.meta["hypothesis_columns"] = list(dl.colnames)
    return top


def _no_model_check(top: Table, z_lens: float) -> Table:
    """Deep-field version of ``lens_check``: |mu| from the catalogue, no parity information."""
    zs = np.maximum(
        0.5 * (np.asarray(top["z_a"], float) + np.asarray(top["z_b"], float)), z_lens + 0.05
    )
    top["z_pair"] = zs
    big = np.fmax(np.abs(np.asarray(top["mu_cat_a"], float)), np.abs(np.asarray(top["mu_cat_b"])))
    top["parity_flip"] = np.zeros(len(top), bool)
    top["model_covered"] = np.zeros(len(top), bool)
    top["cluster_explains"] = big > MAX_ORDINARY_MU
    dl = dark_lens_numbers(top["sep"], z_lens, zs)
    for c in dl.colnames:
        top[c] = dl[c]
    top.meta["source"] = (
        f"{top.meta.get('source', '')}; no cluster lens model (deep field): catalogue mu only; "
        f"dark-lens numbers at the notional z_l = {z_lens}"
    )
    top.meta["model_prediction_columns"] = ["z_pair", "cluster_explains"]
    top.meta["hypothesis_columns"] = list(dl.colnames)
    return top


def cutouts_and_sheet(top: Table, field: str, out_dir: Path, title: str) -> Path:
    """F150W/F277W/F444W cutouts centred on each pair's midpoint, members circled."""
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.wcs.utils import proj_plane_pixel_scales
    from matplotlib.figure import Figure
    from matplotlib.patches import Circle

    from jwst_anomaly import cutouts, viz

    tgt = Table(
        {
            "source_uid": [f"r{r['rank']}_{r['id_a']}_{r['id_b']}" for r in top],
            "ra": 0.5 * (np.asarray(top["ra_a"]) + np.asarray(top["ra_b"])),
            "dec": 0.5 * (np.asarray(top["dec_a"]) + np.asarray(top["dec_b"])),
        }
    )
    spec = FIELDS[field]
    if spec.get("deep"):
        uris = deep_i2d_uris(spec["mast"], tgt)
    else:
        i2d_prefix = spec["i2d"]
        man = Table.read(paths.manifests_dir() / f"{field}_products.ecsv")
        names = [str(n) for n in man["productFilename"]]
        uris = {
            b: [str(man["cloud_uri"][names.index(f"{i2d_prefix}_nircam_clear-{b}_i2d.fits")])]
            * len(tgt)
            for b in CUTOUT_BANDS
        }
    cut: dict[str, list[str]] = {}
    for b, per_target in uris.items():  # one byte-range reader per i2d file
        cut[b] = [""] * len(tgt)
        for uri in sorted(set(per_target) - {""}):
            k = np.flatnonzero(np.asarray(per_target) == uri)
            res = _with_retries(
                lambda uri=uri, k=k: cutouts.make_cutouts(
                    uri, tgt[k], CUTOUT_ARCSEC, out_dir / "cutouts"
                )
            )
            for kk, path in zip(k, res["path"], strict=True):
                cut[b][kk] = str(path)
    top.meta["i2d_files"] = {
        b: sorted({u.rsplit("/", 1)[-1] for u in v if u}) for b, v in uris.items()
    }
    per_row = 3  # pairs per sheet row, each pair one panel per band
    nrow = -(-len(tgt) // per_row)
    ncol = per_row * len(cut)
    fig = Figure(figsize=(2.0 * ncol, 2.3 * nrow + 0.6), layout="constrained")
    axes = fig.subplots(nrow, ncol, squeeze=False)
    for ax in axes.ravel():
        ax.set_axis_off()
    for r, row in enumerate(top):
        for c, b in enumerate(cut):
            ax = axes[r // per_row, (r % per_row) * len(cut) + c]
            ax.set_axis_on()
            ax.set_xticks([])
            ax.set_yticks([])
            path = cut[b][r]
            if not path:
                ax.set_title(f"{tgt['source_uid'][r]} {b}: outside", fontsize=7)
                continue
            with fits.open(path) as h:
                data, w = h["SCI"].data, WCS(h["SCI"].header)
            ax.imshow(data, norm=viz.normalize(data), origin="lower", cmap="gray")
            r_px = 0.15 / (proj_plane_pixel_scales(w)[0] * 3600.0)  # 0.15" circles
            for tag, col in (("a", "tab:red"), ("b", "tab:cyan")):
                px, py = w.world_to_pixel_values(row[f"ra_{tag}"], row[f"dec_{tag}"])
                ax.add_patch(Circle((px, py), r_px, fill=False, ec=col, lw=0.8))
            ax.set_title(
                f"#{row['rank']} {row['id_a']}/{row['id_b']} {b}\n"
                f'sep {row["sep"]:.2f}" chi2nu {row["chi2_nu"]:.2f} '
                f"z {row['z_a']:.2f}/{row['z_b']:.2f}",
                fontsize=6,
            )
    fig.suptitle(title, fontsize=9)
    out = out_dir / "contact_sheet.png"
    fig.savefig(out, dpi=110)
    return out


def _with_retries(fn, attempts: int = 4, wait: float = 10.0):
    """Call ``fn``; retry on OSError (S3 range reads through a proxy fail intermittently with
    spurious NoSuchBucket / FileNotFoundError), waiting ``wait`` x attempt seconds."""
    import time

    for attempt in range(1, attempts + 1):
        try:
            return fn()
        except OSError:
            if attempt == attempts:
                raise
            time.sleep(wait * attempt)


def _in_region(s_region: str, ra: float, dec: float) -> bool:
    """Whether (ra, dec) lies inside a MAST ``s_region`` (one or more ``POLYGON``s)."""
    from matplotlib.path import Path as MplPath

    for poly in re.findall(r"POLYGON\s+(?:[A-Z]\w*\s+)?([-0-9.E\s]+)", s_region.upper()):
        v = np.array(poly.split(), float).reshape(-1, 2)
        x = (v[:, 0] - ra + 180.0) % 360.0 - 180.0  # wrap-safe RA offsets
        x *= np.cos(np.deg2rad(dec))
        if MplPath(np.column_stack([x, v[:, 1] - dec])).contains_point((0.0, 0.0)):
            return True
    return False


def deep_i2d_uris(mast: tuple[str, list[str]], tgt: Table) -> dict[str, list[str]]:
    """Per band, the S3 URI of a level-3 ``_i2d`` covering each target ("" if none).

    One MAST observation query and one product listing per band (batched over targets); the
    observation is the first, by obs_id, whose footprint contains the target."""
    from jwst_anomaly import query

    program, prefixes = mast
    out = {}
    for b in CUTOUT_BANDS:
        obs = query.query_observations(
            proposal_id=program, instrument_name="NIRCAM/IMAGE", filters=b.upper(), calib_level=3
        )
        ids = query.str_values(obs["obs_id"])
        obs = obs[[any(i.startswith(p) for p in prefixes) for i in ids]]
        obs.sort("obs_id")  # "first by obs_id", independent of the archive's row order
        chosen = []
        for r in tgt:
            hit = [
                o for o in obs if _in_region(str(o["s_region"]), float(r["ra"]), float(r["dec"]))
            ]
            chosen.append(str(hit[0]["obs_id"]) if hit else "")
        need = obs[[str(i) in set(chosen) for i in obs["obs_id"]]]
        prods = query.list_products(need, subgroups=("I2D",), cloud_uris=True) if len(need) else []
        by_obs = {}
        for p in prods:
            name = str(p["productFilename"])
            if name.endswith(f"_nircam_clear-{b}_i2d.fits") and str(p["cloud_uri"]):
                by_obs[name.split("_nircam_")[0]] = str(p["cloud_uri"])
        out[b] = [by_obs.get(c.split("_nircam_")[0], "") if c else "" for c in chosen]
    return out


def fetch_tar_member(spec: dict, cache_dir: Path | None = None) -> Path:
    """One member of a remote ``.tar.gz``, verified; the archive itself is never stored.

    ``spec``: ``url``, ``sha256`` (archive), ``member``, ``member_sha256``, ``max_bytes``. The
    archive is streamed once: its sha256 is computed on the fly over every byte and the member
    is kept only when both checksums match (the cache name carries the member's sha256 prefix,
    as ``fetch_catalog`` does)."""
    cache = cache_dir or paths.cache_dir() / "external"
    cache.mkdir(parents=True, exist_ok=True)
    msha = spec["member_sha256"].lower()
    target = cache / f"{msha[:12]}_{Path(spec['member']).name}"
    if target.exists():
        h = hashlib.sha256()
        with open(target, "rb") as fh:
            while chunk := fh.read(1 << 20):
                h.update(chunk)
        if h.hexdigest() == msha:
            return target
        target.unlink()

    class _Tee:  # hashes and counts every byte the tar reader pulls
        def __init__(self, raw):
            self.raw, self.h, self.n = raw, hashlib.sha256(), 0

        def read(self, k=-1):
            chunk = self.raw.read(k)
            self.h.update(chunk)
            self.n += len(chunk)
            if self.n > spec["max_bytes"]:
                raise ValueError(f"{spec['url']}: more than max_bytes {spec['max_bytes']}")
            return chunk

    fd, tmp = tempfile.mkstemp(dir=cache, suffix=".part")
    os.close(fd)
    try:
        mh, found = hashlib.sha256(), False
        with urllib.request.urlopen(spec["url"], timeout=120) as resp:  # noqa: S310
            tee = _Tee(resp)
            with tarfile.open(fileobj=tee, mode="r|gz") as tar:
                for m in tar:
                    if m.name != spec["member"]:
                        continue
                    found = True
                    with open(tmp, "wb") as out, tar.extractfile(m) as src:
                        while chunk := src.read(1 << 20):
                            mh.update(chunk)
                            out.write(chunk)
            while tee.read(1 << 20):  # hash the rest of the archive
                pass
        if not found:
            raise ValueError(f"{spec['url']}: no member {spec['member']}")
        if tee.h.hexdigest() != spec["sha256"].lower():
            raise ValueError(f"{spec['url']}: archive sha256 {tee.h.hexdigest()} does not match")
        if mh.hexdigest() != msha:
            raise ValueError(f"{spec['member']}: sha256 {mh.hexdigest()} does not match")
        os.replace(tmp, target)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)
    return target


def load_catalogue(field: str) -> Table:
    """The field's catalogue in the standard layout (downloaded once, sha256-verified)."""
    spec = FIELDS[field]
    url, sha = spec["catalogue"]
    path = fetch_catalog(url, sha, max_bytes=spec.get("max_bytes", 200_000_000))
    if spec.get("format", "canucs") == "canucs":
        return canucs_standard(Table.read(path))
    zcols = ["id", "z_spec", "z_ml", "z025", "z160", "z840"]
    zout = Table.read(fetch_tar_member(spec["photoz"]))[zcols]
    std = dja_standard(Table.read(path), zout)
    std.meta["source"] += f" ({Path(url).name}; {spec['photoz']['member']})"
    return std


def searched_footprint(cat: Table, step: float = 1.0, radius: float = 4.0):
    """Grid points (arcsec offsets from the catalogue median) of the searched footprint.

    A point is searched when a catalogue row with valid photometry in every S/N band and in at
    least ``MIN_BANDS`` bands lies within ``radius`` arcsec (ASSUMPTION; the D-049 convention).
    Masked regions (bright stars, edges, no coverage) have no such rows. Returns
    ``(gx, gy, area_arcsec2, (ra0, dec0))``."""
    std = as_standard(cat)
    ra0, dec0 = float(np.nanmedian(std["ra"])), float(np.nanmedian(std["dec"]))
    x, y = tangent_xy(std["ra"], std["dec"], ra0, dec0)
    fs, _ = flux_matrix(std, snr_bands_of(std))
    nvalid = np.isfinite(np.asarray(std["flux"], float)).sum(1)
    ok = np.isfinite(fs).all(1) & (nvalid >= MIN_BANDS) & np.isfinite(x) & np.isfinite(y)
    tree = cKDTree(np.column_stack([x[ok], y[ok]]))
    gx, gy = np.meshgrid(
        np.arange(np.nanmin(x[ok]), np.nanmax(x[ok]) + step, step),
        np.arange(np.nanmin(y[ok]), np.nanmax(y[ok]) + step, step),
    )
    gx, gy = gx.ravel(), gy.ravel()
    near = tree.query(np.column_stack([gx, gy]), distance_upper_bound=radius)[0] <= radius
    return gx[near], gy[near], float(near.sum() * step * step), (ra0, dec0)


def load_images(field: str) -> Table:
    """All published multiple images for ``field``: CATS arcs.txt (JWST frame) + CANUCS lists.

    Deep fields have none (an empty table)."""
    import lens_consistency as lc

    spec = FIELDS[field]
    if spec.get("deep"):
        out = Table({"image_id": [], "ra": [], "dec": []}, dtype=[str, float, float])
        out.meta["n_per_list"] = []
        return out
    lists = []
    cats, _ = lc.image_list(spec["cats"], {}, None)
    lc.shift_images(spec["cats"], cats)
    lists.append(Table({"image_id": cats["image_id"], "ra": cats["ra"], "dec": cats["dec"]}))
    for url, sha in spec["image_lists"]:
        lists.append(read_image_list(fetch_catalog(url, sha)))
    for name in spec.get("image_models", []):
        t, _ = lc.image_list(name, {}, None)
        lc.shift_images(name, t)
        lists.append(Table({"image_id": t["image_id"], "ra": t["ra"], "dec": t["dec"]}))
    from astropy.table import vstack

    out = vstack(lists)
    out.meta["n_per_list"] = [len(t) for t in lists]
    return out


def _finite(obj):
    """NaN/inf -> None, recursively, so summary.json is strict JSON."""
    if obj is np.ma.masked:
        return None
    if isinstance(obj, dict):
        return {k: _finite(v) for k, v in obj.items()}
    if isinstance(obj, np.ndarray):
        obj = np.ma.filled(np.ma.asarray(obj).astype(object), None).tolist()
    if isinstance(obj, (list, tuple)):
        return [_finite(v) for v in obj]
    if isinstance(obj, np.generic):
        obj = obj.item()
    if isinstance(obj, float) and not np.isfinite(obj):
        return None
    return obj


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--field", choices=sorted(FIELDS), required=True)
    ap.add_argument("--out", type=Path, default=paths.outputs_dir() / "orphan_pairs")
    ap.add_argument("--n-shift", type=int, default=20)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--top", type=int, default=15)
    ap.add_argument("--cutouts", action="store_true", help="S3 cutouts + contact sheet of the top")
    args = ap.parse_args(argv)
    spec = FIELDS[args.field]
    out = args.out / args.field
    out.mkdir(parents=True, exist_ok=True)
    # never leave an earlier run's outputs next to this run's
    for stale in ("summary.json", "matched_pairs.ecsv", "top_orphans.ecsv", "contact_sheet.png"):
        (out / stale).unlink(missing_ok=True)
    if (out / "cutouts").exists():
        shutil.rmtree(out / "cutouts")  # fail loudly (e.g. a locked file on Windows)
    cat = load_catalogue(args.field)
    images = load_images(args.field)
    z_lens = z_lens_of(spec)
    pairs, summary = search(cat, z_lens, images, args.n_shift, args.seed)
    summary["n_published_images"] = images.meta["n_per_list"]
    if spec.get("deep"):
        summary["footprint_arcmin2"] = searched_footprint(cat)[2] / 3600.0
    pairs.write(out / "matched_pairs.ecsv", overwrite=True)
    top = rank_orphans(pairs)[: args.top]
    if len(top):
        top = lens_check(top, spec["cats"], z_lens)
        if args.cutouts:
            summary["contact_sheet"] = str(
                cutouts_and_sheet(top, args.field, out, f"{args.field} orphan pairs (unvetted)")
            )
            summary["i2d_files"] = top.meta["i2d_files"]
        top.write(out / "top_orphans.ecsv", overwrite=True)
    summary["top"] = [
        {k: _finite(row[k]) for k in top.colnames if k not in ("i", "j", "ci", "cj")} for row in top
    ]
    (out / "summary.json").write_text(json.dumps(_finite(summary), indent=1, default=str))
    print(
        json.dumps(_finite({k: v for k, v in summary.items() if k != "top"}), indent=1, default=str)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
