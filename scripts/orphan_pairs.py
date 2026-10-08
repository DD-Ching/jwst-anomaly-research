"""Blind search for "orphan image pairs" behind a lensing cluster (docs/orphan_pairs.md).

An orphan pair is two background sources 0.3-3" apart whose SEDs and photo-z are statistically
identical, with no catalogued galaxy between them and no published multiple image nearby. It is a
*possible* signature of an unseen compact deflector (a "dark lens") - a hypothesis, never a
conclusion. The search is lens-model independent; a lens model is used afterwards only to test
the ordinary explanation (a cluster-scale image pair near a critical curve).

Steps (every threshold is an ASSUMPTION; defaults below):

1. Sources: ``USE_PHOT_APER03``, not ``FLAG_BCG``, S/N >= ``MIN_SNR`` in the summed F277W+F356W+
   F444W colour-aperture flux, S/N >= ``HI_SNR`` in at least ``MIN_HI_BANDS`` bands (without it
   the chi^2 test has no power: 25 % of random faint pairs "match"), and behind the cluster
   (``Z025 > z_cluster + Z_MARGIN``, or ``Z_SPEC > z_cluster + Z_MARGIN`` when a spectroscopic
   redshift exists).
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
   clustering (physical neighbours share SEDs).
5. ``--cutouts``: F150W/F277W/F444W cutouts of the top orphans (S3 byte ranges) on one contact
   sheet, plus a lens-model check (CATS map magnification and parity at both members).

Outputs are ``derived`` (from observed catalogues) with ``model_prediction`` columns from the lens
model; results go to ``outputs/orphan_pairs/<field>/``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table
from scipy.spatial import cKDTree
from scipy.stats import chi2 as chi2_dist

sys.path.insert(0, str(Path(__file__).resolve().parent))

from jwst_anomaly import paths, schema  # noqa: E402

SEP_MIN, SEP_MAX = 0.3, 3.0  # arcsec, the orphan-pair annulus (task definition)
FAR_MIN, FAR_MAX = 10.0, 30.0  # arcsec, null annulus (b) (ASSUMPTION)
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


def canucs_dir() -> Path:
    return paths.cache_dir() / "external" / "canucs"


_CAT = "hlsp_canucs_jwst-hst_multi_{f}-clu_multi_v1_photometry-cat.fits.gz"
FIELDS = {
    "macs0416": {
        "canucs": "macs0416",
        "z_cluster": 0.396,
        "cats": "macs0416-cats",
        "image_lists": [
            "macs0416/hlsp_canucs_jwst-hst_multi_macs0416-allmultim-cat_multi_v1_model.txt"
        ],
        "i2d": "jw01208-o004_t002",
    },
    "macs1149": {
        "canucs": "macs1149",
        "z_cluster": 0.543,
        "cats": "macs1149-cats",
        "image_lists": [],  # CANUCS DR1 publishes no MACS1149 image list (readme: with v2)
        "i2d": "jw01208-o008_t004",
    },
    "abell370": {
        "canucs": "a370",
        "z_cluster": 0.375,
        "cats": "abell370-cats",
        "image_lists": ["a370/hlsp_canucs_jwst-hst_multi_a370-lenstool-multim_multi_v1_model.txt"],
        "i2d": "jw01208-o002_t001",
    },
}


# ------------------------------------------------------------------------------- catalogue


def bands_of(cat: Table) -> list[str]:
    """Bands with colour-aperture fluxes, in catalogue order."""
    pre = "FLUX_COLOR03_TOTAL_"
    return [c[len(pre) :] for c in cat.colnames if c.startswith(pre)]


def flux_matrix(cat: Table, bands: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """(N, B) colour-aperture fluxes and errors (nJy); invalid entries are NaN."""
    f = np.column_stack([np.asarray(cat[f"FLUX_COLOR03_TOTAL_{b}"], float) for b in bands])
    e = np.column_stack([np.asarray(cat[f"FLUXERR_COLOR03_TOTAL_{b}"], float) for b in bands])
    bad = ~np.isfinite(f) | ~np.isfinite(e) | (e <= 0)
    f[bad], e[bad] = np.nan, np.nan
    return f, e


def select_sources(cat: Table, z_cluster: float) -> np.ndarray:
    """Boolean mask of step 1 (see module docstring)."""
    num = sum(np.nan_to_num(np.asarray(cat[f"FLUX_COLOR03_TOTAL_{b}"], float)) for b in SNR_BANDS)
    var = sum(
        np.nan_to_num(np.asarray(cat[f"FLUXERR_COLOR03_TOTAL_{b}"], float), nan=np.inf) ** 2
        for b in SNR_BANDS
    )
    snr = num / np.sqrt(var)
    zspec = np.asarray(cat["Z_SPEC"], float)
    has_spec = np.isfinite(zspec) & (zspec > 0)
    z_low = np.where(has_spec, zspec, np.asarray(cat["Z025"], float))
    behind = z_low > z_cluster + Z_MARGIN
    keep = np.asarray(cat["USE_PHOT_APER03"], bool) & ~np.asarray(cat["FLAG_BCG"], bool)
    f, e = flux_matrix(cat, bands_of(cat))
    with np.errstate(invalid="ignore"):
        n_hi = (f / e >= HI_SNR).sum(1)
    return keep & (snr >= MIN_SNR) & (n_hi >= MIN_HI_BANDS) & behind


def tangent_xy(ra, dec, ra0: float, dec0: float) -> tuple[np.ndarray, np.ndarray]:
    """Flat-sky offsets in arcsec (x towards East, y North); fine over a cluster field."""
    ra, dec = np.asarray(ra, float), np.asarray(dec, float)
    return (ra - ra0) * np.cos(np.deg2rad(dec0)) * 3600.0, (dec - dec0) * 3600.0


def aperture_radius(cat: Table) -> np.ndarray:
    """Kron aperture semi-major axis (arcsec); the isophotal radius where KRON_RADIUS is 0."""
    kron = KRON_SCALE * np.asarray(cat["KRON_RADIUS"], float) * np.asarray(cat["A"], float)
    iso = np.sqrt(np.clip(np.asarray(cat["AREA_ISO"], float), 0, None) / np.pi)
    return np.where(np.isfinite(kron) & (kron > 0), kron, iso) * PIXSCALE


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
    x, y = tangent_xy(cat["RA"], cat["DEC"], ra0, dec0)
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
    deb = np.asarray(cat["FLAG_DEBLEND"], bool)

    def overlap(lo, hi):
        a_lo, a_hi = np.asarray(cat[lo], float), np.asarray(cat[hi], float)
        return (a_lo[ci] <= a_hi[cj]) & (a_lo[cj] <= a_hi[ci])

    boxes = overlap("X_MIN", "X_MAX") & overlap("Y_MIN", "Y_MAX")
    same = (sep < SAME_GALAXY_FACTOR * (rad[ci] + rad[cj])) | (deb[ci] & deb[cj] & boxes)
    # nearest other catalogued source (any redshift, BCGs included) to the midpoint, and to
    # the segment joining the members (a lens need not sit at the midpoint if the fluxes differ)
    good = np.isfinite(x) & np.isfinite(y)
    idx = np.flatnonzero(good)
    tree = cKDTree(np.column_stack([x[idx], y[idx]]))
    mx, my = 0.5 * (x[ci] + x[cj]), 0.5 * (y[ci] + y[cj])
    lens_dist = np.full(len(pairs), np.inf)
    line_dist = np.full(len(pairs), np.inf)
    lens_id = np.zeros(len(pairs), int)
    n_between = np.zeros(len(pairs), int)  # sources inside the circle on the pair's diameter
    source = np.asarray(cat["SOURCE"])
    near_mid = tree.query_ball_point(np.column_stack([mx, my]), 0.5 * sep + LENS_RADIUS)
    for p, cand in enumerate(near_mid):
        cand = idx[np.asarray(cand, int)]
        cand = cand[(cand != ci[p]) & (cand != cj[p])]
        if len(cand) == 0:
            continue
        dm = np.hypot(x[cand] - mx[p], y[cand] - my[p])
        k = int(np.argmin(dm))
        lens_dist[p], lens_id[p] = dm[k], source[cand[k]]
        # distance to the segment, for sources projecting inside it and clear of both members
        ux, uy = x[cj[p]] - x[ci[p]], y[cj[p]] - y[ci[p]]
        t = ((x[cand] - x[ci[p]]) * ux + (y[cand] - y[ci[p]]) * uy) / (ux * ux + uy * uy)
        dl = np.abs((x[cand] - x[ci[p]]) * uy - (y[cand] - y[ci[p]]) * ux) / np.hypot(ux, uy)
        clear = np.minimum(
            np.hypot(x[cand] - x[ci[p]], y[cand] - y[ci[p]]),
            np.hypot(x[cand] - x[cj[p]], y[cand] - y[cj[p]]),
        )
        inside = (t > 0) & (t < 1) & (clear >= LENS_RADIUS)
        if inside.any():
            line_dist[p] = dl[inside].min()
        n_between[p] = int(((dm <= 0.5 * sep[p]) & (clear >= LENS_RADIUS)).sum())
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
    X, Y, pa = (np.asarray(cat[c], float) for c in ("X", "Y", "PA"))
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
    """Steps 1-4 on one catalogue: ``(pairs, summary)``; ``pairs`` holds SED-matched close pairs."""
    sel = select_sources(cat, z_cluster)
    rows = np.flatnonzero(sel)
    sub = cat[rows]
    ra0, dec0 = float(np.nanmedian(cat["RA"])), float(np.nanmedian(cat["DEC"]))
    x, y = tangent_xy(sub["RA"], sub["DEC"], ra0, dec0)
    bands = bands_of(cat)
    f, e = flux_matrix(sub, bands)
    z16, z84 = np.asarray(sub["Z160"], float), np.asarray(sub["Z840"], float)

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
        cols = (("SOURCE", "id"), ("RA", "ra"), ("DEC", "dec"), ("Z_ML", "z"), ("Z_SPEC", "zspec"))
        cols += (("MU", "mu_canucs"), ("FLAG_POINTSRC", "pointsrc"))
        for col, out in cols:
            pairs[f"{out}_{tag}"] = np.asarray(cat[col])[pairs[k]]
    pairs["mirror_deg"] = mirror_angle(cat, pairs["ci"], pairs["cj"])
    pairs.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source="CANUCS DR1 photometry catalogue (colour apertures, EAzY photo-z)",
    )
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
    }
    return pairs, summary


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


def lens_check(top: Table, model_name: str, z_cluster: float) -> Table:
    """CATS map magnification (signed) at both members at their photo-z (model_prediction)."""
    import lens_consistency as lc

    model, _, _ = lc.load_model(model_name)
    lc.apply_frame_offset(model_name, model)
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
        np.abs(np.asarray(top["mu_canucs_a"], float)), np.abs(np.asarray(top["mu_canucs_b"], float))
    )
    big = np.where(np.isfinite(big), big, canucs)
    top["model_covered"] = known
    top["cluster_explains"] = top["parity_flip"] | (big > MAX_ORDINARY_MU)
    dl = dark_lens_numbers(top["sep"], z_cluster, zs)
    for c in dl.colnames:
        top[c] = dl[c]
    return top


def cutouts_and_sheet(top: Table, i2d_prefix: str, out_dir: Path, title: str) -> Path:
    """F150W/F277W/F444W cutouts centred on each pair's midpoint, members circled."""
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.wcs.utils import proj_plane_pixel_scales
    from matplotlib.figure import Figure
    from matplotlib.patches import Circle

    from jwst_anomaly import cutouts, viz

    man = Table.read(paths.repo_root() / "data/manifests" / FIELD_MANIFEST[i2d_prefix])
    names = [str(n) for n in man["productFilename"]]
    uris = {
        b: str(man["cloud_uri"][names.index(f"{i2d_prefix}_nircam_clear-{b}_i2d.fits")])
        for b in CUTOUT_BANDS
    }
    tgt = Table(
        {
            "source_uid": [f"r{r['rank']}_{r['id_a']}_{r['id_b']}" for r in top],
            "ra": 0.5 * (np.asarray(top["ra_a"]) + np.asarray(top["ra_b"])),
            "dec": 0.5 * (np.asarray(top["dec_a"]) + np.asarray(top["dec_b"])),
        }
    )
    cut = {b: cutouts.make_cutouts(uris[b], tgt, CUTOUT_ARCSEC, out_dir / "cutouts") for b in uris}
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
            path = cut[b]["path"][r]
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


FIELD_MANIFEST = {v["i2d"]: f"{k}_products.ecsv" for k, v in FIELDS.items()}


def load_images(field: str) -> Table:
    """All published multiple images for ``field``: CATS arcs.txt (JWST frame) + CANUCS lists."""
    import lens_consistency as lc

    spec = FIELDS[field]
    lists = []
    cats, _ = lc.image_list(spec["cats"], {}, None)
    lc.shift_images(spec["cats"], cats)
    lists.append(Table({"image_id": cats["image_id"], "ra": cats["ra"], "dec": cats["dec"]}))
    for rel in spec["image_lists"]:
        lists.append(read_image_list(canucs_dir() / rel))
    from astropy.table import vstack

    out = vstack(lists)
    out.meta["n_per_list"] = [len(t) for t in lists]
    return out


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
    cat = Table.read(canucs_dir() / _CAT.format(f=spec["canucs"]))
    images = load_images(args.field)
    pairs, summary = search(cat, spec["z_cluster"], images, args.n_shift, args.seed)
    summary["n_published_images"] = images.meta["n_per_list"]
    pairs.write(out / "matched_pairs.ecsv", overwrite=True)
    top = rank_orphans(pairs)[: args.top]
    if len(top):
        top = lens_check(top, spec["cats"], spec["z_cluster"])
        if args.cutouts:
            summary["contact_sheet"] = str(
                cutouts_and_sheet(top, spec["i2d"], out, f"{args.field} orphan pairs (unvetted)")
            )
        top.write(out / "top_orphans.ecsv", overwrite=True)
    summary["top"] = [
        {
            k: (row[k].item() if hasattr(row[k], "item") else row[k])
            for k in top.colnames
            if k not in ("i", "j", "ci", "cj")
        }
        for row in top
    ]
    (out / "summary.json").write_text(json.dumps(summary, indent=1, default=str))
    print(json.dumps({k: v for k, v in summary.items() if k != "top"}, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
