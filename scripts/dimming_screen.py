"""Multi-epoch dimming / vanished-source screen for NIRCam imaging (W3, D-052).

Signature W3 of docs/exotic_lensing.md (inverted microlensing by a negative mass) makes a compact
source vanish for 2 t_E sqrt(4 - u0²) between two caustic spikes. Ordinary sources fade too
(supernovae, AGN, microlensing by intracluster stars), and most catalogue-level "fades" are
artefacts (deblending, edges, diffraction spikes, persistence). This screen lists dimming events and
the cheapest ordinary flags; it never interprets them. Every threshold is an ASSUMPTION
(``Params``); outputs are ``derived``, injections ``simulated``.

Subcommands (configs/dimming_screen.yaml lists the fields and epochs):
- ``fetch``: one MAST query per program for the fields' level-3 ``_cat.ecsv`` (downloaded) and
  ``_i2d.fits`` (listed with S3 URIs, never downloaded); manifests under ``data/manifests/``.
- ``screen``: per-epoch catalogues -> master source list -> light curves (aper50 flux per band and
  epoch, frame and zero-point tied, errors scaled by the observed epoch-to-epoch scatter) -> flags:
  ``vanish`` (S/N >= 10 in one epoch, < 3 sigma in another), ``dim_achromatic`` (the same fractional
  drop in >= 2 bands within errors) and ``rise_dip_rise`` (>= 3 epochs, a significant dip bracketed
  by brighter epochs). Ordinary-explanation columns: Gaia star proximity (D-027 mask), point-like,
  blended, edge proxy, single-epoch detection (persistence suspect, D-039).
- ``inject``: injection-recovery on the monitored compact sources with
  ``exotic_sim.inject_light_curve`` (negative point mass, n = 1, sign = -1) over the real cadence,
  plus plain achromatic dimming of 20/50/100 %; writes the efficiency table and the rate limits.
- ``forced``: re-measures every flagged source (and random controls, for the noise scale) by forced
  aperture photometry on S3 byte-range cutouts of every epoch (``transient_forced.measure``),
  re-applies the same flag logic to the forced light curves, adds cutout quality (edge / no data,
  low weight, hexagonal spike statistic, D-018/D-021) and a SIMBAD/NED match, and draws an epoch
  contact sheet.

    python scripts/dimming_screen.py fetch --field nexus
    python scripts/dimming_screen.py screen --field nexus --out outputs/dimming/nexus
    python scripts/dimming_screen.py inject --field nexus --out outputs/dimming/nexus
    python scripts/dimming_screen.py forced --field nexus --out outputs/dimming/nexus
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
from epoch_compare import mutual_matches  # noqa: E402
from transient_combine import near_bright  # noqa: E402
from transient_forced import robust_std  # noqa: E402
from transient_search import _neighbour_counts  # noqa: E402

from jwst_anomaly import exotic_sim, paths, pipeline, schema  # noqa: E402

DAYS_PER_YEAR = 365.25
DEFAULT_CONFIG = Path(__file__).resolve().parents[1] / "configs" / "dimming_screen.yaml"
FLAGS = ("vanish", "dim_achromatic", "rise_dip_rise")


@dataclass(frozen=True)
class Params:
    """Screen thresholds; every value is an ASSUMPTION."""

    match_radius_arcsec: float = 0.3  # same as D-027
    ref_snr: float = 30.0  # frame / zero-point reference pairs (D-017)
    master_snr: float = 5.0  # a catalogue row enters the master list at this aper50 S/N
    detect_snr: float = 10.0  # "detected" in the vanish test
    gone_sigma: float = 3.0  # "below the limit" in the vanish test
    min_drop: float = 0.2  # fractional drop for dim / dip flags
    min_sigma: float = 5.0  # significance of a drop (scaled errors)
    achromatic_sigma: float = 3.0  # band drops agree within this many sigma
    footprint_radius_arcsec: float = 5.0  # coverage proxy (D-027)
    footprint_min: int = 3
    edge_fraction: float = 0.5  # neighbour count below this x the epoch median -> edge proxy
    blend_arcsec: float = 0.5  # nearest neighbour closer than this -> blended
    # Reference-epoch F200W CI_70_30 of the stellar locus is 2.0-2.5 (NEXUS o001: 1,280 of the S/N
    # >= 20 rows); below 1.9 a source is sharper than the PSF (hot pixel, cosmic ray: an artefact).
    point_ci_min: float = 1.9
    point_ci_max: float = 2.7
    gaia_self_arcsec: float = 0.3  # a Gaia DR3 source this close is the source itself (D-012)
    gaia_saturated_g: float = 18.0  # brighter Gaia stars saturate NIRCam and grow spikes
    sys_floor: float = 0.03  # fractional flux error floor added in quadrature (D-027: 0.03 mag)


DEFAULT_PARAMS = Params()


# --------------------------------------------------------------------------------------------
# configuration and data access


def load_field(config: Path, field: str) -> dict[str, Any]:
    import yaml

    cfg = yaml.safe_load(Path(config).read_text(encoding="utf-8"))
    if field not in cfg.get("fields", {}):
        raise SystemExit(f"error: field {field!r} not in {config}")
    f = dict(cfg["fields"][field])
    f["name"] = field
    f["bands"] = [b.upper() for b in cfg["bands"]]
    f["cloud"] = cfg["cloud"]
    f["epochs"] = sorted(f["epochs"], key=lambda e: float(e["mjd"]))
    return f


def obs_id(epoch: dict[str, Any], band: str) -> str:
    return f"{epoch['obs']}_nircam_clear-{band.lower()}"


def catalog_path(epoch: dict[str, Any], band: str) -> Path:
    o = obs_id(epoch, band)
    return paths.data_root() / "cache" / "mast" / o / f"{o}_cat.ecsv"


def cmd_fetch(args) -> dict:
    from jwst_anomaly import acquire, query

    field = load_field(args.config, args.field)
    by_program: dict[str, list[str]] = {}
    for e in field["epochs"]:
        for b in e["bands"]:
            by_program.setdefault(e["obs"][2:7], []).append(obs_id(e, b))
    prods = []
    for prog, ids in sorted(by_program.items()):  # one MAST query per program
        obs = query.query_observations(proposal_id=str(int(prog)), obs_id=sorted(ids))
        prods.append(query.list_products(obs, subgroups=("CAT", "I2D"), cloud_uris=True))
    from astropy.table import vstack

    allp = vstack(prods, metadata_conflicts="silent")
    mdir = paths.manifests_dir()
    acquire.write_ecsv(
        allp[[c for c in allp.colnames if c in _PRODUCT_KEEP]],
        mdir / f"dimming_{args.field}_products.ecsv",
    )
    cats = allp[np.char.endswith(np.asarray(allp["productFilename"]).astype(str), "_cat.ecsv")]
    rows = acquire.fetch_products(cats, manifest_path=mdir / f"dimming_{args.field}.ecsv")
    return {
        "field": args.field,
        "n_products": len(allp),
        "n_catalogues": len(rows),
        "catalogue_bytes": int(np.sum(rows["size"])),
        "missing": sorted(
            {obs_id(e, b) for e in field["epochs"] for b in e["bands"]}
            - {str(f).rsplit("_cat", 1)[0] for f in rows["productFilename"]}
        ),
    }


_PRODUCT_KEEP = (
    "obs_id",
    "dataURI",
    "productFilename",
    "size",
    "prvversion",
    "productSubGroupDescription",
    "cloud_uri",
)


def load_catalogs(field: dict[str, Any]) -> dict[tuple[int, str], Table]:
    from jwst_anomaly import catalog

    cats = {}
    for k, e in enumerate(field["epochs"]):
        for b in e["bands"]:
            cats[(k, b)] = catalog.load_pipeline_catalog(catalog_path(e, b))
    return cats


# --------------------------------------------------------------------------------------------
# light curves


def _col(t: Table, name: str) -> np.ndarray:
    return np.asarray(np.ma.filled(np.ma.asarray(t[name], float), np.nan), float)


def _flux(t: Table) -> tuple[np.ndarray, np.ndarray]:
    """aper50 flux and error in µJy (the catalogues give Jy)."""
    return _col(t, "aper50_flux") * 1e6, _col(t, "aper50_flux_err") * 1e6


def frame_shift(ref: SkyCoord, c: SkyCoord, snr_ref, snr, p: Params) -> np.ndarray:
    """Median (dRA cos δ, dDec) in arcsec of ``c`` relative to ``ref`` from bright mutual matches
    (the D-027 global tie); zero with fewer than 10 references."""
    i, j = mutual_matches(ref, c, p.match_radius_arcsec)
    bright = (snr_ref[i] >= p.ref_snr) & (snr[j] >= p.ref_snr)
    if bright.sum() < 10:
        return np.zeros(2)
    dra, ddec = ref[i[bright]].spherical_offsets_to(c[j[bright]])
    return np.array([np.median(dra.to_value(u.arcsec)), np.median(ddec.to_value(u.arcsec))])


def build_light_curves(
    cats: dict[tuple[int, str], Table],
    n_epochs: int,
    bands: list[str],
    p: Params = DEFAULT_PARAMS,
) -> dict[str, Any]:
    """Master source list and per-band, per-epoch fluxes.

    The deepest catalogue of the first (detection) band defines the frame. Every catalogue is
    shifted onto it (``frame_shift``). Master sources: the reference catalogue's rows with S/N >=
    ``master_snr``, then each other detection-band epoch's rows of that S/N without a master source
    within the match radius. Each (epoch, band) catalogue is matched to the master list by mutual
    nearest neighbours. Zero points: per band, every epoch is scaled to the band's reference epoch
    by the median flux ratio of S/N >= ``ref_snr`` pairs (unchanged with < 10). Coverage proxy: >=
    ``footprint_min`` catalogue rows within ``footprint_radius_arcsec``. A covered master source
    without a match gets flux 0 with the catalogue's typical faint-source error (a non-detection;
    ``forced`` measures it). Uncovered: NaN.
    """
    det_band = bands[0]
    keys_det = [k for k in range(n_epochs) if (k, det_band) in cats]
    ref_k = max(keys_det, key=lambda k: len(cats[(k, det_band)]))
    coords, snrs, shifts = {}, {}, {}
    ref_c = SkyCoord(cats[(ref_k, det_band)]["ra"], cats[(ref_k, det_band)]["dec"], unit="deg")
    f, e = _flux(cats[(ref_k, det_band)])
    ref_snr = f / e
    for key, t in cats.items():
        c = SkyCoord(t["ra"], t["dec"], unit="deg")
        fl, er = _flux(t)
        s = fl / er
        sh = frame_shift(ref_c, c, ref_snr, s, p) if key != (ref_k, det_band) else np.zeros(2)
        coords[key] = c.spherical_offsets_by(-sh[0] * u.arcsec, -sh[1] * u.arcsec)
        snrs[key], shifts[key] = s, sh
    # master list
    with np.errstate(invalid="ignore"):
        keep = snrs[(ref_k, det_band)] >= p.master_snr
    mra = list(coords[(ref_k, det_band)].ra.deg[keep])
    mdec = list(coords[(ref_k, det_band)].dec.deg[keep])
    for k in keys_det:
        if k == ref_k:
            continue
        c = coords[(k, det_band)]
        with np.errstate(invalid="ignore"):
            sel = np.flatnonzero(snrs[(k, det_band)] >= p.master_snr)
        if sel.size == 0:
            continue
        m = SkyCoord(mra, mdec, unit="deg")
        _, sep, _ = c[sel].match_to_catalog_sky(m)
        new = sel[sep.arcsec > p.match_radius_arcsec]
        mra += list(c.ra.deg[new])
        mdec += list(c.dec.deg[new])
    master = SkyCoord(mra, mdec, unit="deg")
    n, nb = len(master), len(bands)
    flux = np.full((n, n_epochs, nb), np.nan)
    err = np.full((n, n_epochs, nb), np.nan)
    det = np.zeros((n, n_epochs, nb), bool)
    nbr = np.zeros((n, n_epochs, nb))
    row = np.full((n, n_epochs, nb), -1, int)
    zp, depth = {}, {}
    for (k, b), t in cats.items():
        jb = bands.index(b)
        fl, er = _flux(t)
        i, j = mutual_matches(master, coords[(k, b)], p.match_radius_arcsec)
        _, counts = _neighbour_counts(master, coords[(k, b)], p.footprint_radius_arcsec)
        nbr[:, k, jb] = counts
        covered = counts >= p.footprint_min
        s = fl / er
        with np.errstate(invalid="ignore"):
            faint = np.isfinite(er) & (s < 20)
        d = float(np.median(er[faint])) if faint.any() else float(np.nanmedian(er))
        depth[f"{b}_e{k}"] = d
        flux[covered, k, jb] = 0.0
        err[covered, k, jb] = d
        flux[i, k, jb], err[i, k, jb] = fl[j], er[j]
        det[i, k, jb] = np.isfinite(fl[j])
        row[i, k, jb] = j
    for jb, b in enumerate(bands):
        ks = [k for k in range(n_epochs) if (k, b) in cats]
        rk = max(ks, key=lambda k: len(cats[(k, b)]))
        for k in ks:
            both = det[:, k, jb] & det[:, rk, jb]
            with np.errstate(invalid="ignore", divide="ignore"):
                good = both & (flux[:, k, jb] / err[:, k, jb] >= p.ref_snr)
                good &= flux[:, rk, jb] / err[:, rk, jb] >= p.ref_snr
                ratio = flux[good, rk, jb] / flux[good, k, jb]
            factor = float(np.median(ratio)) if good.sum() >= 10 else 1.0
            zp[f"{b}_e{k}"] = factor
            flux[:, k, jb] *= factor
            err[:, k, jb] *= factor
    return {
        "ra": master.ra.deg,
        "dec": master.dec.deg,
        "flux": flux,
        "err_raw": err,
        "det": det,
        "neighbours": nbr,
        "row": row,
        "ref_epoch": ref_k,
        "zero_point_factor": zp,
        "depth_1sigma_ujy": depth,
        "frame_shift_arcsec": {f"{b}_e{k}": [float(x) for x in s] for (k, b), s in shifts.items()},
    }


def noise_scale(flux: np.ndarray, err: np.ndarray, det: np.ndarray, min_n: int = 100) -> np.ndarray:
    """Per band, the robust std of pairwise epoch differences in units of their quoted error.

    For every epoch pair, ``(f_i - f_j) / sqrt(e_i² + e_j²)`` over sources detected in both; the
    median over pairs of the 1.4826 MAD is the band's scale (never below 1; 1 with fewer than
    ``min_n`` sources). Catalogue and ERR-based errors underestimate the scatter by 1.2-1.5x
    (D-027).
    """
    _, ne, nb = flux.shape
    out = np.ones(nb)
    for jb in range(nb):
        vals = []
        for i in range(ne):
            for j in range(i + 1, ne):
                both = det[:, i, jb] & det[:, j, jb]
                z = (flux[both, i, jb] - flux[both, j, jb]) / np.hypot(
                    err[both, i, jb], err[both, j, jb]
                )
                s = robust_std(z, min_n)
                if np.isfinite(s):
                    vals.append(s)
        if vals:
            out[jb] = max(1.0, float(np.median(vals)))
    return out


def scaled_errors(flux, err_raw, scale, sys_floor: float) -> np.ndarray:
    """Quoted errors times the band noise scale, plus a fractional floor in quadrature."""
    return np.hypot(err_raw * scale[None, None, :], sys_floor * np.abs(flux))


# --------------------------------------------------------------------------------------------
# flags


def mask_shallow_nondetections(
    flux: np.ndarray, err: np.ndarray, det: np.ndarray, p: Params = DEFAULT_PARAMS
) -> np.ndarray:
    """Copy of ``flux`` with uninformative catalogue non-detections set to NaN (uncovered).

    A covered epoch without a catalogue match counts only where the source's brightest detected flux
    in that band would have been >= ``detect_snr`` times that epoch's error: catalogues are assumed
    complete at S/N >= 10 (ASSUMPTION), and below that a missing row says nothing about the flux.
    """
    with _quiet():
        peak = np.nanmax(np.where(det, flux, np.nan), axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        informative = peak[:, None, :] / err >= p.detect_snr
    out = flux.copy()
    out[~det & ~informative] = np.nan
    return out


def classify(
    flux: np.ndarray, err: np.ndarray, p: Params = DEFAULT_PARAMS
) -> dict[str, np.ndarray]:
    """Dimming flags from light curves ``flux``/``err`` of shape (n, epochs, bands), epochs in time
    order; NaN marks an uncovered epoch. Returns per-source arrays:

    - ``vanish``: an epoch pair (i, j) where every band with S/N_i >= ``detect_snr`` and coverage in
      j has S/N_j < ``gone_sigma`` and a drop significant at ``min_sigma`` (at least one band);
      ``vanish_bands`` counts those bands.
    - ``dim_bands`` (max over epochs of the bands where that epoch is >= ``min_drop`` fainter than
      the median of the other epochs at >= ``min_sigma``), ``dim_achromatic`` when >= 2 such bands
      in one epoch have drops consistent within ``achromatic_sigma``, ``dim_epoch`` (that epoch,
      else -1), ``max_drop`` (largest fractional drop in a flagged band, else of any band).
    - ``rdr_bands``: bands where an epoch is fainter than an earlier and a later one, each by
      >= ``min_sigma`` and >= ``min_drop`` of the fainter of the two; ``rise_dip_rise`` when any
      band qualifies.
    """
    n, ne, nb = flux.shape
    cov = np.isfinite(flux) & np.isfinite(err) & (err > 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        snr = np.where(cov, flux / err, np.nan)
    # vanish: some epoch pair (i, j) with S/N_i >= detect_snr, S/N_j < gone_sigma and a drop
    # f_i - f_j significant at min_sigma (so a shallow epoch's non-detection does not count)
    # (achromatic: every band detected at >= detect_snr in epoch i and covered in j must qualify)
    vanish = np.zeros(n, bool)
    vanish_nb = np.zeros(n, int)
    for i in range(ne):
        for j in range(ne):
            if i == j:
                continue
            with np.errstate(invalid="ignore", divide="ignore"):
                drop_sig = (flux[:, i] - flux[:, j]) / np.hypot(err[:, i], err[:, j])
                testable = (snr[:, i] >= p.detect_snr) & cov[:, j]
                q = testable & (snr[:, j] < p.gone_sigma) & (drop_sig >= p.min_sigma)
            ok = (q | ~testable).all(axis=1) & q.any(axis=1)
            vanish |= ok
            vanish_nb = np.where(ok, np.maximum(vanish_nb, q.sum(axis=1)), vanish_nb)
    # dim relative to the median of the other epochs
    drop = np.full((n, ne, nb), np.nan)
    sig = np.full((n, ne, nb), np.nan)
    sdrop = np.full((n, ne, nb), np.nan)
    for k in range(ne):
        others = np.delete(np.arange(ne), k)
        fo = np.where(cov[:, others], flux[:, others], np.nan)
        eo = np.where(cov[:, others], err[:, others], np.nan)
        no = cov[:, others].sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"), _quiet():
            base = np.nanmedian(fo, axis=1)
            ebase = 1.2533 * np.sqrt(np.nanmean(eo**2, axis=1) / np.maximum(no, 1))
            ok = cov[:, k] & (no >= 1) & (base > 0)
            drop[:, k] = np.where(ok, 1.0 - flux[:, k] / base, np.nan)
            sig[:, k] = np.where(ok, (base - flux[:, k]) / np.hypot(err[:, k], ebase), np.nan)
            sdrop[:, k] = np.where(
                ok, np.hypot(err[:, k] / base, flux[:, k] * ebase / base**2), np.nan
            )
    with np.errstate(invalid="ignore"):
        dim = (drop >= p.min_drop) & (sig >= p.min_sigma)
    dim_bands = dim.sum(axis=2)
    ach = np.zeros((n, ne), bool)
    for a in range(nb):
        for b in range(a + 1, nb):
            with np.errstate(invalid="ignore"):
                agree = np.abs(drop[:, :, a] - drop[:, :, b]) <= p.achromatic_sigma * np.hypot(
                    sdrop[:, :, a], sdrop[:, :, b]
                )
            ach |= dim[:, :, a] & dim[:, :, b] & agree
    dim_ach = ach.any(axis=1)
    dim_epoch = np.where(dim_ach, np.argmax(ach, axis=1), -1)
    with _quiet():
        max_drop = np.nanmax(np.where(dim, drop, np.nan), axis=(1, 2))
        max_drop = np.where(np.isfinite(max_drop), max_drop, np.nanmax(drop, axis=(1, 2)))
    # rise - dip - rise
    rdr_b = np.zeros((n, nb), bool)
    for k in range(1, ne - 1):
        for i in range(k):
            for j in range(k + 1, ne):
                with np.errstate(invalid="ignore", divide="ignore"):
                    s1 = (flux[:, i] - flux[:, k]) / np.hypot(err[:, i], err[:, k])
                    s2 = (flux[:, j] - flux[:, k]) / np.hypot(err[:, j], err[:, k])
                    lower = np.fmin(flux[:, i], flux[:, j])
                    d = 1.0 - flux[:, k] / lower
                    hit = (s1 >= p.min_sigma) & (s2 >= p.min_sigma) & (d >= p.min_drop)
                rdr_b |= hit & cov[:, i] & cov[:, k] & cov[:, j]
    return {
        "vanish": vanish,
        "vanish_bands": vanish_nb,
        "dim_achromatic": dim_ach,
        "dim_bands": dim_bands.max(axis=1),
        "dim_epoch": dim_epoch,
        "max_drop": max_drop,
        "rise_dip_rise": rdr_b.any(axis=1),
        "rdr_bands": rdr_b.sum(axis=1),
        "max_snr": np.nanmax(np.where(cov, snr, -np.inf), axis=(1, 2)),
    }


class _quiet:
    """Silence numpy's all-NaN-slice RuntimeWarnings inside a block."""

    def __enter__(self):
        import warnings

        self._w = warnings.catch_warnings()
        self._w.__enter__()
        warnings.simplefilter("ignore", RuntimeWarning)

    def __exit__(self, *exc):
        self._w.__exit__(*exc)


def flagged(res: dict[str, np.ndarray]) -> np.ndarray:
    return res["vanish"] | res["dim_achromatic"] | res["rise_dip_rise"]


def ordinary_columns(
    lc: dict[str, Any],
    cats: dict[tuple[int, str], Table],
    bands: list[str],
    p: Params = DEFAULT_PARAMS,
    gaia: Table | None = None,
) -> dict[str, np.ndarray]:
    """Catalogue-level ordinary-explanation columns (each a flag; True = ordinary explanation open).

    - ``near_star``: within the D-027 Gaia exclusion radius (``transient_combine.near_bright``).
    - ``point_like``: ``point_ci_min`` <= CI_70_30 <= ``point_ci_max`` in the reference-epoch
      detection-band row (not an ordinary flag; it selects compact sources).
    - ``blended``: nearest catalogue neighbour within ``blend_arcsec`` there.
    - ``edge_proxy``: in some covered epoch-band the neighbour count is below ``edge_fraction`` x
      that catalogue's median (mosaic edge or low weight; cutouts check it in ``forced``).
    - ``gaia_star``: a Gaia DR3 source within ``gaia_self_arcsec`` (a Galactic star: variable stars
      are the ordinary explanation; not a veto).
    - ``sharp_artifact``: CI_70_30 below ``point_ci_min`` (sharper than the PSF).
    - ``single_epoch``: detected (catalogue match) in exactly one detection-band epoch (persistence
      or artefact suspect, D-039).
    """
    n = len(lc["ra"])
    det = lc["det"]
    rk = lc["ref_epoch"]
    t = cats[(rk, bands[0])]
    rows = lc["row"][:, rk, 0]
    has = rows >= 0
    ci = np.full(n, np.nan)
    nn = np.full(n, np.nan)
    if "CI_70_30" in t.colnames:
        ci[has] = _col(t, "CI_70_30")[rows[has]]
    if "nn_dist" in t.colnames:
        nn[has] = _col(t, "nn_dist")[rows[has]] * float(t.meta.get("pixel_scale_arcsec") or np.nan)
    with np.errstate(invalid="ignore"):
        point = has & (ci >= p.point_ci_min) & (ci <= p.point_ci_max)
        sharp = has & (ci < p.point_ci_min)
        blended = nn < p.blend_arcsec
    nbr = lc["neighbours"]
    edge = np.zeros(n, bool)
    for (k, b), _ in cats.items():
        jb = bands.index(b)
        c = nbr[:, k, jb]
        covered = c >= p.footprint_min
        if covered.any():
            edge |= covered & (c < p.edge_fraction * np.median(c[covered]))
    near = np.zeros(n, bool)
    star = np.zeros(n, bool)
    if gaia is not None and len(gaia):
        hit = near_bright(lc["ra"], lc["dec"], gaia)
        pos = SkyCoord(lc["ra"], lc["dec"], unit="deg")
        g = SkyCoord(gaia["ra"], gaia["dec"], unit="deg")
        _, sep, _ = pos.match_to_catalog_sky(g)
        star = sep.arcsec <= p.gaia_self_arcsec
        gmag = np.asarray(gaia["gmag"], float)
        # the source itself is the Gaia match: only a saturating star counts as "near a star"
        own = star & (hit >= 0)
        bright_self = own & (gmag[np.clip(hit, 0, None)] < p.gaia_saturated_g)
        near = (hit >= 0) & (~own | bright_self)
    return {
        "near_star": near,
        "gaia_star": star,
        "sharp_artifact": sharp,
        "point_like": point,
        "ci_70_30": ci,
        "blended": blended,
        "edge_proxy": edge,
        "single_epoch": det[:, :, 0].sum(axis=1) == 1,
        "n_epochs_detected": det[:, :, 0].sum(axis=1),
    }


def catalogue_veto(ordinary: dict[str, np.ndarray]) -> np.ndarray:
    """Rows with an open catalogue-level ordinary explanation (bright star, edge, blend, single
    epoch, sharper than the PSF)."""
    return (
        ordinary["near_star"]
        | ordinary["edge_proxy"]
        | ordinary["blended"]
        | ordinary["single_epoch"]
        | ordinary["sharp_artifact"]
    )


def monitored(
    lc: dict[str, Any], err: np.ndarray, ordinary: dict[str, np.ndarray], p: Params = DEFAULT_PARAMS
):
    """Compact sources the limit is computed for: point-like, no catalogue veto, detected at S/N >=
    ``detect_snr`` in the detection band in >= 2 covered epochs (so a later vanish is testable)."""
    f = lc["flux"][:, :, 0]
    with np.errstate(invalid="ignore", divide="ignore"):
        good = (lc["det"][:, :, 0]) & (f / err[:, :, 0] >= p.detect_snr)
    return ordinary["point_like"] & ~catalogue_veto(ordinary) & (good.sum(axis=1) >= 2)


# --------------------------------------------------------------------------------------------
# injection-recovery


def reference_flux(flux: np.ndarray, det: np.ndarray) -> np.ndarray:
    """Per source and band, the median catalogue flux over detected epochs (the unlensed flux)."""
    with _quiet():
        return np.nanmedian(np.where(det, flux, np.nan), axis=1)


# Peak magnification of the n = 1 negative-mass caustic spike for a uniform disk of radius rho
# (simulated; D-047 as merged in PR #67, checked there against inverse ray shooting). The pre-merge
# simulator gave 9.2 / 4.9 for rho = 0.01 / 0.1; D-047's first draft quoted x7.5 / x3.4.
SPIKE_PEAK = {0.01: 7.0, 0.1: 2.35, 0.3: 1.53}


def w3_factor(
    times_yr: np.ndarray,
    t0: np.ndarray,
    t_e: float,
    u0: np.ndarray,
    rho: float,
    blend: float = 1.0,
    spike_peak: float | None = None,
) -> np.ndarray:
    """Flux factor ``blend A(t) + 1 - blend`` per row and epoch for a negative point mass (n = 1,
    sign = -1), from ``exotic_sim.inject_light_curve`` (``simulated``).

    ``spike_peak`` optionally caps ``A`` (a sensitivity check on the spike height; None = the
    simulator's value, ``SPIKE_PEAK``). The umbra (A = 0) is the robust part of the signal: the
    spike only matters for epochs within about rho t_E of a caustic crossing.
    """
    cap = np.inf if spike_peak is None else spike_peak
    out = np.empty((len(t0), len(times_yr)))
    for s in range(len(t0)):
        lcv = exotic_sim.inject_light_curve(
            times_yr, 1.0, float(t0[s]), t_e, float(u0[s]), n=1.0, sign=-1, rho=rho
        )
        a = np.minimum(np.asarray(lcv["magnification"]), cap)
        out[s] = blend * a + 1.0 - blend
    return out


def apply_factor(
    flux: np.ndarray, err: np.ndarray, factor: np.ndarray, rng, sys_floor: float = 0.03
) -> tuple[np.ndarray, np.ndarray]:
    """Injected light curves ``F f_obs + (1 - F) n``, the same factor ``F`` (rows x epochs) in every
    band (achromatic), and their errors (``simulated``).

    ``n`` ~ N(0, sigma_noise), where sigma_noise is the noise part of the error (the fractional
    floor removed): the observed epoch-to-epoch scatter scales with the flux, and a vanished source
    (F = 0) leaves sky noise only. The errors keep the noise part and recompute the floor at the new
    flux.
    """
    with np.errstate(invalid="ignore"):
        noise = np.sqrt(np.clip(err**2 - (sys_floor * flux) ** 2, 0, None))
    f = factor[:, :, None]
    out = f * flux + (1.0 - f) * rng.normal(0.0, 1.0, flux.shape) * noise
    return out, np.sqrt(noise**2 + (sys_floor * out) ** 2)


def dimming_factor(n_epochs: int, epoch: np.ndarray, depth: float) -> np.ndarray:
    """Factor ``1 - depth`` in one epoch per row, 1 elsewhere (plain achromatic dimming)."""
    out = np.ones((len(epoch), n_epochs))
    out[np.arange(len(epoch)), epoch] = 1.0 - depth
    return out


MAG_BINS = (15.0, 22.0, 24.0, 25.0, 26.0, 27.0, 29.0)


def ab_mag(flux_ujy: np.ndarray) -> np.ndarray:
    with np.errstate(invalid="ignore", divide="ignore"):
        return 23.9 - 2.5 * np.log10(flux_ujy)


def efficiency_table(
    lc: dict[str, Any],
    err: np.ndarray,
    veto: np.ndarray,
    sel: np.ndarray,
    times_yr: np.ndarray,
    t_e_grid: tuple[float, ...],
    rho_grid: tuple[float, ...],
    depths: tuple[float, ...],
    p: Params = DEFAULT_PARAMS,
    seed: int = 0,
    n_rep: int = 20,
    forced_inflation: np.ndarray | None = None,
) -> Table:
    """Recovery fraction per magnification model and magnitude bin (``simulated``).

    Injections go into selected sources that are not flagged without an injection (so a recovery is
    the injection's), ``n_rep`` times each. W3: every copy gets one event per (t_E, rho) with u0 ~
    U[0, 2) (umbra crossings) and t0 ~ U[t_first - 2 t_E, t_last + 2 t_E], so the umbra (half-length
    <= 2 t_E) can overlap the baseline; recovered = any flag and no catalogue veto. ``window_yr`` =
    T + 4 t_E is the t0 range, so ``efficiency x window_yr`` is the exposure per source in years.
    Dimming: one random covered detection-band epoch per source is dimmed by ``depth`` in every
    band.
    """
    rng = np.random.default_rng(seed)
    base = classify(mask_shallow_nondetections(lc["flux"], err, lc["det"], p), err, p)
    idx = np.repeat(np.flatnonzero(sel & ~flagged(base)), n_rep)
    flux = lc["flux"][idx]
    det = lc["det"][idx]
    e = err[idx]
    v = veto[idx]
    mag = ab_mag(reference_flux(flux, det)[:, 0])
    t_first, t_last = float(times_yr.min()), float(times_yr.max())
    rows = []

    def record(kind, t_e, rho, depth, rec, window):
        for lo, hi in zip(MAG_BINS[:-1], MAG_BINS[1:], strict=True):
            m = (mag >= lo) & (mag < hi)
            rows.append(
                {
                    "model": kind,
                    "t_e_yr": t_e,
                    "rho": rho,
                    "depth": depth,
                    "mag_lo": lo,
                    "mag_hi": hi,
                    "n_sources": int(m.sum()) // n_rep,
                    "n_injected": int(m.sum()),
                    "n_recovered": int(rec[m].sum()),
                    "efficiency": float(rec[m].mean()) if m.any() else np.nan,
                    "window_yr": window,
                }
            )

    for t_e in t_e_grid:
        for rho in rho_grid:
            t0 = rng.uniform(t_first - 2 * t_e, t_last + 2 * t_e, len(flux))
            u0 = rng.uniform(0.0, 2.0, len(flux))
            fi, ei = apply_factor(flux, e, w3_factor(times_yr, t0, t_e, u0, rho), rng, p.sys_floor)
            rec = _recovered(fi, ei, det, p, forced_inflation) & ~v
            record("w3", t_e, rho, np.nan, rec, (t_last - t_first) + 4 * t_e)
    covered = np.isfinite(flux[:, :, 0])
    for depth in depths:
        ep = np.array([rng.choice(np.flatnonzero(c)) for c in covered])
        fac = dimming_factor(flux.shape[1], ep, depth)
        fi, ei = apply_factor(flux, e, fac, rng, p.sys_floor)
        rec = _recovered(fi, ei, det, p, forced_inflation) & ~v
        record("dimming", np.nan, np.nan, depth, rec, np.nan)
    out = Table(rows=rows)
    out.meta.update(
        provenance=schema.Provenance.SIMULATED.value,
        source="dimming_screen.efficiency_table: exotic_sim.inject_light_curve(n=1, sign=-1) "
        "on real light curves",
        mag_band="detection band (reference flux)",
        seed=seed,
        n_rep=n_rep,
        params=asdict(p),
    )
    return out


def _recovered(fi, ei, det, p: Params, inflation: np.ndarray | None) -> np.ndarray:
    """Flagged at the catalogue stage and, with ``inflation`` (per band, >= 1), still flagged when
    the errors are inflated to the forced-photometry noise (a proxy for the forced confirmation,
    which only needs some flag to reappear)."""
    rec = flagged(_classify_injected(fi, ei, det, p))
    if inflation is not None:
        rec &= flagged(_classify_injected(fi, ei * np.asarray(inflation)[None, None, :], det, p))
    return rec


def _classify_injected(fi, ei, det, p: Params) -> dict[str, np.ndarray]:
    """Classify injected light curves as the screen does: rows whose injected flux falls below
    ``master_snr`` count as catalogue non-detections (kept at their noisy injected flux)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        det_i = det & (fi / ei >= p.master_snr)
    return classify(mask_shallow_nondetections(fi, ei, det_i, p), ei, p)


def rate_limits(
    eff: Table, n_by_bin: dict[tuple[float, float], int], n_pairs: int, area_deg2: float
) -> Table:
    """95 % upper limits for zero surviving events, per W3 model (t_E, rho).

    Poisson: 3.0 events. Per source per year: 3 / Σ_bins N_bin ε_bin window. Per source per epoch
    pair (the owner's form): 3 / (Σ_bins N_bin ε_bin x n_pairs). Per deg² per epoch pair:
    3 / (area x n_pairs x ε̄), ε̄ the source-weighted mean efficiency. All ``derived``.
    """
    rows = []
    w3 = eff[eff["model"] == "w3"]
    for t_e in np.unique(w3["t_e_yr"]):
        for rho in np.unique(w3["rho"]):
            s = w3[(w3["t_e_yr"] == t_e) & (w3["rho"] == rho)]
            ne = sum(
                n_by_bin.get((float(r["mag_lo"]), float(r["mag_hi"])), 0)
                * (r["efficiency"] if np.isfinite(r["efficiency"]) else 0.0)
                for r in s
            )
            ntot = sum(n_by_bin.values())
            window = float(s["window_yr"][0])
            rows.append(
                {
                    "t_e_yr": float(t_e),
                    "rho": float(rho),
                    "n_monitored": ntot,
                    "n_effective": float(ne),
                    "mean_efficiency": float(ne / ntot) if ntot else np.nan,
                    "window_yr": window,
                    "exposure_source_yr": float(ne * window),
                    "limit_per_source_per_yr": 3.0 / (ne * window) if ne > 0 else np.inf,
                    "limit_per_source_per_pair": 3.0 / (ne * n_pairs) if ne > 0 else np.inf,
                    "limit_per_deg2_per_pair": (
                        3.0 / (area_deg2 * n_pairs * ne / ntot) if ne > 0 else np.inf
                    ),
                }
            )
    return Table(rows=rows, meta={"provenance": "derived", "poisson_95": 3.0})


def area_deg2(ra, dec, cov_epochs: np.ndarray, cell_arcsec: float = 5.0) -> float:
    """Sky area of cells holding a source covered in >= 2 epochs (proxy; ``cell_arcsec`` grid)."""
    if len(ra) == 0:
        return 0.0
    keep = cov_epochs >= 2
    d0 = float(np.median(dec))
    x = np.asarray(ra)[keep] * np.cos(np.deg2rad(d0)) * 3600 / cell_arcsec
    y = np.asarray(dec)[keep] * 3600 / cell_arcsec
    cells = {(int(np.floor(a)), int(np.floor(b))) for a, b in zip(x, y, strict=True)}
    return len(cells) * (cell_arcsec / 3600) ** 2


# --------------------------------------------------------------------------------------------
# commands


def _times_yr(field: dict[str, Any]) -> np.ndarray:
    return np.array([float(e["mjd"]) for e in field["epochs"]]) / DAYS_PER_YEAR


def _screen(field: dict[str, Any], p: Params, gaia: Table | None):
    cats = load_catalogs(field)
    bands = field["bands"]
    lc = build_light_curves(cats, len(field["epochs"]), bands, p)
    scale = noise_scale(lc["flux"], lc["err_raw"], lc["det"])
    err = scaled_errors(lc["flux"], lc["err_raw"], scale, p.sys_floor)
    res = classify(mask_shallow_nondetections(lc["flux"], err, lc["det"], p), err, p)
    ordn = ordinary_columns(lc, cats, bands, p, gaia)
    return cats, lc, scale, err, res, ordn


def _gaia(field: dict[str, Any], out: Path) -> Table:
    from transient_combine import fetch_gaia

    path = out / "gaia.ecsv"
    if path.exists():
        return Table.read(path)
    g = fetch_gaia(*field["centre"], radius_arcmin=6.0)
    out.mkdir(parents=True, exist_ok=True)
    g.write(path, overwrite=True)
    return g


def lightcurve_table(field, lc, err, res, ordn) -> Table:
    bands = field["bands"]
    t = Table({"uid": [f"{field['name']}-d{k:05d}" for k in range(len(lc["ra"]))]})
    t["ra"], t["dec"] = lc["ra"], lc["dec"]
    for jb, b in enumerate(bands):
        for k in range(len(field["epochs"])):
            t[f"{b}_e{k}_flux"] = lc["flux"][:, k, jb]
            t[f"{b}_e{k}_err"] = err[:, k, jb]
            t[f"{b}_e{k}_det"] = lc["det"][:, k, jb]
    for name, v in {**res, **ordn}.items():
        t[name] = v
    t["flagged"] = flagged(res)
    t["catalogue_veto"] = catalogue_veto(ordn)
    return t


def cmd_screen(args) -> dict:
    p = Params()
    field = load_field(args.config, args.field)
    t_start = time.time()
    gaia = None if args.no_gaia else _gaia(field, args.out)
    cats, lc, scale, err, res, ordn = _screen(field, p, gaia)
    t = lightcurve_table(field, lc, err, res, ordn)
    times = _times_yr(field)
    t.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source=f"scripts/dimming_screen.py screen --field {args.field}",
        epochs=[{"obs": e["obs"], "mjd": float(e["mjd"])} for e in field["epochs"]],
        bands=field["bands"],
        noise_scale={b: float(s) for b, s in zip(field["bands"], scale, strict=True)},
        zero_point_factor=lc["zero_point_factor"],
        frame_shift_arcsec=lc["frame_shift_arcsec"],
        depth_1sigma_ujy=lc["depth_1sigma_ujy"],
        catalogues={f"{b}_e{k}": c.meta.get("input_sha256") for (k, b), c in cats.items()},
        params=asdict(p),
        thresholds_provenance="assumption",
    )
    args.out.mkdir(parents=True, exist_ok=True)
    t.write(args.out / "lightcurves.ecsv", overwrite=True)
    fl = t[t["flagged"]]
    fl.write(args.out / "flags.ecsv", overwrite=True)
    mon = monitored(lc, err, ordn, p)
    summary = {
        "field": args.field,
        "n_sources": len(t),
        "n_epochs": len(field["epochs"]),
        "baseline_yr": float(times.max() - times.min()),
        "noise_scale": t.meta["noise_scale"],
        "n_point_like": int(ordn["point_like"].sum()),
        "n_monitored_compact": int(mon.sum()),
        "flags": {k: int(res[k].sum()) for k in FLAGS},
        "n_flagged": len(fl),
        "n_flagged_no_catalogue_veto": int((t["flagged"] & ~t["catalogue_veto"]).sum()),
        "n_flagged_compact_no_veto": int((t["flagged"] & ~t["catalogue_veto"] & mon).sum()),
        "veto_counts_among_flagged": {
            k: int(np.sum(fl[k]))
            for k in (
                "near_star",
                "gaia_star",
                "sharp_artifact",
                "edge_proxy",
                "blended",
                "single_epoch",
                "point_like",
            )
        },
        "wall_s": round(time.time() - t_start, 1),
    }
    (args.out / "screen_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


T_E_GRID = (0.01, 0.03, 0.1, 0.3, 1.0, 3.0)  # yr; ASSUMPTION grid
RHO_GRID = (0.01, 0.1)  # source radius / Einstein radius (D-047 recommended pair)
DEPTHS = (0.2, 0.5, 1.0)


def cmd_inject(args) -> dict:
    p = Params()
    field = load_field(args.config, args.field)
    t_start = time.time()
    gaia = None if args.no_gaia else _gaia(field, args.out)
    cats, lc, scale, err, res, ordn = _screen(field, p, gaia)
    sel = monitored(lc, err, ordn, p)
    times = _times_yr(field)
    veto = catalogue_veto(ordn)
    inflation = None
    fs = args.out / "forced.ecsv"
    if fs.exists():  # forced noise scale over the catalogue-stage one, per band (never below 1)
        forced_scale = Table.read(fs).meta["noise_scale"]
        inflation = np.array(
            [max(1.0, forced_scale[b] / s) for b, s in zip(field["bands"], scale, strict=True)]
        )
    eff = efficiency_table(
        lc,
        err,
        veto,
        sel,
        times,
        T_E_GRID,
        RHO_GRID,
        DEPTHS,
        p,
        seed=args.seed,
        forced_inflation=inflation,
    )
    eff.meta["forced_inflation"] = None if inflation is None else [float(x) for x in inflation]
    mag = ab_mag(reference_flux(lc["flux"][sel], lc["det"][sel])[:, 0])
    n_by_bin = {
        (lo, hi): int(((mag >= lo) & (mag < hi)).sum())
        for lo, hi in zip(MAG_BINS[:-1], MAG_BINS[1:], strict=True)
    }
    covd = np.isfinite(lc["flux"][:, :, 0]).sum(axis=1)
    area = area_deg2(lc["ra"], lc["dec"], covd)
    ne = len(field["epochs"])
    n_pairs = ne * (ne - 1) // 2
    lim = rate_limits(eff, n_by_bin, n_pairs, area)
    lim.meta.update(field=args.field, area_deg2=area, n_pairs=n_pairs)
    args.out.mkdir(parents=True, exist_ok=True)
    eff.write(args.out / "efficiency.ecsv", overwrite=True)
    lim.write(args.out / "limits.ecsv", overwrite=True)
    summary = {
        "field": args.field,
        "n_monitored": int(sel.sum()),
        "n_by_mag_bin": {f"{lo}-{hi}": n for (lo, hi), n in n_by_bin.items()},
        "area_deg2": area,
        "n_epoch_pairs": n_pairs,
        "baseline_yr": float(times.max() - times.min()),
        "forced_inflation": eff.meta["forced_inflation"],
        "wall_s": round(time.time() - t_start, 1),
    }
    (args.out / "inject_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


# --------------------------------------------------------------------------------------------
# forced photometry of flags (network: S3 byte-range cutouts)


def cmd_forced(args) -> dict:
    from transient_forced import measure, select_controls

    from jwst_anomaly import crossmatch, cutouts

    p = Params()
    field = load_field(args.config, args.field)
    t_start = time.time()
    lct = Table.read(args.out / "lightcurves.ecsv")
    fl = lct[lct["flagged"] & ~lct["catalogue_veto"]]  # cheapest ordinary tests first
    if args.max_flags and len(fl) > args.max_flags:
        # every point-like flag (the W3 population), plus a random audit sample of the others
        compact = np.asarray(fl["point_like"], bool)
        rest = np.flatnonzero(~compact)
        n_audit = max(0, args.max_flags - int(compact.sum()))
        audit = np.random.default_rng(0).choice(rest, min(n_audit, rest.size), replace=False)
        fl = fl[np.sort(np.concatenate([np.flatnonzero(compact), audit]))]
    # controls: unflagged master sources covered in >= 2 detection-band epochs (so epoch pairs
    # have enough of them for the noise scale), drawn by transient_forced.select_controls
    det_band = field["bands"][0]
    ne_all = len(field["epochs"])
    fcols = [
        f"{det_band}_e{k}_flux" for k in range(ne_all) if f"{det_band}_e{k}_flux" in lct.colnames
    ]
    fl_det = np.column_stack([np.asarray(lct[c], float) for c in fcols])
    det_cols = [c.replace("_flux", "_det") for c in fcols]
    dmask = np.column_stack([np.asarray(lct[c], bool) for c in det_cols])
    with _quiet():
        ref_mag = ab_mag(np.nanmedian(np.where(dmask, fl_det, np.nan), axis=1))
    pool = lct[~np.asarray(lct["flagged"], bool) & (np.isfinite(fl_det).sum(axis=1) >= 2)]
    pool_tab = Table(
        {
            "ra": pool["ra"],
            "dec": pool["dec"],
            "aper_total_abmag": ref_mag[
                ~np.asarray(lct["flagged"], bool) & (np.isfinite(fl_det).sum(axis=1) >= 2)
            ],
        }
    )
    ctl = select_controls(
        pool_tab,
        args.n_controls,
        tuple(args.control_mag),
        np.asarray(fl["ra"]),
        np.asarray(fl["dec"]),
        seed=0,
    )
    uids = list(fl["uid"]) + [f"ctl{k:04d}" for k in range(len(ctl))]
    targets = Table(
        {
            "source_uid": uids,
            "ra": np.concatenate([np.asarray(fl["ra"]), np.asarray(ctl["ra"])]),
            "dec": np.concatenate([np.asarray(fl["dec"]), np.asarray(ctl["dec"])]),
        }
    )
    n, ne, bands = len(targets), len(field["epochs"]), field["bands"]
    is_ctl = np.arange(n) >= len(fl)
    flux = np.full((n, ne, len(bands)), np.nan)
    ferr = np.full((n, ne, len(bands)), np.nan)
    qual = {k: np.zeros(n, bool) for k in ("cut_edge", "cut_low_weight", "cut_spike")}
    files: dict[tuple[int, str], list[str]] = {}
    skipped: list[str] = []
    for jb, b in enumerate(bands):
        scales = {}
        for k, e in enumerate(field["epochs"]):
            if b not in e["bands"]:
                continue
            uri = pipeline.l3_image_uri(field["cloud"], obs_id(e, b))
            ct = None
            for _attempt in range(2):  # S3 reads through a proxy fail now and then
                try:
                    ct = cutouts.make_cutouts(
                        uri,
                        targets,
                        size_arcsec=2.0,
                        extensions=["ERR", "WHT"],
                        spike_radii_arcsec=(0.2, 0.8),
                        out_dir=args.out / "cutouts" / f"{b}_e{k}",
                    )
                    break
                except OSError as exc:
                    print(f"warning: {b} e{k}: {exc}", file=sys.stderr)
            if ct is None:  # the epoch stays unmeasured (NaN), recorded in meta
                skipped.append(f"{b}_e{k}")
                continue
            scales[k] = float(np.nanmean(ct.meta["pixel_scale_arcsec"]))
            m = {str(a): r for a, r in zip(ct["source_uid"], ct, strict=True)}
            files[(k, b)] = [str(m[a]["path"]) if a in m else "" for a in uids]
            inside = np.array([a in m and str(m[a]["quality_flag"]) != "outside" for a in uids])
            for key, test in (
                ("cut_edge", lambda r: bool(r["on_edge"]) or float(r["frac_nan"]) > 0.1),
                ("cut_low_weight", lambda r: float(r["wht_rel"]) < 0.5),
                ("cut_spike", lambda r: float(r["spike_s6"]) >= 3.0),
            ):
                qual[key] |= inside & np.array(
                    [
                        bool(test(m[a])) if a in m and inside[i] else False
                        for i, a in enumerate(uids)
                    ]
                )
        # centroid in the epoch where the source is brightest (catalogue S/N), then measure all
        ks = sorted(files_k for files_k in scales)
        best = np.zeros(n, int)
        snr_cols = (
            np.column_stack(
                [
                    np.asarray(fl[f"{b}_e{k}_flux"]) / np.asarray(fl[f"{b}_e{k}_err"])
                    if f"{b}_e{k}_flux" in fl.colnames
                    else np.full(len(fl), np.nan)
                    for k in ks
                ]
            )
            if len(fl)
            else np.zeros((0, len(ks)))
        )
        with _quiet():
            best[: len(fl)] = (
                np.asarray(ks)[np.nanargmax(np.nan_to_num(snr_cols, nan=-1e9), axis=1)]
                if len(fl)
                else []
            )
        best[len(fl) :] = ks[0]
        ra_c, dec_c = np.asarray(targets["ra"], float), np.asarray(targets["dec"], float)
        for k in ks:
            sel = best == k
            if sel.any():
                _, _, r_, d_ = measure(
                    [f for f, s in zip(files[(k, b)], sel, strict=True) if s],
                    0.15,
                    scales[k],
                    ra_c[sel],
                    dec_c[sel],
                    recentre_arcsec=0.1,
                )
                ok = np.isfinite(r_) & np.isfinite(d_)
                idx = np.flatnonzero(sel)
                ra_c[idx[ok]], dec_c[idx[ok]] = r_[ok], d_[ok]
        for k in ks:
            f_, e_, _, _ = measure(files[(k, b)], 0.15, scales[k], ra_c, dec_c)
            flux[:, k, jb], ferr[:, k, jb] = f_, e_
    # zero points from the controls, per band against the epoch with most measurable controls
    for jb in range(len(bands)):
        good = np.isfinite(flux[is_ctl, :, jb]).sum(axis=0)
        if good.max() == 0:
            continue
        rk = int(np.argmax(good))
        for k in range(ne):
            with np.errstate(invalid="ignore", divide="ignore"):
                c = is_ctl & (flux[:, k, jb] > 0) & (flux[:, rk, jb] > 0)
                ratio = flux[c, rk, jb] / flux[c, k, jb]
            if c.sum() >= 10:
                flux[:, k, jb] *= np.median(ratio)
                ferr[:, k, jb] *= np.median(ratio)
    scale = noise_scale(flux[is_ctl], ferr[is_ctl], np.isfinite(flux[is_ctl]), min_n=50)
    err = scaled_errors(flux, ferr, scale, p.sys_floor)
    res = classify(flux, err, p)
    out = Table({"uid": uids, "ra": targets["ra"], "dec": targets["dec"], "control": is_ctl})
    for jb, b in enumerate(bands):
        for k in range(ne):
            out[f"{b}_e{k}_fflux"] = flux[:, k, jb]
            out[f"{b}_e{k}_ferr"] = err[:, k, jb]
    for name, v in res.items():
        out[f"forced_{name}"] = v
    out["forced_flagged"] = flagged(res)
    for name, v in qual.items():
        out[name] = v
    # same flag type as the catalogue stage must reappear in forced photometry
    same = np.zeros(n, bool)
    for kname in FLAGS:
        cat_k = np.concatenate([np.asarray(fl[kname], bool), np.zeros(int(is_ctl.sum()), bool)])
        same |= cat_k & res[kname]
    out["forced_confirmed"] = same & ~is_ctl
    for name in ("near_star", "edge_proxy", "blended", "single_epoch", "point_like"):
        out[name] = np.concatenate([np.asarray(fl[name], bool), np.zeros(int(is_ctl.sum()), bool)])
    cut_ok = ~(out["cut_edge"] | out["cut_low_weight"] | out["cut_spike"])
    cat_ok = ~(out["near_star"] | out["edge_proxy"] | out["blended"] | out["single_epoch"])
    out["survives_artefact_tests"] = out["forced_confirmed"] & cut_ok & cat_ok
    surv = np.flatnonzero(out["survives_artefact_tests"])
    xm_note = "not run (no survivors)"
    out["known_object"] = np.zeros(n, "U40")
    if surv.size and not args.no_crossmatch:
        tg = Table({"source_uid": out["uid"][surv], "ra": out["ra"][surv], "dec": out["dec"][surv]})
        xm = crossmatch.crossmatch(tg, radius_arcsec=1.0, services=("simbad", "ned"))
        xm_note = "SIMBAD and NED, 1 arcsec (crossmatch.crossmatch, batched)"
        for i, r in zip(surv, xm, strict=True):
            names = [str(r[c]) for c in xm.colnames if c.endswith("_otype") or c == "otype"]
            out["known_object"][i] = ";".join(x for x in names if x and x != "--")[:40]
    out.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source=f"scripts/dimming_screen.py forced --field {args.field}",
        noise_scale={b: float(s) for b, s in zip(bands, scale, strict=True)},
        aperture_arcsec=0.15,
        recentre_arcsec=0.1,
        crossmatch=xm_note,
        unread_images=skipped,
        params=asdict(p),
    )
    out.write(args.out / "forced.ecsv", overwrite=True)
    show = np.flatnonzero(out["forced_confirmed"])
    if show.size:
        _epoch_sheet(out, show, files, field, args.out / "forced_sheet.png")
    summary = {
        "field": args.field,
        "n_flags_measured": int((~is_ctl).sum()),
        "n_controls": int(is_ctl.sum()),
        "noise_scale": out.meta["noise_scale"],
        "forced_confirmed": int(out["forced_confirmed"].sum()),
        "confirmed_by_flag": {
            k: int((out["forced_confirmed"] & out[f"forced_{k}"]).sum()) for k in FLAGS
        },
        "confirmed_and_cutout_clean": int((out["forced_confirmed"] & cut_ok).sum()),
        "survives_artefact_tests": int(out["survives_artefact_tests"].sum()),
        "controls_flagged": int((out["forced_flagged"] & is_ctl).sum()),
        "crossmatch": xm_note,
        "unread_images": skipped,
        "wall_s": round(time.time() - t_start, 1),
    }
    (args.out / "forced_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


def _epoch_sheet(out: Table, rows: np.ndarray, files, field, png: Path, max_rows: int = 40):
    """One row per source, one column per (band, epoch): the forced-photometry cutouts."""
    from astropy.io import fits
    from matplotlib.figure import Figure

    from jwst_anomaly.viz import normalize

    keys = sorted(files, key=lambda kb: (field["bands"].index(kb[1]), kb[0]))
    rows = rows[:max_rows]
    fig = Figure(figsize=(1.3 * len(keys), 1.4 * len(rows) + 0.5), dpi=90, layout="constrained")
    axes = fig.subplots(len(rows), len(keys), squeeze=False)
    for i, r in enumerate(rows):
        for j, kb in enumerate(keys):
            ax = axes[i, j]
            ax.set_xticks([])
            ax.set_yticks([])
            path = files[kb][r]
            if path:
                d = fits.getdata(path, extname="SCI")
                ax.imshow(d, norm=normalize(d), origin="lower", cmap="gray")
            if i == 0:
                ax.set_title(f"{kb[1]} e{kb[0]}", fontsize=7)
            if j == 0:
                ax.set_ylabel(str(out["uid"][r])[-6:], fontsize=7)
    fig.savefig(png)


def combine_limits(per_field: dict[str, Table]) -> Table:
    """Joint 95 % limits over fields with zero surviving events: exposures add (``derived``).

    Per source per year: 3 / Σ exposure_source_yr. Per source per epoch pair: 3 / Σ n_effective
    n_pairs. Per deg² per epoch pair: 3 / Σ area n_pairs ε̄.
    """
    rows = []
    first = next(iter(per_field.values()))
    for r0 in first:
        t_e, rho = float(r0["t_e_yr"]), float(r0["rho"])
        exp_yr = pair = area = 0.0
        n_mon = 0
        for lim in per_field.values():
            r = lim[(lim["t_e_yr"] == t_e) & (lim["rho"] == rho)][0]
            n_pairs = float(lim.meta["n_pairs"])
            exp_yr += float(r["exposure_source_yr"])
            pair += float(r["n_effective"]) * n_pairs
            area += float(lim.meta["area_deg2"]) * n_pairs * float(r["mean_efficiency"])
            n_mon += int(r["n_monitored"])
        rows.append(
            {
                "t_e_yr": t_e,
                "rho": rho,
                "n_monitored": n_mon,
                "exposure_source_yr": exp_yr,
                "limit_per_source_per_yr": 3.0 / exp_yr if exp_yr > 0 else np.inf,
                "limit_per_source_per_pair": 3.0 / pair if pair > 0 else np.inf,
                "limit_per_deg2_per_pair": 3.0 / area if area > 0 else np.inf,
            }
        )
    return Table(rows=rows, meta={"provenance": "derived", "fields": list(per_field)})


def cmd_combine(args) -> dict:
    fields = args.field.split(",")
    root = args.out.parent
    per = {f: Table.read(root / f / "limits.ecsv") for f in fields}
    out = combine_limits(per)
    out.write(root / "limits_combined.ecsv", overwrite=True)
    return {"fields": fields, "rows": [dict(zip(out.colnames, r, strict=True)) for r in out]}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("command", choices=("fetch", "screen", "inject", "forced", "combine"))
    ap.add_argument("--field", required=True, help="a field; combine: comma-separated fields")
    ap.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--no-gaia", action="store_true", help="skip the Gaia star mask (offline)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument(
        "--max-flags",
        type=int,
        default=0,
        help="forced: measure every point-like flag and a random audit sample of the rest, "
        "this many in total (0 = all)",
    )
    ap.add_argument("--n-controls", type=int, default=300)
    ap.add_argument("--control-mag", nargs=2, type=float, default=(24.0, 27.5))
    ap.add_argument("--no-crossmatch", action="store_true")
    args = ap.parse_args(argv)
    args.out = args.out or paths.outputs_dir() / "dimming" / args.field
    fn = {
        "fetch": cmd_fetch,
        "screen": cmd_screen,
        "inject": cmd_inject,
        "forced": cmd_forced,
        "combine": cmd_combine,
    }
    print(json.dumps(fn[args.command](args), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
