"""Star/galaxy separation so stars are ranked apart from galaxies (decision D-012).

On program 2736, bright stars (saturated cores, diffraction spikes) filled the top of the
anomaly ranking, and the pipeline's ``is_extended`` flag calls them extended. External
catalogs identify them: like ``crossmatch``'s ``is_star``, the *nearest* SIMBAD/Gaia DR3
match decides. A source is a ``star`` when that match is a star (SIMBAD otype or Gaia
astrometry), or when it is a Gaia DR3 source within ``gaia_radius_arcsec`` (in extragalactic
fields Gaia detections are mostly stars; a nearer SIMBAD galaxy wins). Everything else is
``other``. Stars are not discarded: the runner ranks them as their own stratum.
"""

from __future__ import annotations

import numpy as np
from astropy.table import Table

from jwst_anomaly import schema
from jwst_anomaly.crossmatch import STAR_SERVICES

# Query used by the runner: one bulk cross-match of every source.
DEFAULT_QUERY_RADIUS_ARCSEC = 0.5
DEFAULT_SERVICES = ("gaia", "simbad")
# ASSUMPTION (D-012): a Gaia DR3 source this close is the JWST source itself.
DEFAULT_GAIA_RADIUS_ARCSEC = 0.3


def classify_sources(
    sources: Table,
    matches: Table,
    *,
    gaia_radius_arcsec: float = DEFAULT_GAIA_RADIUS_ARCSEC,
) -> Table:
    """Population per merged source from ``crossmatch.query_matches`` output (long format).

    Returns ``schema.CLASSIFY_COLUMNS``: ``population`` is ``star`` or ``other`` and
    ``star_basis`` names the evidence (``simbad_or_gaia_astrometry``, ``gaia_position``, or "").
    Raises ``ValueError`` when a service that feeds the star rule (``crossmatch.STAR_SERVICES``)
    failed, so a partial cross-match never passes for a complete classification, or when
    ``gaia_radius_arcsec`` exceeds the query radius recorded in ``matches.meta``.
    """
    schema.validate(sources, schema.SOURCE_COLUMNS, name="sources")
    failed = sorted(set(matches.meta.get("services_failed") or []) & set(STAR_SERVICES))
    if failed:
        raise ValueError(f"cross-match services failed: {failed}")
    query_radius = matches.meta.get("radius_arcsec")
    if query_radius is not None and gaia_radius_arcsec > float(query_radius):
        raise ValueError(
            f"gaia_radius_arcsec {gaia_radius_arcsec} exceeds the query radius {query_radius}"
        )
    uids = [str(u) for u in sources["source_uid"]]
    nearest: dict[str, tuple[float, str, bool]] = {}  # uid -> (sep, service, is_star)
    if len(matches):
        is_star = np.asarray(np.ma.filled(matches["is_star"], False), dtype=bool)
        sep = np.asarray(np.ma.filled(matches["sep_arcsec"], np.inf), dtype=float)
        service = np.asarray(matches["service"]).astype(str)
        uid_col = np.asarray(matches["source_uid"]).astype(str)
        for uid, star, s, svc in zip(uid_col, is_star, sep, service, strict=True):
            if svc not in STAR_SERVICES:
                continue
            if uid not in nearest or s < nearest[uid][0]:
                nearest[uid] = (float(s), svc, bool(star))
    basis = []
    for uid in uids:
        s, svc, star = nearest.get(uid, (np.inf, "", False))
        if star:
            basis.append("simbad_or_gaia_astrometry")
        elif svc == "gaia" and s <= gaia_radius_arcsec:
            basis.append("gaia_position")
        else:
            basis.append("")
    out = Table(
        {
            "source_uid": sources["source_uid"],
            "population": np.array(["star" if b else "other" for b in basis], dtype="U8"),
            # Fixed width so later bases (e.g. "stellar_locus", D-015) are never truncated.
            "star_basis": np.array(basis, dtype="U32"),
        }
    )
    services = matches.meta.get("services_requested") or sorted(
        set(np.asarray(matches["service"]).astype(str))
    )
    out.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source=f"merged sources ({len(uids)}) cross-matched with {', '.join(map(str, services))}",
        thresholds={
            "gaia_radius_arcsec": gaia_radius_arcsec,
            "query_radius_arcsec": query_radius,
            "rule": "nearest SIMBAD/Gaia match decides (as crossmatch is_star)",
            "provenance": schema.Provenance.ASSUMPTION.value,
        },
    )
    return schema.validate(out, schema.CLASSIFY_COLUMNS, name="populations")


# ASSUMPTIONS (D-015): the stellar locus in a matched-photometry catalog's detection image.
DEFAULT_LOCUS = {
    "calib_mag_range": (20.0, 22.5),  # unsaturated catalogued stars calibrate the locus
    "r50_tolerance": 0.2,  # |r50 / r50_psf - 1| <= this
    "mag_max": 24.0,  # fainter, point-like galaxies and stars cannot be told apart by size
    "colours": (("f150w", "f444w"), ("f200w", "f356w")),
    "colour_pad": 0.3,  # mag beyond the catalogued stars' 5-95% colour range
    "min_ref_stars": 10,
}


def stellar_locus(
    sources: Table, label: str, known_stars: np.ndarray, **overrides
) -> tuple[np.ndarray, dict]:
    """Point-like sources with ordinary stellar colours, calibrated on catalogued stars (D-015).

    ``sources`` carries the matched-photometry columns of ``label`` (``<label>_r50_pix``,
    ``<label>_mag_auto``, ``<band>_<label>_abmag``). Catalogued stars (``known_stars``) with
    ``calib_mag_range`` magnitudes give the point-source half-light radius ``r50_psf`` and the
    stellar colour range. A source is in the locus when its r50 is within ``r50_tolerance`` of
    ``r50_psf``, it is brighter than ``mag_max``, and every configured colour lies within the
    stars' 5-95% range widened by ``colour_pad``. Point-like sources with unusual colours (brown
    dwarfs, compact high-z galaxies) are deliberately left out: they stay rankable as galaxies.
    Raises ``ValueError`` with fewer than ``min_ref_stars`` calibration stars.
    """
    cfg = {**DEFAULT_LOCUS, **overrides}
    r50 = np.asarray(sources[f"{label}_r50_pix"], dtype=float)
    mag = np.asarray(sources[f"{label}_mag_auto"], dtype=float)
    lo_m, hi_m = cfg["calib_mag_range"]
    with np.errstate(invalid="ignore"):
        ref = np.asarray(known_stars, dtype=bool) & (mag > lo_m) & (mag < hi_m) & np.isfinite(r50)
    if ref.sum() < cfg["min_ref_stars"]:
        raise ValueError(
            f"stellar locus needs >= {cfg['min_ref_stars']} catalogued stars with "
            f"{lo_m} < mag < {hi_m}; found {int(ref.sum())}"
        )
    r50_psf = float(np.median(r50[ref]))
    with np.errstate(invalid="ignore"):
        member = (np.abs(r50 / r50_psf - 1.0) <= cfg["r50_tolerance"]) & (mag < cfg["mag_max"])
    ranges = {}
    for blue, red in cfg["colours"]:
        colour = np.asarray(sources[f"{blue}_{label}_abmag"], float) - np.asarray(
            sources[f"{red}_{label}_abmag"], float
        )
        lo, hi = np.nanpercentile(colour[ref], [5, 95])
        pad = cfg["colour_pad"]
        with np.errstate(invalid="ignore"):
            member &= (colour >= lo - pad) & (colour <= hi + pad)
        ranges[f"{blue}-{red}"] = [round(float(lo), 3), round(float(hi), 3)]
    info = {
        "r50_psf_pix": round(r50_psf, 3),
        "n_calibration_stars": int(ref.sum()),
        "colour_ranges": ranges,
        "thresholds": {k: (list(v) if isinstance(v, tuple) else v) for k, v in cfg.items()},
        "provenance": schema.Provenance.ASSUMPTION.value,
    }
    return member, info
