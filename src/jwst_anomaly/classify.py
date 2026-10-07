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
            "population": ["star" if b else "other" for b in basis],
            "star_basis": basis,
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
