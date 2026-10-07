"""Star/galaxy separation so stars are ranked apart from galaxies (decision D-012).

On program 2736, bright stars (saturated cores, diffraction spikes) filled the top of the
anomaly ranking, and the pipeline's ``is_extended`` flag calls them extended. External
catalogs identify them reliably: a source is a ``star`` when SIMBAD or Gaia DR3 astrometry
says so (``crossmatch`` rules), or when any Gaia DR3 source lies within ``gaia_radius_arcsec``
(in extragalactic fields Gaia detections are mostly stars). Everything else is ``other``.
Stars are not discarded: the runner ranks them as their own stratum.
"""

from __future__ import annotations

import numpy as np
from astropy.table import Table

from jwst_anomaly import schema

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
    Raises ``ValueError`` when a requested service failed, so a partial cross-match never
    passes for a complete classification.
    """
    schema.validate(sources, schema.SOURCE_COLUMNS, name="sources")
    failed = list(matches.meta.get("services_failed") or [])
    if failed:
        raise ValueError(f"cross-match services failed: {failed}")
    uids = [str(u) for u in sources["source_uid"]]
    basis = dict.fromkeys(uids, "")
    if len(matches):
        is_star = np.asarray(np.ma.filled(matches["is_star"], False), dtype=bool)
        sep = np.asarray(np.ma.filled(matches["sep_arcsec"], np.inf), dtype=float)
        service = np.asarray(matches["service"]).astype(str)
        for uid, star, s, svc in zip(
            np.asarray(matches["source_uid"]).astype(str), is_star, sep, service, strict=True
        ):
            if uid not in basis:
                continue
            if star:
                basis[uid] = "simbad_or_gaia_astrometry"
            elif svc == "gaia" and s <= gaia_radius_arcsec and not basis[uid]:
                basis[uid] = "gaia_position"
    out = Table(
        {
            "source_uid": sources["source_uid"],
            "population": ["star" if basis[u] else "other" for u in uids],
            "star_basis": [basis[u] for u in uids],
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
            "provenance": schema.Provenance.ASSUMPTION.value,
        },
    )
    return schema.validate(out, schema.CLASSIFY_COLUMNS, name="populations")
