"""Cross-check targets against external catalogs. Owner: bootstrap unit 5."""

from __future__ import annotations

from collections.abc import Sequence

from astropy.table import Table


def crossmatch(
    targets: Table,
    radius_arcsec: float = 1.0,
    services: Sequence[str] = ("simbad", "ned", "gaia"),
) -> Table:
    """Match ``targets`` (``schema.TARGET_COLUMNS``) against external services.

    Returns one summary row per target with ``schema.XMATCH_COLUMNS``.
    ``meta["provenance"] == "observed"`` (external catalog content) with services and
    query dates recorded in ``meta``.
    """
    raise NotImplementedError("bootstrap unit 5: external cross-check")
