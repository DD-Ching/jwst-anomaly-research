"""Ingest JWST pipeline source catalogs (``_cat.ecsv``). Owner: bootstrap unit 2."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

from astropy.table import Table


def load_pipeline_catalog(path: str | Path) -> Table:
    """Read one level-3 ``_cat.ecsv`` into a normalized single-band table.

    Keeps the pipeline columns, adds ``schema.BAND_CATALOG_COLUMNS`` (``ra``/``dec`` in deg),
    and records the band, pipeline/photutils versions and input file in ``meta``.
    ``meta["provenance"] == "observed"``.
    """
    raise NotImplementedError("bootstrap unit 2: catalog ingestion")


def merge_bands(
    catalogs: Mapping[str, Table],
    ref_band: str,
    radius_arcsec: float = 0.1,
) -> Table:
    """Positionally merge single-band catalogs ``{band: table}`` into one source table.

    Returns ``schema.SOURCE_COLUMNS`` plus per-band columns named with
    ``schema.band_column(band, quantity)`` and per-band match separations.
    ``meta["provenance"] == "derived"``.
    """
    raise NotImplementedError("bootstrap unit 2: catalog ingestion")
