"""Reproducible, cached acquisition with checksummed manifests. Owner: bootstrap unit 1."""

from __future__ import annotations

from pathlib import Path

from astropy.table import Table


def fetch_products(
    products: Table,
    data_root: Path | None = None,
    manifest_path: Path | None = None,
) -> Table:
    """Download ``products`` (``schema.PRODUCT_COLUMNS``) into the data cache.

    Idempotent: files already present with a matching sha256 are not re-downloaded.
    Writes/updates a manifest (``schema.MANIFEST_COLUMNS``) and returns its rows for
    ``products``. ``data_root`` defaults to ``paths.data_root()``; ``manifest_path``
    defaults to a file under ``paths.manifests_dir()``.
    """
    raise NotImplementedError("bootstrap unit 1: archive access")
