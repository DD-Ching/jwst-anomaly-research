"""Image cutouts and quality flags for ranked targets. Owner: bootstrap unit 4."""

from __future__ import annotations

from pathlib import Path

from astropy.table import Table


def make_cutouts(
    image_uri: str,
    targets: Table,
    size_arcsec: float = 3.0,
    out_dir: Path | None = None,
) -> Table:
    """Cut ``targets`` (``schema.TARGET_COLUMNS``) out of one level-3 ``_i2d.fits``.

    ``image_uri`` may be a local path or a cloud URI (``s3://stpubdata/...``); cloud
    reads should fetch only the needed bytes. Returns ``schema.CUTOUT_COLUMNS``.
    """
    raise NotImplementedError("bootstrap unit 4: cutouts + visualization")
