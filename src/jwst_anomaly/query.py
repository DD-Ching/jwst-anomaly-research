"""MAST metadata/query layer. Owner: bootstrap unit 1 (archive access)."""

from __future__ import annotations

from collections.abc import Sequence

from astropy.table import Table


def query_observations(**criteria) -> Table:
    """Query public JWST observations on MAST.

    ``criteria`` are MAST CAOM constraints, e.g. ``proposal_id="2736"``,
    ``instrument_name="NIRCAM/IMAGE"``, ``filters="F200W"``, ``calib_level=3``,
    ``dataproduct_type="image"``. Returns one row per observation with at least
    ``schema.OBSERVATION_COLUMNS``, restricted to ``dataRights == "PUBLIC"`` unless the
    caller explicitly asks otherwise. ``meta["provenance"] == "observed"``.
    """
    raise NotImplementedError("bootstrap unit 1: archive access")


def list_products(
    observations: Table,
    subgroups: Sequence[str] = ("CAT", "I2D"),
    calib_level: int = 3,
) -> Table:
    """List data products of ``observations`` filtered by product subgroup and calib level.

    MAST product lists for a level-3 observation also contain its level-2 members;
    those must be excluded. Returns ``schema.PRODUCT_COLUMNS``.
    """
    raise NotImplementedError("bootstrap unit 1: archive access")
