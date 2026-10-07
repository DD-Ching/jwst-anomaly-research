"""MAST metadata/query layer. Owner: bootstrap unit 1 (archive access).

A thin layer over ``astroquery.mast.Observations`` (DECISIONS.md D-002). It adds what this
project needs on top: PUBLIC-only results by default, level-3 product lists without the
level-2 member files MAST returns alongside them, and optional S3 URIs for byte-range reads.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Sequence
from typing import Any

import astroquery
import numpy as np
from astropy.table import Table, unique
from astropy.utils.introspection import minversion
from astroquery.mast import Observations
from astroquery.mast.utils import mast_relative_path

from . import schema

log = logging.getLogger(__name__)

COLLECTION = "JWST"
PUBLIC = "PUBLIC"
TOKEN_ENV = "MAST_API_TOKEN"  # optional; read from the environment only, never persisted
S3_BUCKET = "stpubdata"  # STScI public bucket (anonymous access)

_logged_in = False


def mast_client():
    """Return the shared ``astroquery.mast.Observations`` client.

    If ``$MAST_API_TOKEN`` is set, the client logs in with it once per process. The token is
    never logged, printed or written anywhere by this package.
    """
    global _logged_in
    token = os.environ.get(TOKEN_ENV, "").strip()
    if token and not _logged_in:
        Observations.login(token=token)
        _logged_in = True
    return Observations


def query_observations(*, public_only: bool = True, **criteria: Any) -> Table:
    """Query public JWST observations on MAST.

    ``criteria`` are MAST CAOM constraints, e.g. ``proposal_id="2736"``,
    ``instrument_name="NIRCAM/IMAGE"``, ``filters="F200W"``, ``calib_level=3``,
    ``dataproduct_type="image"``. Returns one row per observation with at least
    ``schema.OBSERVATION_COLUMNS``, restricted to ``dataRights == "PUBLIC"`` unless the
    caller explicitly asks otherwise. ``meta["provenance"] == "observed"``.

    ``obs_collection`` defaults to ``"JWST"``. Pass ``public_only=False`` to include
    exclusive-access observations (downloading those needs ``$MAST_API_TOKEN``).
    """
    if not set(criteria) - {"obs_collection", "dataRights"}:  # would list a whole collection
        raise ValueError("query_observations needs a constraint, e.g. proposal_id='2736'")
    criteria.setdefault("obs_collection", COLLECTION)
    if public_only:
        rights = criteria.setdefault("dataRights", PUBLIC)
        if set(np.atleast_1d(rights).tolist()) != {PUBLIC}:
            raise ValueError(f"dataRights={rights!r} contradicts public_only=True")

    obs = mast_client().query_criteria(**criteria)
    obs = _with_columns(obs, schema.OBSERVATION_COLUMNS)
    if public_only:  # the server already filtered; this guards against criteria being ignored
        obs = obs[str_values(obs["dataRights"]) == PUBLIC]
    obs.sort("obs_id")

    shown = ", ".join(f"{k}={v!r}" for k, v in sorted(criteria.items()))
    obs.meta["provenance"] = schema.Provenance.OBSERVED.value
    obs.meta["source"] = (
        f"MAST CAOM Observations.query_criteria({shown}); astroquery {astroquery.__version__}"
    )
    return schema.validate(obs, schema.OBSERVATION_COLUMNS, name="observations")


def list_products(
    observations: Table,
    subgroups: Sequence[str] = ("CAT", "I2D"),
    calib_level: int = 3,
    *,
    public_only: bool = True,
    cloud_uris: bool = False,
) -> Table:
    """List data products of ``observations`` filtered by product subgroup and calib level.

    MAST product lists for a level-3 observation also contain its level-2 members;
    those must be excluded. Returns ``schema.PRODUCT_COLUMNS``.

    A product is kept only if it belongs to one of ``observations`` (its own ``obsID``;
    members carry their own obsID with ``parent_obsid`` pointing at the level-3 one), its
    file name starts with that observation's ``obs_id + "_"``, and it matches ``calib_level``,
    ``subgroups`` and (by default) ``dataRights == "PUBLIC"``. All MAST product columns are
    kept (``prvversion`` is MAST's record of the pipeline version). With ``cloud_uris=True``
    a ``cloud_uri`` column holds ``s3://stpubdata/...`` keys from MAST's path-lookup service
    (one request per 50 products; existence on S3 is not checked; ``""`` if unknown).
    Rows are unique and sorted by ``dataURI``.
    """
    subgroups = [subgroups] if isinstance(subgroups, str) else list(subgroups)
    owner = dict(
        zip(str_values(observations["obsid"]), str_values(observations["obs_id"]), strict=True)
    )
    source = (
        f"MAST Observations.get_product_list({len(owner)} observations); "
        f"filter calib_level={calib_level}, productSubGroupDescription={subgroups}"
        f"{', dataRights=PUBLIC' if public_only else ''}; astroquery {astroquery.__version__}"
    )

    if owner:
        client = mast_client()
        products = client.get_product_list(sorted(owner))
        products = _with_columns(products, schema.PRODUCT_COLUMNS + ("dataRights",))
        filters: dict[str, Any] = {
            "calib_level": [int(calib_level)],  # astroquery's numeric filter needs a plain int
            "productSubGroupDescription": subgroups,
        }
        if public_only:
            filters["dataRights"] = [PUBLIC]
        products = client.filter_products(products, **filters)
        keep = [
            obsid in owner and fname.startswith(owner[obsid] + "_")
            for obsid, fname in zip(
                str_values(products["obsID"]),
                str_values(products["productFilename"]),
                strict=True,
            )
        ]
        products = products[np.asarray(keep, dtype=bool)]
    else:
        products = _with_columns(Table(), schema.PRODUCT_COLUMNS + ("dataRights",))

    if len(products):
        products = unique(products, keys="dataURI")
    products.sort("dataURI")
    if cloud_uris:
        uris = resolve_cloud_uris(products["dataURI"]) if len(products) else []
        public = str_values(products["dataRights"]) == PUBLIC
        products["cloud_uri"] = [u if ok else "" for u, ok in zip(uris, public, strict=True)]

    products.meta["provenance"] = schema.Provenance.OBSERVED.value
    products.meta["source"] = source
    return schema.validate(products, schema.PRODUCT_COLUMNS, name="products")


def resolve_cloud_uris(data_uris: Sequence[str]) -> list[str]:
    """Map MAST ``dataURI`` values to ``s3://stpubdata/...`` URIs (``""`` where MAST has none).

    Uses MAST's ``path_lookup`` service via astroquery (the lookup ``get_cloud_uris`` does
    internally) without its per-file boto3 existence check; see D-002. Needs astroquery
    >= 0.4.11: older versions take one URI (<= 0.4.9) or drop unknown ones from the result
    (0.4.10), which would misalign URIs.
    """
    data_uris = [str(u) for u in data_uris]
    if not data_uris:
        return []
    if not minversion("astroquery", "0.4.11"):
        raise RuntimeError(f"cloud URIs need astroquery>=0.4.11, found {astroquery.__version__}")
    keys = mast_relative_path(data_uris, verbose=False)
    if len(keys) != len(data_uris):
        raise RuntimeError(f"path_lookup returned {len(keys)} paths for {len(data_uris)} URIs")
    return [f"s3://{S3_BUCKET}/{key}" if key else "" for key in keys]


def str_values(col) -> np.ndarray:
    """Column values as a plain ``str`` array; masked entries (any dtype) become ``""``."""
    values = np.asarray(getattr(col, "data", col)).astype(str)
    return np.where(np.ma.getmaskarray(col), "", values)


def _with_columns(table: Table, columns: Sequence[str]) -> Table:
    """Return ``table``; if it is empty, add any missing ``columns`` (MAST omits them then)."""
    if len(table) == 0:
        table = Table(table, copy=False)
        for name in columns:
            if name not in table.colnames:
                dtype = int if name in ("calib_level", "size") else str
                table[name] = np.zeros(0, dtype=dtype)
    return table
