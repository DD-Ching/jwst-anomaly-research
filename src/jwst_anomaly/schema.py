"""Column contracts and provenance labels shared by all pipeline stages.

Every stage returns an ``astropy.table.Table`` that contains at least the
columns listed for it here (extra columns are allowed) and sets:

* ``meta["provenance"]`` -- a :class:`Provenance` value (its ``.value`` string);
* ``meta["source"]`` -- short description of the inputs it was derived from
  (e.g. MAST data URIs, manifest path, upstream run id).

Per-band quantities in wide tables are named with :func:`band_column`,
e.g. ``f200w_aper50_abmag``.
"""

from __future__ import annotations

from collections.abc import Iterable
from enum import StrEnum

from astropy.table import Table


class Provenance(StrEnum):
    """What kind of knowledge a table or value represents (see docs/methodology.md)."""

    OBSERVED = "observed"  # archive data products as delivered (incl. pipeline catalogs)
    DERIVED = "derived"  # deterministic computations from observed data
    SIMULATED = "simulated"  # injected/synthetic data
    MODEL_PREDICTION = "model_prediction"  # outputs of fitted/learned models, incl. anomaly scores
    ASSUMPTION = "assumption"  # inputs taken as given without verification
    HYPOTHESIS = "hypothesis"  # proposed interpretations awaiting tests


# query.query_observations -> one row per MAST observation (CAOM column names kept as-is).
OBSERVATION_COLUMNS = (
    "obsid",
    "obs_id",
    "proposal_id",
    "instrument_name",
    "filters",
    "calib_level",
    "t_exptime",
    "s_ra",
    "s_dec",
    "dataRights",
)

# query.list_products -> one row per data product (MAST product-list column names).
PRODUCT_COLUMNS = (
    "obsID",
    "obs_id",
    "productFilename",
    "productSubGroupDescription",
    "dataURI",
    "size",
    "calib_level",
)

# acquire.fetch_products -> one row per local file; also persisted under data/manifests/.
MANIFEST_COLUMNS = (
    "dataURI",
    "productFilename",
    "local_path",  # relative to paths.data_root()
    "size",
    "sha256",
    "retrieved_utc",
    "pipeline_version",  # e.g. jwst CAL_VER / catalog meta; "" if unknown
)

# catalog.load_pipeline_catalog -> one band; pipeline columns kept, plus these.
BAND_CATALOG_COLUMNS = ("label", "ra", "dec")

# catalog.merge_bands -> one row per merged source; per-band columns via band_column().
SOURCE_COLUMNS = ("source_uid", "ra", "dec", "ref_band", "n_bands")

# features.build_features -> source_uid + numeric features; meta["feature_spec"] maps
# each feature name to a one-line definition.
FEATURE_ID_COLUMNS = ("source_uid",)

# quality.assess_sources -> one row per merged source; quality_ok gates ranking (D-011).
QUALITY_COLUMNS = (
    "source_uid",
    "rel_weight",
    "edge_dist_arcsec",
    "sharper_than_psf",
    "quality_ok",
    "quality_reason",
)

# rank.score_anomalies -> one row per source. score: higher = more anomalous;
# rank: 1 = most anomalous. Per-method scores as score_<method>; attributions in top_features.
SCORE_COLUMNS = ("source_uid", "score", "rank")

# Inputs to cutouts/crossmatch.
TARGET_COLUMNS = ("source_uid", "ra", "dec")

# cutouts.make_cutouts -> one row per (target, image).
CUTOUT_COLUMNS = ("source_uid", "band", "path", "frac_nan", "on_edge", "quality_flag")

# crossmatch.crossmatch -> one summary row per target.
XMATCH_COLUMNS = (
    "source_uid",
    "n_matches",
    "is_known_object",
    "is_star",
    "best_match_id",
    "best_match_type",
    "best_match_service",
    "best_match_sep_arcsec",
)


def band_column(band: str, quantity: str) -> str:
    """Name of a per-band column, e.g. ``band_column("F200W", "aper50_abmag")``."""
    return f"{band.lower()}_{quantity}"


def validate(table: Table, required: Iterable[str], *, name: str = "table") -> Table:
    """Check required columns and provenance metadata; return the table unchanged."""
    missing = [c for c in required if c not in table.colnames]
    if missing:
        raise ValueError(f"{name}: missing required columns {missing}")
    prov = table.meta.get("provenance")
    valid = {p.value for p in Provenance}
    if prov not in valid:
        raise ValueError(f"{name}: meta['provenance'] must be one of {sorted(valid)}, got {prov!r}")
    if not table.meta.get("source"):
        raise ValueError(f"{name}: meta['source'] must describe the inputs")
    return table
