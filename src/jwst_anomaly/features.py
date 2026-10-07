"""Derived per-source features for anomaly ranking. Owner: bootstrap unit 3."""

from __future__ import annotations

from astropy.table import Table


def build_features(sources: Table) -> Table:
    """Compute numeric features from a merged source table (``schema.SOURCE_COLUMNS``).

    Returns ``schema.FEATURE_ID_COLUMNS`` plus feature columns, with
    ``meta["feature_spec"]`` mapping each feature to a one-line definition.
    ``meta["provenance"] == "derived"``.
    """
    raise NotImplementedError("bootstrap unit 3: features + baseline ranking")
