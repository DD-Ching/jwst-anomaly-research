"""Baseline anomaly scoring and ranking. Owner: bootstrap unit 3."""

from __future__ import annotations

from collections.abc import Sequence

from astropy.table import Table


def score_anomalies(
    features: Table,
    methods: Sequence[str] = ("robust_z", "isolation_forest", "lof"),
    random_state: int = 0,
) -> Table:
    """Score each source; higher ``score`` = more anomalous, ``rank`` 1 = most anomalous.

    Returns ``schema.SCORE_COLUMNS`` plus ``score_<method>`` columns and a
    ``top_features`` attribution string. ``meta["provenance"] == "model_prediction"``.
    """
    raise NotImplementedError("bootstrap unit 3: features + baseline ranking")
