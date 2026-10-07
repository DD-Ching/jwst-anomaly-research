"""Baseline anomaly scoring and ranking. Owner: bootstrap unit 3.

Design and rejected alternatives: DECISIONS.md D-004. All methods see the same input: per-feature
robust z-scores (median / normal-consistent MAD), NaN imputed as 0 (the median) and clipped at
``clip``. Methods are scikit-learn ``IsolationForest`` and ``LocalOutlierFactor`` plus a robust
diagonal-Mahalanobis score. The ensemble is the mean of per-method percentile ranks. The evaluation
helpers (:func:`inject_outliers`, :func:`injection_recovery`, :func:`seed_stability`) are how any
replacement model must show it beats this baseline (docs/methodology.md).

A high score means "unusual in this sample under this feature set", never a discovery.
"""

from __future__ import annotations

import operator
from collections.abc import Mapping, Sequence
from itertools import combinations
from typing import Any

import numpy as np
from astropy.table import Table
from scipy.stats import median_abs_deviation, rankdata
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor

from . import schema
from .features import column_as_float

METHODS = ("robust_z", "isolation_forest", "lof")
ENSEMBLE = "ensemble"


def feature_names(features: Table) -> list[str]:
    """Feature columns: ``meta["feature_spec"]`` keys if present, else numeric non-ID columns."""
    spec = features.meta.get("feature_spec")
    if spec:
        return [n for n in spec if n in features.colnames]
    skip = set(schema.FEATURE_ID_COLUMNS)
    return [c for c in features.colnames if c not in skip and features[c].dtype.kind in "fiu"]


def feature_matrix(features: Table, names: Sequence[str]) -> np.ndarray:
    """``(n_sources, n_features)`` float array; masked/non-finite values are NaN."""
    if not names:
        return np.empty((len(features), 0))
    return np.column_stack([column_as_float(features, n) for n in names])


def robust_location_scale(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per-column median and robust sigma, ignoring NaN.

    Sigma is the normal-consistent MAD (1.4826 x MAD). Where that is 0 (e.g. a mostly-constant
    integer feature) it falls back to IQR / 1.349, then the standard deviation, then 1.
    All-NaN columns get location 0 and scale 1.
    """
    p = x.shape[1]
    loc = np.zeros(p)
    scale = np.ones(p)
    for j in range(p):
        col = x[:, j][np.isfinite(x[:, j])]
        if col.size == 0:
            continue
        loc[j] = np.median(col)
        s = float(median_abs_deviation(col, scale="normal"))
        if not s > 0:
            q75, q25 = np.percentile(col, [75, 25])
            s = float(q75 - q25) / 1.349
        if not s > 0:
            s = float(np.std(col))
        scale[j] = s if s > 0 else 1.0
    return loc, scale


def robust_zscores(x: np.ndarray) -> np.ndarray:
    """``(x - median) / robust sigma`` per column (see :func:`robust_location_scale`); NaN kept."""
    loc, scale = robust_location_scale(x)
    return (x - loc) / scale


def _percentile_rank(score: np.ndarray) -> np.ndarray:
    """Fraction of sources with a score <= this one (ties averaged); 1.0 = most anomalous."""
    return rankdata(score, method="average") / score.size


def _attributions(z: np.ndarray, names: Sequence[str], n_top: int) -> list[str]:
    absz = np.where(np.isfinite(z), np.abs(z), -np.inf)
    order = np.argsort(-absz, axis=1, kind="stable")[:, :n_top]
    out = []
    for i, idx in enumerate(order):
        parts = [f"{names[j]}={z[i, j]:+.1f}" for j in idx if np.isfinite(z[i, j])]
        out.append("; ".join(parts))
    return out


def _lof_scores(zc: np.ndarray, n_neighbors: int) -> np.ndarray:
    """Local outlier factor of each row, with exact duplicates capped at ``n_neighbors`` copies.

    Median imputation turns sources that share a missing pattern and carry no other information
    into exact duplicates (groups of hundreds in real data). With more than ``n_neighbors`` copies
    LOF is undefined (zero k-distance, infinite density; neighbours get LOF ~ 1e7). Keeping at most
    ``n_neighbors`` copies of each distinct row keeps such a group maximally dense (an inlier) and
    every LOF finite; without such groups this is plain LOF. Copies share their distinct row's LOF.
    """
    k = min(n_neighbors, len(zc) - 1)
    uniq, inverse, counts = np.unique(zc, axis=0, return_inverse=True, return_counts=True)
    if len(uniq) < 2:
        return np.ones(len(zc))
    kept = np.minimum(counts, k)
    fit_rows = np.repeat(uniq, kept, axis=0)
    lof = LocalOutlierFactor(n_neighbors=min(k, len(fit_rows) - 1)).fit(fit_rows)
    first_copy = np.cumsum(kept) - kept
    return -lof.negative_outlier_factor_[first_copy][inverse.reshape(-1)]


def _check_k(k: int, n: int) -> int:
    k = operator.index(k)  # accepts numpy integers, rejects floats
    if not 1 <= k <= n:
        raise ValueError(f"k must be between 1 and the number of sources ({n}), got {k}")
    return k


def _check_methods(methods: Sequence[str]) -> list[str]:
    chosen = list(dict.fromkeys(methods))
    unknown = [m for m in chosen if m not in METHODS]
    if not chosen or unknown:
        raise ValueError(f"methods must be a non-empty subset of {METHODS}, got {list(methods)}")
    return chosen


def score_anomalies(
    features: Table,
    methods: Sequence[str] = ("robust_z", "isolation_forest", "lof"),
    random_state: int = 0,
    *,
    clip: float | None = 10.0,
    n_estimators: int = 1000,
    n_neighbors: int = 20,
    n_top_features: int = 3,
) -> Table:
    """Score each source; higher ``score`` = more anomalous, ``rank`` 1 = most anomalous.

    Returns ``schema.SCORE_COLUMNS`` plus ``score_<method>`` columns and a
    ``top_features`` attribution string. ``meta["provenance"] == "model_prediction"``.

    Per-method scores (higher = more anomalous):

    * ``robust_z``: root-mean-square of the clipped robust z-scores (diagonal robust Mahalanobis);
    * ``isolation_forest``: ``-IsolationForest.score_samples`` (``n_estimators`` trees);
    * ``lof``: local outlier factor, ``-negative_outlier_factor_`` (``n_neighbors`` neighbours),
      with exact duplicates capped at ``n_neighbors`` copies so imputed duplicates stay finite.

    ``score`` is the mean over methods of each method's percentile rank, in (0, 1]. Ties in
    ``score`` are broken by input order. ``top_features`` lists the ``n_top_features`` largest
    unclipped |robust z| as ``name=+z``: a model-agnostic attribution, not an exact decomposition
    of the forest or LOF scores. ``n_missing`` counts NaN (median-imputed) features per source.
    Rows are returned sorted by ``rank``. Deterministic for a given ``random_state``.

    If ``features`` contains simulated rows (``meta["provenance"] == "simulated"``, e.g. from
    :func:`inject_outliers`), the output keeps that label and its ``is_injected``/``provenance``
    columns, so synthetic rows can never pass as real candidates.
    """
    schema.validate(features, schema.FEATURE_ID_COLUMNS, name="features")
    chosen = _check_methods(methods)
    names = feature_names(features)
    if not names:
        raise ValueError("features: no numeric feature columns")
    n = len(features)
    if n < 2:
        raise ValueError(f"features: need at least 2 sources to rank, got {n}")

    x = feature_matrix(features, names)
    z = robust_zscores(x)
    zc = np.nan_to_num(z, nan=0.0)  # median imputation
    if clip is not None:
        zc = np.clip(zc, -clip, clip)

    per_method: dict[str, np.ndarray] = {}
    for m in chosen:
        if m == "robust_z":
            per_method[m] = np.sqrt(np.mean(zc**2, axis=1))
        elif m == "isolation_forest":
            forest = IsolationForest(n_estimators=n_estimators, random_state=random_state)
            per_method[m] = -forest.fit(zc).score_samples(zc)
        else:  # lof
            per_method[m] = _lof_scores(zc, n_neighbors)

    score = np.mean([_percentile_rank(s) for s in per_method.values()], axis=0)
    order = np.argsort(-score, kind="stable")
    rank = np.empty(n, dtype=int)
    rank[order] = np.arange(1, n + 1)

    out = Table()
    out["source_uid"] = features["source_uid"]
    out["score"] = score
    out["rank"] = rank
    for m, s in per_method.items():
        out[f"score_{m}"] = s
    out["n_missing"] = np.isnan(x).sum(axis=1)
    out["top_features"] = _attributions(z, names, n_top_features)
    simulated = features.meta["provenance"] == schema.Provenance.SIMULATED.value
    for passthrough in ("is_injected", "provenance"):
        if passthrough in features.colnames:
            out[passthrough] = features[passthrough]
    out = out[order]
    out.meta.update(
        provenance=(
            schema.Provenance.SIMULATED if simulated else schema.Provenance.MODEL_PREDICTION
        ).value,
        source=f"score_anomalies({', '.join(chosen)}) of: {features.meta.get('source', 'unknown')}",
        methods=chosen,
        random_state=random_state,
        params={
            "clip": clip,
            "n_estimators": n_estimators,
            "n_neighbors": n_neighbors,
            "n_top_features": n_top_features,
        },
        feature_names=names,
        ensemble_rule="mean of per-method percentile ranks",
    )
    return schema.validate(out, schema.SCORE_COLUMNS, name="scores")


def _score_columns(methods: Sequence[str]) -> dict[str, str]:
    return {**{m: f"score_{m}" for m in methods}, ENSEMBLE: "score"}


def _top_k(scores: Table, column: str, k: int) -> np.ndarray:
    """``source_uid`` of the k highest scores; ties broken by ``source_uid``, not by row order."""
    uid = np.asarray(scores["source_uid"]).astype(str)
    order = np.lexsort((uid, -np.asarray(scores[column], dtype=float)))
    return uid[order[:k]]


def _expected_hits(score: np.ndarray, positive: np.ndarray, k: int) -> float:
    """Positives among the k highest scores; ties at the k-th score count pro rata.

    This equals the expectation under random tie-breaking, so neither real nor injected rows
    win ties systematically (clipping makes exact ties common).
    """
    threshold = np.sort(score)[::-1][k - 1]
    above = score > threshold
    tied = score == threshold
    return float(positive[above].sum() + (k - above.sum()) * positive[tied].sum() / tied.sum())


def inject_outliers(
    features: Table,
    n_inject: int = 50,
    *,
    mode: str = "shift",
    shift_sigma: float = 6.0,
    n_shift: int = 2,
    random_state: int = 0,
) -> Table:
    """Append ``n_inject`` synthetic outliers to a feature table (for evaluation only).

    ``mode="shift"``: copy a random real source and set ``n_shift`` random continuous features to
    ``median +/- shift_sigma`` robust sigma (random sign). ``mode="shuffle"``: draw every feature
    from an independently chosen real source, which keeps the marginals and breaks correlations.

    Returns ``source_uid`` + feature columns + ``is_injected`` + a per-row ``provenance``
    (``derived`` or ``simulated``). The table's ``meta["provenance"]`` is ``simulated`` because it
    contains synthetic rows; ``meta["feature_spec"]`` is carried over.
    """
    schema.validate(features, schema.FEATURE_ID_COLUMNS, name="features")
    names = feature_names(features)
    x = feature_matrix(features, names)
    n, p = x.shape
    if n == 0 or p == 0:
        raise ValueError("features: need at least one source and one feature")
    if n_inject < 1:
        raise ValueError("n_inject must be >= 1")
    if mode == "shift" and not (n_shift >= 1 and shift_sigma > 0):
        raise ValueError(f"need n_shift >= 1 and shift_sigma > 0, got {n_shift}, {shift_sigma}")
    rng = np.random.default_rng(random_state)
    if mode == "shift":
        loc, scale = robust_location_scale(x)
        finite = np.where(np.isfinite(x), x, 0.0)
        continuous = [j for j in range(p) if np.any(finite[:, j] != np.round(finite[:, j]))]
        if not continuous:
            raise ValueError("mode='shift' needs at least one non-integer feature")
        fake = x[rng.integers(0, n, n_inject)].copy()
        k = min(n_shift, len(continuous))
        for row in fake:
            cols = rng.choice(continuous, size=k, replace=False)
            row[cols] = loc[cols] + rng.choice([-1.0, 1.0], size=k) * shift_sigma * scale[cols]
    elif mode == "shuffle":
        fake = x[rng.integers(0, n, size=(n_inject, p)), np.arange(p)]
    else:
        raise ValueError(f"mode must be 'shift' or 'shuffle', got {mode!r}")

    out = Table()
    uids = [str(u) for u in features["source_uid"]] + [f"injected-{i:05d}" for i in range(n_inject)]
    out["source_uid"] = uids
    for j, name in enumerate(names):
        out[name] = np.concatenate([x[:, j], fake[:, j]])
    out["is_injected"] = np.r_[np.zeros(n, bool), np.ones(n_inject, bool)]
    out["provenance"] = [schema.Provenance.DERIVED.value] * n + [
        schema.Provenance.SIMULATED.value
    ] * n_inject
    out.meta.update(
        provenance=schema.Provenance.SIMULATED.value,
        source=(
            f"inject_outliers(n_inject={n_inject}, mode={mode}, shift_sigma={shift_sigma}, "
            f"n_shift={n_shift}, random_state={random_state}) into: "
            f"{features.meta.get('source', 'unknown')}"
        ),
        feature_spec={nm: features.meta.get("feature_spec", {}).get(nm, nm) for nm in names},
    )
    return schema.validate(out, schema.FEATURE_ID_COLUMNS, name="injected features")


def injection_recovery(
    features: Table,
    methods: Sequence[str] = METHODS,
    *,
    n_inject: int = 50,
    k: int | Sequence[int] | None = None,
    mode: str = "shift",
    shift_sigma: float = 6.0,
    n_shift: int = 2,
    random_state: int = 0,
    score_kwargs: Mapping[str, Any] | None = None,
) -> Table:
    """Inject synthetic outliers, rank everything, report precision@k and recall@k per method.

    One row per (method, k), methods plus ``ensemble``. ``k`` defaults to ``n_inject``.
    ``random_precision`` is the expected precision@k of a random ranking; ties at the k-th score
    count pro rata. Real sources that are genuinely unusual count as false positives here, so
    these are lower bounds on usefulness.
    ``meta["provenance"] == "simulated"``: the numbers describe synthetic injections, not the sky.
    """
    chosen = _check_methods(methods)
    ks = [n_inject] if k is None else [k] if np.ndim(k) == 0 else list(k)
    ks = [_check_k(kk, len(features) + n_inject) for kk in ks]
    injected = inject_outliers(
        features,
        n_inject,
        mode=mode,
        shift_sigma=shift_sigma,
        n_shift=n_shift,
        random_state=random_state,
    )
    scores = score_anomalies(injected, chosen, random_state, **dict(score_kwargs or {}))
    positive = np.asarray(scores["is_injected"], dtype=bool)
    n_total = len(injected)
    rows = []
    for method, column in _score_columns(chosen).items():
        values = np.asarray(scores[column], dtype=float)
        for kk in ks:
            hits = _expected_hits(values, positive, kk)
            rows.append((method, kk, hits / kk, hits / n_inject))
    out = Table(rows=rows, names=("method", "k", "precision_at_k", "recall_at_k"))
    out["n_injected"] = n_inject
    out["n_total"] = n_total
    out["random_precision"] = n_inject / n_total
    out.meta.update(
        provenance=schema.Provenance.SIMULATED.value,
        source=f"injection_recovery: {injected.meta['source']}",
        mode=mode,
    )
    return out


def seed_stability(
    features: Table,
    methods: Sequence[str] = METHODS,
    *,
    seeds: Sequence[int] = (0, 1, 2, 3, 4),
    k: int = 20,
    score_kwargs: Mapping[str, Any] | None = None,
) -> Table:
    """Overlap of the top-``k`` sets across ``random_state`` seeds, per method and ensemble.

    ``mean_jaccard``/``min_jaccard`` are over all seed pairs; 1.0 means identical top-k sets.
    Ties are broken by ``source_uid``, so deterministic methods (``robust_z``, ``lof``) give 1.0.
    """
    chosen = _check_methods(methods)
    if len(seeds) < 2:
        raise ValueError("seed_stability needs at least two seeds")
    k = _check_k(k, len(features))
    columns = _score_columns(chosen)
    tops: dict[str, list[set[str]]] = {m: [] for m in columns}
    for seed in seeds:
        scores = score_anomalies(features, chosen, seed, **dict(score_kwargs or {}))
        for method, column in columns.items():
            tops[method].append({str(u) for u in _top_k(scores, column, k)})
    rows = []
    for method, sets in tops.items():
        jac = [len(a & b) / len(a | b) for a, b in combinations(sets, 2)]
        rows.append((method, k, len(seeds), float(np.mean(jac)), float(np.min(jac))))
    out = Table(rows=rows, names=("method", "k", "n_seeds", "mean_jaccard", "min_jaccard"))
    out.meta.update(
        provenance=schema.Provenance.MODEL_PREDICTION.value,
        source=f"seed_stability(seeds={list(seeds)}) of: {features.meta.get('source', 'unknown')}",
    )
    return out
