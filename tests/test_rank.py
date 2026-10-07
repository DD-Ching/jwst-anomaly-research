"""Offline tests for rank.py: scoring contract, determinism, NaN handling, evaluation helpers."""

import inspect
import warnings

import numpy as np
import pytest
from astropy.table import Table

from jwst_anomaly import schema
from jwst_anomaly.rank import (
    METHODS,
    _expected_hits,
    feature_names,
    inject_outliers,
    injection_recovery,
    robust_location_scale,
    robust_zscores,
    score_anomalies,
    seed_stability,
)


def make_features(n=300, p=5, seed=0):
    rng = np.random.default_rng(seed)
    t = Table()
    t["source_uid"] = [f"s{i:04d}" for i in range(n)]
    for j in range(p):
        t[f"x{j}"] = rng.normal(0.0, 1.0, n)
    t.meta.update(
        provenance="derived",
        source="synthetic features",
        feature_spec={f"x{j}": f"standard normal {j}" for j in range(p)},
    )
    return t


def test_contract():
    f = make_features()
    s = score_anomalies(f)
    schema.validate(s, schema.SCORE_COLUMNS)
    assert s.meta["provenance"] == schema.Provenance.MODEL_PREDICTION.value
    assert "synthetic features" in s.meta["source"]
    for m in METHODS:
        assert f"score_{m}" in s.colnames
    assert "top_features" in s.colnames and "n_missing" in s.colnames
    assert sorted(s["rank"]) == list(range(1, len(f) + 1))
    assert list(s["rank"]) == list(range(1, len(f) + 1))  # rows sorted by rank
    assert np.all(np.diff(s["score"]) <= 0)
    assert np.all((s["score"] > 0) & (s["score"] <= 1))
    assert sorted(s["source_uid"]) == sorted(f["source_uid"])
    assert s.meta["feature_names"] == [f"x{j}" for j in range(5)]


def test_planted_outlier_ranks_at_the_top():
    f = make_features()
    f["x2"][17] = 12.0
    s = score_anomalies(f, random_state=0)
    assert s["source_uid"][0] == "s0017"
    assert s["top_features"][0].startswith("x2=+1")
    for m in ("robust_z", "lof"):
        assert s["source_uid"][np.argmax(s[f"score_{m}"])] == "s0017", m
    # Axis-parallel isolation is weaker for a single-feature outlier after clipping (D-004).
    by_if = np.asarray(s["source_uid"])[np.argsort(-s["score_isolation_forest"])]
    assert "s0017" in by_if[:3]


def test_deterministic_given_random_state():
    f = make_features()
    a = score_anomalies(f, random_state=3)
    b = score_anomalies(f, random_state=3)
    for c in a.colnames:
        assert np.array_equal(a[c], b[c]), c
    c = score_anomalies(f, random_state=4)
    a.sort("source_uid")
    c.sort("source_uid")
    assert np.array_equal(a["score_robust_z"], c["score_robust_z"])
    assert np.array_equal(a["score_lof"], c["score_lof"])
    assert not np.array_equal(a["score_isolation_forest"], c["score_isolation_forest"])


def test_method_subset_and_validation():
    f = make_features(n=50)
    s = score_anomalies(f, methods=("robust_z",))
    assert "score_robust_z" in s.colnames and "score_lof" not in s.colnames
    assert s.meta["methods"] == ["robust_z"]
    with pytest.raises(ValueError, match="methods"):
        score_anomalies(f, methods=("robust_z", "magic"))
    with pytest.raises(ValueError, match="methods"):
        score_anomalies(f, methods=())
    with pytest.raises(ValueError, match="at least 2"):
        score_anomalies(f[:1])


def test_nan_is_median_imputed_and_counted():
    f = make_features(n=100)
    f["x0"][:10] = np.nan
    f["x1"][0] = np.nan
    s = score_anomalies(f)
    for m in METHODS:
        assert np.all(np.isfinite(s[f"score_{m}"]))
    s.sort("source_uid")
    assert list(s["n_missing"][:3]) == [2, 1, 1]
    assert "x0=" not in s["top_features"][0] and "x1=" not in s["top_features"][0]


def test_lof_is_finite_with_many_duplicates():
    # Sources with the same missing pattern become identical after imputation.
    f = make_features(n=300)
    for name in feature_names(f):
        f[name][:150] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # sklearn warns on duplicates when LOF is ill-defined
        s = score_anomalies(f, methods=("lof",))
    assert np.max(s["score_lof"]) < 100


def test_lof_large_duplicate_group_is_an_inlier():
    rng = np.random.default_rng(0)
    x = np.vstack([rng.normal(0, 1, (300, 3)), np.full((200, 3), 4.0)])
    f = Table({"source_uid": [f"s{i:03d}" for i in range(500)]})
    for j in range(3):
        f[f"x{j}"] = x[:, j]
    f.meta.update(provenance="derived", source="duplicates probe")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        s = score_anomalies(f, methods=("lof",))
    top50 = set(s["source_uid"][:50])
    assert not top50 & {f"s{i:03d}" for i in range(300, 500)}


def test_ties_do_not_depend_on_the_seed():
    f = make_features(n=300)
    for i, (a, b) in enumerate([(1e3, 1e3), (-1e3, 1e3), (1e3, -1e3), (-1e3, -1e3)]):
        f["x0"][i], f["x1"][i] = a, b  # clipped at 10 -> exact ties in robust_z
    st = seed_stability(f, seeds=(0, 1, 2), k=2, score_kwargs={"n_estimators": 20})
    row = {m: r for m, r in zip(st["method"], st, strict=True)}
    assert row["robust_z"]["min_jaccard"] == 1.0 and row["lof"]["min_jaccard"] == 1.0


def test_expected_hits_shares_ties_pro_rata():
    score = np.array([3.0, 2.0, 2.0, 2.0, 1.0])
    positive = np.array([False, True, False, False, True])
    assert _expected_hits(score, positive, 2) == pytest.approx(1 / 3)
    assert _expected_hits(score, positive, 4) == pytest.approx(1.0)
    assert _expected_hits(score, positive, 5) == pytest.approx(2.0)


def test_evaluation_argument_validation():
    f = make_features(n=60)
    rec = injection_recovery(f, n_inject=5, k=np.int64(5), score_kwargs={"n_estimators": 20})
    assert list(rec["k"]) == [5] * 4
    with pytest.raises(ValueError, match="k must be"):
        injection_recovery(f, n_inject=5, k=0)
    with pytest.raises(ValueError, match="k must be"):
        seed_stability(f, k=0)
    with pytest.raises(ValueError, match="n_shift"):
        inject_outliers(f, 5, n_shift=0)


def test_inputs_need_provenance():
    f = make_features(n=20)
    f.meta.clear()
    with pytest.raises(ValueError, match="provenance"):
        score_anomalies(f)
    with pytest.raises(ValueError, match="provenance"):
        inject_outliers(f, 2)


def test_scores_of_injected_tables_stay_simulated():
    inj = inject_outliers(make_features(n=50), 3)
    s = score_anomalies(inj)
    assert s.meta["provenance"] == schema.Provenance.SIMULATED.value
    assert int(np.sum(s["is_injected"])) == 3
    assert set(s["provenance"][s["is_injected"]]) == {"simulated"}


def test_default_methods_match_methods_constant():
    default = inspect.signature(score_anomalies).parameters["methods"].default
    assert tuple(default) == METHODS


def test_feature_names_from_spec_or_numeric_columns():
    f = make_features(n=10)
    f["extra"] = np.arange(10.0)
    f["is_flag"] = np.zeros(10, bool)
    assert feature_names(f) == [f"x{j}" for j in range(5)]
    del f.meta["feature_spec"]
    assert feature_names(f) == [f"x{j}" for j in range(5)] + ["extra"]


def test_robust_location_scale_and_fallbacks():
    rng = np.random.default_rng(1)
    gauss = rng.normal(5.0, 2.0, 20000)
    mostly_zero = np.r_[np.zeros(980), np.ones(20)]  # MAD = IQR = 0 -> std
    x = np.column_stack(
        [gauss, np.full(20000, 3.0), np.tile(mostly_zero, 20), np.full(20000, np.nan)]
    )
    loc, scale = robust_location_scale(x)
    assert loc[0] == pytest.approx(5.0, abs=0.05) and scale[0] == pytest.approx(2.0, rel=0.03)
    assert scale[1] == 1.0  # constant column
    assert scale[2] == pytest.approx(np.std(mostly_zero))
    assert loc[3] == 0.0 and scale[3] == 1.0  # all-NaN column
    z = robust_zscores(x)
    assert np.all(np.isnan(z[:, 3]))


def test_inject_outliers_shift():
    f = make_features(n=200)
    inj = inject_outliers(f, 20, mode="shift", shift_sigma=7.0, n_shift=2, random_state=1)
    assert len(inj) == 220
    assert inj.meta["provenance"] == schema.Provenance.SIMULATED.value
    assert inj.meta["feature_spec"] == f.meta["feature_spec"]
    assert int(np.sum(inj["is_injected"])) == 20
    assert set(inj["provenance"][inj["is_injected"]]) == {"simulated"}
    assert set(inj["provenance"][~inj["is_injected"]]) == {"derived"}
    x = np.column_stack([np.asarray(f[n]) for n in feature_names(f)])
    loc, scale = robust_location_scale(x)
    fake = np.column_stack([np.asarray(inj[n])[200:] for n in feature_names(f)])
    z = (fake - loc) / scale
    assert np.all(np.sum(np.isclose(np.abs(z), 7.0), axis=1) == 2)
    again = inject_outliers(f, 20, mode="shift", shift_sigma=7.0, n_shift=2, random_state=1)
    assert np.array_equal(np.asarray(inj["x0"]), np.asarray(again["x0"]))


def test_inject_outliers_shuffle_keeps_marginals():
    f = make_features(n=100)
    inj = inject_outliers(f, 30, mode="shuffle", random_state=0)
    for n in feature_names(f):
        assert set(np.asarray(inj[n])[100:]) <= set(np.asarray(f[n]))
    with pytest.raises(ValueError, match="mode"):
        inject_outliers(f, 5, mode="bogus")


def test_injection_recovery_finds_strong_injections():
    f = make_features(n=500)
    rec = injection_recovery(f, n_inject=10, k=(10, 20), shift_sigma=9.0, random_state=0)
    assert rec.meta["provenance"] == schema.Provenance.SIMULATED.value
    assert set(rec["method"]) == {*METHODS, "ensemble"}
    assert len(rec) == 8
    assert np.all(rec["random_precision"] == pytest.approx(10 / 510))
    at10 = rec[rec["k"] == 10]
    assert np.all(at10["precision_at_k"] >= 0.8), at10
    assert np.all(rec[rec["k"] == 20]["recall_at_k"] >= 0.9)


def test_seed_stability():
    f = make_features(n=200)
    st = seed_stability(f, seeds=(0, 1, 2), k=10, score_kwargs={"n_estimators": 50})
    assert st.meta["provenance"] == schema.Provenance.MODEL_PREDICTION.value
    row = {m: r for m, r in zip(st["method"], st, strict=True)}
    assert row["robust_z"]["mean_jaccard"] == 1.0 and row["lof"]["min_jaccard"] == 1.0
    assert 0.0 <= row["isolation_forest"]["mean_jaccard"] <= 1.0
    assert row["ensemble"]["n_seeds"] == 3
    with pytest.raises(ValueError, match="two seeds"):
        seed_stability(f, seeds=(0,))
