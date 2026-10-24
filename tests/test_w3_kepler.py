"""Offline tests of the frozen rules in scripts/w3_kepler.py (results/w3_kepler/PLAN.md)."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import w3_kepler as wk  # noqa: E402

from jwst_anomaly import exotic_sim as es  # noqa: E402


def _star(n_q=2, n=400, quality=None):
    t = np.concatenate([q * 100.0 + np.arange(n) * wk.CAD for q in range(n_q)])
    st = {
        "TIME": t,
        "CADENCENO": np.arange(len(t), dtype=float),
        "SAP_FLUX": np.full(len(t), 1000.0),
        "SAP_FLUX_ERR": np.full(len(t), 0.1),
        "PDCSAP_FLUX": np.full(len(t), 1000.0),
        "SAP_QUALITY": np.zeros(len(t)) if quality is None else quality,
        "QUARTER": np.repeat(np.arange(1, n_q + 1, dtype=float), n),
    }
    return st


def test_good_mask_bits():
    q = np.zeros(800)
    q[:5] = [1, 16, 512, 2048, 8192]
    m = wk.good_mask(_star(quality=q))
    assert list(m[:5]) == [False, False, False, True, True]


def test_w3_table_matches_direct_average():
    t = np.linspace(-2.2, 2.2, 301)
    for rho in (0.01, 0.1):
        f = wk.w3_flux_avg(t, 0.0, 1.0, 0.3, rho, 1.0)
        sub = t[:, None] + np.linspace(-wk.CAD / 2, wk.CAD / 2, wk.NSUB)[None, :]
        d = es.light_curve(sub.ravel(), 0.0, 1.0, 0.3, 1.0, -1, rho).reshape(sub.shape).mean(1)
        assert np.abs(f - d).max() < 2e-3


def test_inject_sap_noise_rule():
    rng = np.random.default_rng(0)
    sap, err = np.full(20000, 1000.0), np.full(20000, 1.0)
    s1, e1 = wk.inject_sap(sap, err, np.ones(20000), rng)
    assert np.array_equal(s1, sap) and np.array_equal(e1, err)
    s0, e0 = wk.inject_sap(sap, err, np.zeros(20000), rng)
    assert abs(s0.mean()) < 0.05 and np.allclose(e0, 0.3)
    assert abs(s0.std() - 0.3) < 0.01


def _screen_with_dip(depth_rel, at_quarter_end=False, shoulders=False):
    st = _star(n_q=2, n=400)
    sap = st["SAP_FLUX"].copy()
    lo, hi = (340, 400) if at_quarter_end else (150, 170)
    sap[lo:hi] = 1000.0 * depth_rel
    if shoulders:
        sap[lo - 2 : lo] = 1300.0
        sap[hi : hi + 2] = 1300.0
    return wk.screen_star(st, sap, st["SAP_FLUX_ERR"], {})


def test_box_dip_fast_edges_but_no_shoulder():
    res = _screen_with_dip(0.1)
    assert len(res) == 1
    r = res[0]
    assert r["edge_class"] == "fast" and r["pass_no_shoulder"] and not r["pass"]
    assert r["achromatic"] == "not tested (single band)"


def test_shoulders_on_both_sides_pass():
    r = _screen_with_dip(0.1, shoulders=True)[0]
    assert r["shoulder_ok"] and r["pass"]


def test_edge_at_quarter_end_not_bounded():
    r = _screen_with_dip(0.1, at_quarter_end=True)[0]
    assert r["quarter_start"] == 1 and r["edge_pre_class"] == "fast"
    assert r["edge_post_class"] == "insufficient" and not r["pass_no_shoulder"]


def test_common_mode_control_fails_by_construction():
    st = _star(n_q=2, n=400)
    rng = np.random.default_rng(1)
    good = wk.good_mask(st)
    f, over = wk.control_flux("common", st, good, rng)
    sap, err = wk.inject_sap(st["SAP_FLUX"], st["SAP_FLUX_ERR"], f, rng)
    res = wk.screen_star(st, sap, err, {}, over)
    assert res and not any(r["common_mode_ok"] for r in res)


def test_gate_rule():
    summ = {
        "w3": [{"recovered": 50, "n": 100}, {"recovered": 49, "n": 100}],
        "controls": {"a": {"passed": 10, "n": 100}},
    }
    g = wk.gate_of(summ)
    assert not g["all_cells_ge_0.50"] and g["all_controls_le_0.10"] and not g["pass"]
    summ["w3"][1]["recovered"] = 50
    assert wk.gate_of(summ)["pass"]
