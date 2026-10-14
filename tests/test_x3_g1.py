"""Offline tests for scripts/x3_g1.py and scripts/x3_bulk.py (X3 G1)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

_S = Path(__file__).resolve().parents[1] / "scripts"


def _load(name):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _S / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


g1 = _load("x3_g1")
xb = _load("x3_bulk")
w3 = g1.w3

P = {"t0": 500.0, "tE": 60.0, "u0": 0.5, "rho": 0.01}
T = np.linspace(0.0, 1000.0, 400)


def _band(fs, fb, p=P, noise=0.003, seed=0):
    rng = np.random.default_rng(seed)
    a = w3.magnification("N1neg", w3.LightCurve(T, T, np.ones_like(T)), p)
    f = fs * a + fb
    e = np.full_like(T, noise * (fs + fb))
    return T, f + rng.normal(0.0, e), e


def test_common_magnification_allows_any_blend_colour():
    # different blend fractions per band (blend colour differs from source colour): consistent
    for fs, fb in ((0.9, 0.1), (0.5, 0.5)):
        r = g1.band_test(*_band(fs, fb), P, 1.0)
        assert r["verdict"] == "consistent", r
        assert abs(r["fs"] - fs) < 0.05


def test_chromatic_dip_is_inconsistent_and_sparse_band_is_insufficient():
    t, f, e = _band(1.0, 0.0)
    dp = {"t0": 500.0, "width": 40.0, "depth_g": 1.0}
    fd, ed = g1.inject_dust(t, np.full_like(t, 1.0), e, dp, "BP")
    assert g1.band_test(t, fd, ed, P, 1.0)["verdict"] in ("inconsistent", "insufficient")
    far = T < 100.0  # no epoch inside the event
    assert g1.band_test(T[far], f[far], e[far], P, 1.0)["verdict"] == "insufficient"


def test_achromatic_combines_bands():
    bands = {"BP": _band(0.7, 0.3, seed=1), "RP": _band(0.95, 0.05, seed=2)}
    v, res = g1.achromatic(bands, P, {"BP": 1.0, "RP": 1.0})
    assert v == "consistent"
    bands["RP"] = (T[:5], T[:5] * 0 + 1.0, T[:5] * 0 + 0.01)
    assert g1.achromatic(bands, P, {})[0] == "consistent"  # one informative band suffices
    bands["BP"] = bands["RP"]
    assert g1.achromatic(bands, P, {})[0] == "insufficient"


def test_repeat_dip_needs_a_second_dip_not_a_period():
    t, f, e = _band(0.8, 0.2, noise=0.002)
    lc = w3.LightCurve(t, f, e)
    assert not g1.repeat_dip(lc, P, 0.8, 0.2)
    f2 = f.copy()
    f2[(t > 850) & (t < 870)] -= 0.5
    assert g1.repeat_dip(w3.LightCurve(t, f2, e), P, 0.8, 0.2)


def test_wilson_interval():
    lo, hi = g1.wilson(0, 1995)
    assert lo == 0.0 and 0.0015 < hi < 0.0025
    assert g1.wilson(0, 0) == [None, None]


def test_bulk_row_becomes_datalink_csv():
    from jwst_anomaly import gaia_mulens as gm

    row = {
        "source_id": "364105032180748928",
        "transit_id": "[1, 2, 3]",
        "g_transit_time": "[1700.5, NaN, 1720.25]",
        "g_transit_flux": "[100.0, NaN, 90.0]",
        "g_transit_flux_error": "[1.0, NaN, 1.0]",
        "g_transit_flux_over_error": "[100.0, NaN, 90.0]",
        "g_transit_mag": "[18.1, NaN, 18.2]",
        "g_transit_n_obs": "[9, 0, 8]",
        "bp_obs_time": "[1700.6, NaN, NaN]",
        "bp_flux": "[50.0, NaN, NaN]",
        "bp_flux_error": "[2.0, NaN, NaN]",
        "bp_mag": "[18.5, NaN, NaN]",
        "rp_obs_time": "",
        "variability_flag_g_reject": '["False", "True", "False"]',
        "variability_flag_bp_reject": '["False", "True", "False"]',
        "variability_flag_rp_reject": "",
    }
    text = xb.to_datalink_csv(row)
    lines = text.strip().splitlines()
    assert lines[0].split(",") == list(xb.OUT_COLS) and len(lines) == 4
    assert lines[2].endswith(",true")  # g_transit_n_obs == 0 -> rejected_by_photometry
    lc = gm.parse_epoch_csv(text)
    assert len(lc) == 2
    bands = g1.parse_bands(text)
    assert bands["BP"][0].size == 1 and bands["RP"][0].size == 0


def test_bulk_file_names_cover_level8_ranges():
    lo, hi = xb.file_range(xb.FILES[0][0])
    assert (lo, hi) == (40908, 41937)
    assert xb.level8(364105032180748928) == 364105032180748928 >> 43
    ranges = [xb.file_range(n) for n, _, _ in xb.FILES]
    assert all(b + 1 == c for (_, b), (c, _) in zip(ranges, ranges[1:], strict=False))
    assert "MOD(source_id, 10) <> 0" in xb.HELDOUT_QUERY


def test_vetting_order_and_visits():
    xv = _load("x3_vet")
    assert xv.visits(np.array([10.0, 10.1, 10.2])) == 1
    assert xv.visits(np.array([10.0, 10.1, 40.0])) == 2
    base = {
        "dip_visits": 3,
        "spike_visits": 0,
        "G": 16.0,
        "neighbour_dG": 99.0,
        "vari_class": "",
        "simbad_var": False,
        "fit_tE": 60.0,
        "fit_u0": 0.5,
    }
    assert xv.verdict(base) == "unexplained under A-D"
    assert xv.verdict({**base, "dip_visits": 1}) == "A1 single-visit"
    assert xv.verdict({**base, "G": 12.0}) == "A2 bright"
    assert xv.verdict({**base, "neighbour_dG": 1.0}) == "B1 neighbour"
    assert xv.verdict({**base, "vari_class": "ECL"}) == "C1 known variable"
    assert xv.verdict({**base, "fit_u0": 1.95}) == "D1 degenerate fit"
