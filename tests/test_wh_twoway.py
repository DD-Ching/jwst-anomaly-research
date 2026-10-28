"""Offline tests for scripts/wh_twoway.py (small grids, a few seconds in total)."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import wh_twoway as wt  # noqa: E402


def _static_error(n: int) -> float:
    res = wt.evolve(wt.Params(n=n, u_max=3.0, v_max=3.0, a1=0.0, a2=0.0, out_h=0.1))
    r = res["grid"]["r"]
    u = np.arange(r.shape[0]) * res["grid_h"]
    uu, vv = np.meshgrid(u, u, indexing="ij")
    return float(np.nanmax(np.abs(r - np.sqrt(1.0 + ((vv - uu) / 2.0) ** 2))))


def test_static_wormhole_kept_to_second_order():
    e1, e2 = _static_error(20), _static_error(40)
    assert e2 < 1e-3
    assert 3.5 < e1 / e2 < 4.5


def test_pulse_profile_integrates_its_derivative():
    x = np.linspace(0.0, 4.0, 4001)
    d, val = wt.pulse(x, 0.3, 1.0, 2.5)
    num = np.concatenate([[0.0], np.cumsum(0.5 * (d[1:] + d[:-1]) * np.diff(x))])
    assert np.max(np.abs(num - val)) < 1e-6
    assert abs(d[0]) < 1e-20 and abs(d[-1]) < 1e-20


def test_amplitude_for_energy_hits_target():
    a = wt.amplitude_for_energy(1e-3, 3.0, 4.0, 100, 8.0)
    x = np.arange(801) / 100
    assert abs(wt._surface(x, a, 3.0, 4.0, 1.0)["energy"] / 1e-3 - 1.0) < 1e-6


def test_frozen_rule():
    assert wt.passes(1.0, 0.5)
    assert not wt.passes(1.0, 0.4999)
    assert not wt.passes(0.0, 0.0)


def _pulse_run(n: int) -> dict:
    # outbound pulse from the + side, return pulse from the - side; stations at r = 1.5
    p = wt.Params(n=n, u_max=6.0, v_max=8.0, a1=0.01, v1=2.5, v2=3.5, a2=0.01, u1=3.0, u2=4.0,
                  r_station=1.5, out_h=0.1)  # fmt: skip
    return wt.evolve(p)


@pytest.fixture(scope="module")
def pulse_run():
    return _pulse_run(80)


def _bookkeeping_error(res: dict) -> float:
    err = 0.0
    for side in ("A", "B"):
        t = res["stations"][side]
        cum = np.concatenate([[0.0], np.cumsum(t["dm_u"].sum(1) + t["dm_v"].sum(1))])
        err = max(err, float(np.max(np.abs(cum - (t["m"] - t["m"][0])))))
    return err


def test_misner_sharp_bookkeeping_converges(pulse_run):
    # the per-field pieces add up to the change of m along each station track (discretisation error)
    e40, e80 = _bookkeeping_error(_pulse_run(40)), _bookkeeping_error(pulse_run)
    assert e80 < 2e-5 and e80 < e40 / 3.0


def test_station_energies_and_proper_time(pulse_run):
    s = wt.summarize(pulse_run, tau_dead=1e9)
    # outbound crosses A inward then B outward; return crosses B inward then A outward
    assert s["A"]["chi1_in"]["E_kodama"] > 0 and s["B"]["chi1_out"]["E_kodama"] > 0
    assert s["B"]["chi2_in"]["E_kodama"] > 0 and s["A"]["chi2_out"]["E_kodama"] > 0
    assert s["A"]["chi1_in"]["tau50"] < s["A"]["chi2_out"]["tau50"]
    # static era: tau = t - y_station along a constant-r worldline (sigma = 0 there)
    t = pulse_run["stations"]["A"]
    y = math.sqrt(1.5**2 - 1.0)
    early = t["tau"] < 0.5
    tcoord = 0.5 * (t["u"] + t["v"])
    assert np.max(np.abs(t["tau"][early] - (tcoord[early] - y))) < 1e-3


def test_constraint_residual_converges():
    out = []
    for n in (40, 80):
        p = wt.Params(n=n, u_max=4.0, v_max=4.0, a1=0.01, v1=1.0, v2=2.0, out_h=0.1)
        cm = wt.evolve(p)["cons"]
        out.append(max(float(np.max(cm[c])) for c in ("u", "v")))
    assert out[1] < out[0] / 3.0


def _fake_run(frac, tau90, cons11, passed=True):
    return {
        "return": {"fraction": frac, "pass": passed},
        "A": {"chi2_out": {"tau90": tau90}},
        "constraint_tauA_11": {"A": {"u": {"max": cons11}, "v": {"max": 0.0}}},
        "constraint": {"A": {"u": {"max": cons11}, "v": {"max": 0.0}}},
        "bookkeeping": {"A": 0.0},
    }


def test_certify_pass_rule():
    ok = wt.certify_pass(_fake_run(0.70, 10.5, 1e-5), _fake_run(0.71, 10.5, 1e-5))
    assert ok["class"] == "certified"
    # constraint in A's causal past not small -> exploratory
    assert wt.certify_pass(_fake_run(0.70, 10.5, 1e-1), _fake_run(0.71, 10.5, 1e-5))["class"] == (
        "exploratory"
    )
    # verdict needs late energy: 0.9 * 0.55 < 0.5
    assert wt.certify_pass(_fake_run(0.55, 10.5, 1e-5), _fake_run(0.55, 10.5, 1e-5))["class"] == (
        "exploratory"
    )
    # more than 10 % of the energy arrives after tau_A = 11
    assert wt.certify_pass(_fake_run(0.70, 11.2, 1e-5), _fake_run(0.70, 11.2, 1e-5))["class"] == (
        "exploratory"
    )
    assert wt.certify_pass(_fake_run(0.40, 10.5, 1e-5, False), None)["class"] == "fail"
