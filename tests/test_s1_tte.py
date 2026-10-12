"""Offline tests for scripts/s1_tte.py (synthetic events; no network)."""

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))
_spec = importlib.util.spec_from_file_location("s1_tte", ROOT / "scripts" / "s1_tte.py")
st = importlib.util.module_from_spec(_spec)
sys.modules["s1_tte"] = st
_spec.loader.exec_module(st)


def test_common_gti_intersection():
    evs = {"n0": (None, None, [(-30.0, 200.0)]), "n1": (None, None, [(-26.0, 300.0)])}
    assert st.common_gti(evs) == [(-26.0, 200.0)]


def test_background_ignores_no_data_interval():
    rng = np.random.default_rng(1)
    # 200 counts/s background but data only from t = -26 s (TTE start); the alt window reaches -40 s
    t = np.sort(rng.uniform(-26.0, 80.0, int(200 * 106)))
    bands = [t, t]
    gti = [(-26.0, 80.0)]
    f_gti, _ = st.curve(bands, 0.0, 0.5, 0.016, 0.0, bkg="alt", gti=gti)
    f_raw, _ = st.curve(bands, 0.0, 0.5, 0.016, 0.0, bkg="alt", gti=None)
    # with the GTI the net rate is ~0; without it 14 of the 85 background bins are empty, which
    # biases the background low by ~200 * 14 / 85 = 33 counts/s and leaves a fake positive net rate
    assert abs(np.mean(f_gti)) < 15 and np.mean(f_raw) > 25


def test_single_pulse_simulation_is_not_counted_as_two():
    rng = np.random.default_rng(3)
    ev = st.sim_single_pulse(rng, np.array([3000.0, 2000.0]), np.array([300.0, 300.0]), 0.0, 0.5)
    assert st.npulse_multires(ev, 0.0, 0.5) <= 1


def test_epeak_z_sign_convention():
    # b is a's copy stretched in time: y = b stretched by s matches a, so E_b = E_a * s
    lna = math.log(300.0)
    s = 1.5
    assert abs(st.zstat(lna, 0.1, math.log(300.0 * s), 0.1, s)) < 1e-9
    assert abs(st.zstat(lna, 0.1, math.log(300.0 / s), 0.1, s)) > 2
