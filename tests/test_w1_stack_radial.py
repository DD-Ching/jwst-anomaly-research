import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import w1_stack_radial as wsr  # noqa: E402


def test_decide():
    assert wsr.decide(0.005, True) == "vet the centres"
    assert wsr.decide(0.02, True) == "exploratory note"
    assert wsr.decide(0.05, True) == "null"
    assert wsr.decide(0.001, False).startswith("no result")


def test_combine_pools_null_by_centre_count():
    f1 = {
        "centres": [{"s": 1.0, "s_cross": 0.0}],
        "null_sum": np.array([0.0, 2.0]),
        "null_n": np.array([1, 1]),
    }
    f2 = {
        "centres": [{"s": 0.0, "s_cross": 0.0}, {"s": None, "s_cross": None}],
        "null_sum": np.array([0.0, 0.0]),
        "null_n": np.array([1, 1]),
    }
    out = wsr.combine([f1, f2])
    assert out["n_centres"] == 2 and out["mean_s"] == 0.5
    assert out["p_one_sided"] == 0.5  # null means 0.0 and 1.0


def test_radial_stack_recovers_injected_radial_shear():
    es = wsr.es
    rng = np.random.default_rng(0)
    x, y = rng.uniform(-40, 40, (2, 4000))
    e = 0.2 * np.exp(1j * rng.uniform(0, 2 * np.pi, 4000))
    phi = np.arctan2(-x, y)  # position angle seen from the origin (x West)
    e_rad = e + 0.1 * np.exp(2j * phi)  # PA-based: radial alignment (negative tangential shear)
    ap = es.ApertureMass(x, y, [0.0], [0.0], radius=10.0, kind="schirmer", r_min=1.0)
    s_rad, _ = ap.snr(e_rad)
    s_tan, _ = ap.snr(e - 0.1 * np.exp(2j * phi))
    assert s_rad[0] > 3 and s_tan[0] < -3
    tot, n = wsr.stack_null(ap, e_rad, 200, rng)
    assert np.mean(tot / n >= s_rad[0]) < 0.02
