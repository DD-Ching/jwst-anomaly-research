import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import w5_radial_counts as w5  # noqa: E402


def test_decide():
    assert w5.decide(0.005) == "vet the centres"
    assert w5.decide(0.02) == "exploratory note"
    assert w5.decide(0.05) == "null"


def test_mean_ratio_deficit_deep_inside_theta_e_only():
    rng = np.random.default_rng(1)
    mags = 22 + 5 * rng.random(5000) ** 0.5  # rising counts towards the faint end
    m_lim = float(np.max(mags))
    inside = w5.mean_ratio(1.0, 3.0, mags, m_lim)  # r <= theta_E / 3: the dark hole
    far = w5.mean_ratio(3.0, 0.3, mags, m_lim)
    assert inside < 0.9
    assert abs(far - 1.0) < 0.1


def test_combine_power_and_p():
    f = {"radii": {str(r): {"n_obs": 5, "ctrl_sum": list(np.arange(1, 21)),
                            "q_hypothesis": {f"{m:.0e}": 0.0 for m in w5.MASSES}}
                   for r in w5.RADII}}  # fmt: skip
    out = w5.combine([f])
    one = out["radii"]["3.0"]
    assert one["p_low"] == 0.25 and one["p_high"] == 0.8
    assert all(v == 1.0 for v in one["power"].values())  # Q_H = 0 always below the 5th percentile


def test_control_positions_annulus_and_exclusion():
    rng = np.random.default_rng(2)
    pts = rng.uniform(-100, 100, (4000, 2))
    centres = np.array([[30.0, 0.0]])
    pos = w5.control_positions(30.0, centres, cKDTree(pts), 6.0, 300, rng)
    r = np.hypot(pos[:, 0], pos[:, 1])
    assert np.all((r >= 20) & (r <= 40))
    assert np.all(np.hypot(pos[:, 0] - 30, pos[:, 1]) >= 12)
