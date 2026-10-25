import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import sneos_literature as sl  # noqa: E402


def test_ratio_from_magnitudes():
    r, s = sl.ratio(25.0, 0.0, 25.0 - 2.5 * np.log10(2.0), 0.1)
    assert np.isclose(r, 2.0)
    assert np.isclose(s, 0.4 * np.log(10) * 0.1)


def test_weighted_mean_and_chi2():
    w = sl.weighted(np.log([1.0, 4.0]), np.array([0.1, 0.1]))
    assert np.isclose(w["r"], 2.0)
    assert np.isclose(w["sigma_ln"], 0.1 / np.sqrt(2))
    assert w["dof"] == 1 and w["chi2"] > 0


def test_difference_ratio_equals_mu_ratio_for_constant_ratio():
    # 101.2 always 1.5x brighter: the ratio of the flux drops is 1.5 too.
    m = lambda f: -2.5 * np.log10(f)  # noqa: E731
    early = (0, 0, "b", m(10.0), 0.01, m(15.0), 0.01)
    late = (0, 0, "b", m(4.0), 0.01, m(6.0), 0.01)
    r, _ = sl.difference_ratio(early, late)
    assert np.isclose(r, 1.5)


def test_run_writes_summary(tmp_path):
    s = sl.run(tmp_path)
    assert (tmp_path / "literature.json").exists()
    assert set(s["combined"]) == {"HST F814W", "HST F110W", "JWST 6882", "JWST 9493"}
    lo, hi = s["r_model_range"]
    assert lo <= s["r_model"]["pinned map"] <= hi


def test_scaled_sigma_grows_with_excess_scatter():
    w = sl.weighted(np.log([1.0, 4.0]), np.array([0.01, 0.01]))
    assert np.isclose(w["sigma_ln_scaled"], w["sigma_ln"] * np.sqrt(w["chi2"]))


def test_slope():
    b, _ = sl.slope_mag_per_day(np.arange(4.0), 25 + 0.1 * np.arange(4.0), np.full(4, 0.05))
    assert np.isclose(b, 0.1)
