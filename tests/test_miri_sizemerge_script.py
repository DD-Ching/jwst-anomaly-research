import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import miri_sizemerge as m  # noqa: E402


def test_adopt_alpha_takes_largest_passing():
    assert m.adopt_alpha([0.0, 0.25, 0.5, 0.75], [None, 0.04, 0.115, 0.08]) == 0.75
    assert m.adopt_alpha([0.0, 0.25, 0.5], [None, 0.11, 0.2]) == 0.0
    assert m.adopt_alpha([0.0, 0.25], [None, 0.10]) == 0.25  # ceiling is inclusive
    assert m.adopt_alpha([0.0, 0.25], [None, None]) == 0.0


def test_chance_fraction():
    assert m.chance_fraction(31, 1.25) == 1.25 / 31
    assert m.chance_fraction(0, 3.0) is None


def test_shifted_moves_on_sky():
    t = Table({"ra": [10.0], "dec": [60.0]})
    s = m.shifted(t, 3.6, -3.6)
    assert np.isclose((s["dec"][0] - 60.0) * 3600, -3.6)
    assert np.isclose((s["ra"][0] - 10.0) * 3600 * np.cos(np.deg2rad(60.0)), 3.6)
    assert t["ra"][0] == 10.0  # input unchanged


def test_n_pairs():
    assert m.n_pairs(Table({"n_bands": [1, 2, 4]})) == 4


def test_chance_fraction_net_loss_never_passes():
    assert m.chance_fraction(10, -0.5) is None
    assert m.adopt_alpha([0.0, 0.25], [None, m.chance_fraction(10, -0.5)]) == 0.0
