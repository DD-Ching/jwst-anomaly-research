from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest
from astropy.cosmology import Planck18

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import w1_volume as wv  # noqa: E402


def test_constant_efficiency_is_the_shell_volume():
    area = 1.0  # deg^2
    v = wv.v_eff_mpc3([0.1, 1.5], [[0.5, 0.5]], [area])
    shell = (Planck18.comoving_volume(1.5) - Planck18.comoving_volume(0.1)).to_value("Mpc3")
    assert v == pytest.approx(0.5 * shell * area / 41252.96, rel=1e-4)


def test_fields_add_and_zero_efficiency_gives_no_limit():
    one = wv.v_eff_mpc3([0.1, 0.5, 1.0], [[0.1, 0.2, 0.0]], [2.0])
    two = wv.v_eff_mpc3([0.1, 0.5, 1.0], [[0.1, 0.2, 0.0], [0.1, 0.2, 0.0]], [1.0, 1.0])
    assert two == pytest.approx(one)
    assert wv.limit(0.0)["n95_mpc3"] is None


def test_limit_units():
    lim = wv.limit(1e5, h=0.5)
    assert lim["n95_mpc3"] == pytest.approx(-math.log(0.05) / 1e5)
    assert lim["n95_h3_mpc3"] == pytest.approx(lim["n95_mpc3"] / 0.125)
    assert lim["ratio_to_ta13"] == pytest.approx(lim["n95_h3_mpc3"] / 1e-4)


def test_grid_must_increase():
    with pytest.raises(ValueError):
        wv.v_eff_mpc3([0.5, 0.1], [[1, 1]], [1.0])
