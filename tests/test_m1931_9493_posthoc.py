import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from m1931_9493_posthoc import bright_rows, nearest_bright  # noqa: E402


def test_nearest_bright_distance_and_index():
    ra = np.array([10.0, 10.001])
    dec = np.array([0.0, 0.0])
    k, d = nearest_bright(10.0011, 0.0, ra, dec)
    assert k == 1
    assert abs(d - 0.36) < 1e-6


def test_nearest_bright_empty():
    assert nearest_bright(1.0, 1.0, np.array([]), np.array([])) is None


def test_bright_rows_percentile_ignores_nan():
    t = Table(
        {"F444W_flux1": [1.0, 2.0, np.nan, 100.0, 3.0], "ra": np.arange(5.0), "dec": np.zeros(5)}
    )
    b = bright_rows(t, pct=75)
    assert list(b["F444W_flux1"]) == [100.0]


def test_bright_rows_all_nan():
    t = Table({"F444W_flux1": [np.nan, np.nan], "ra": [0.0, 1.0], "dec": [0.0, 0.0]})
    assert len(bright_rows(t)) == 0
