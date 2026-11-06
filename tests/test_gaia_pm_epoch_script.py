import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import gaia_pm_epoch as g  # noqa: E402


def test_propagate_moves_by_pm_times_dt():
    ra, dec = g.propagate(
        np.array([10.0]), np.array([60.0]), np.array([100.0]), np.array([-50.0]), 2.0
    )
    # mu_alpha* includes cos(dec): 200 mas east at dec 60 is 400 mas of RA
    assert np.isclose((ra[0] - 10.0) * 3.6e6, 400.0)
    assert np.isclose((dec[0] - 60.0) * 3.6e6, -100.0)


def test_label_at_radii():
    assert g.label_at(0.3, False) == "position"
    assert g.label_at(0.31, True) == "astrometry_only"
    assert g.label_at(0.31, False) == "missed"
    assert g.label_at(0.5, True) == "astrometry_only"
    assert g.label_at(0.51, True) == "missed"
