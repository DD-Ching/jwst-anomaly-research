import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import smacs_at2024jfz_niriss as m  # noqa: E402

INJ_OK = np.full(m.N_TRACK, 400.0)


def test_outcome_absent_present_untestable():
    quiet = np.array([10.0, -20.0, np.nan, 50.0])
    assert m.outcome(quiet, 20.0, 0.0, 0.0, INJ_OK, 20.0) == "absent"
    # needs both the flux threshold and 5 sigma
    assert m.outcome(np.array([300.0]), 20.0, 0.0, 0.0, INJ_OK, 20.0) == "present"
    assert m.outcome(np.array([300.0]), 70.0, 0.0, 0.0, INJ_OK, 20.0) == "absent"
    assert m.outcome(np.array([150.0]), 5.0, 0.0, 0.0, INJ_OK, 20.0) == "absent"
    # coverage and injection gates
    assert m.outcome(quiet, 20.0, 0.25, 0.0, INJ_OK, 20.0) == "untestable"
    assert m.outcome(quiet, 20.0, 0.0, 0.21, INJ_OK, 20.0) == "untestable"
    inj = INJ_OK.copy()
    inj[:2] = 100.0
    assert m.outcome(quiet, 20.0, 0.0, 0.0, inj, 20.0) == "untestable"
    inj[1] = 400.0
    assert m.outcome(quiet, 20.0, 0.0, 0.0, inj, 20.0) == "absent"


def test_track_and_grid():
    track = m.track_points((110.835165, -73.453608), (110.835136, -73.453496))
    assert track.shape == (m.N_TRACK, 2)
    grid = m.grid_points(track)
    cosd = np.cos(np.radians(track[:, 1].mean()))
    # every grid point within the search radius of some track point (track spacing < grid step)
    d = np.hypot(
        (grid[:, None, 0] - track[None, :, 0]) * cosd * 3600,
        (grid[:, None, 1] - track[None, :, 1]) * 3600,
    ).min(axis=1)
    assert d.max() <= m.SEARCH_ARCSEC + 0.05
    # area of a 0.5" stadium around a ~0.41" segment over 0.01 arcsec^2 per point
    assert 100 < len(grid) < 140


def test_inject_gaussian_integral():
    sci = np.zeros((41, 41))
    out = m.inject_gaussian(sci, 20.0, 20.0, 100.0, 3.0)
    assert abs(out.sum() - 100.0) < 1e-9
    assert np.unravel_index(out.argmax(), out.shape) == (20, 20)
