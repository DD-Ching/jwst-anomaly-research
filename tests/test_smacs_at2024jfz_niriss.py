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


def test_posthoc_outcome_and_point_disc():
    assert m.posthoc_outcome(np.array([5.0, 120.0]), 10.0, 0.0, 100.0) == "present"
    assert m.posthoc_outcome(np.array([5.0, 90.0]), 10.0, 0.0, 100.0) == "absent"
    assert m.posthoc_outcome(np.array([5.0]), 10.0, 0.3, 100.0) == "untestable"
    disc = m.grid_points(np.array([[110.83, -73.45], [110.83, -73.45]]))
    assert 70 < len(disc) < 90  # pi 0.5^2 / 0.01


def test_posthoc_injection_gate():
    assert m.posthoc_outcome(np.array([5.0]), 10.0, 0.0, 100.0, inj_delta=180.0) == "absent"
    assert m.posthoc_outcome(np.array([5.0]), 10.0, 0.0, 100.0, inj_delta=90.0) == "untestable"
    assert m.posthoc_outcome(np.array([5.0]), 30.0, 0.0, 100.0, inj_delta=120.0) == "untestable"
    assert m.posthoc_outcome(np.array([5.0]), 10.0, 0.0, 100.0, inj_delta=np.nan) == "untestable"


def test_classify_flag():
    wing = {
        "flag_to_source_arcsec": 0.2,
        "e_minus_t_offset_arcsec": 0.012,
        "flux_e_minus_t_njy": 207.0,
    }
    assert m.classify_flag("present", wing, 97.0) == "undetermined"
    assert m.classify_flag("present", wing | {"flux_e_minus_t_njy": -50.0}, 97.0) == "absent"
    assert m.classify_flag("present", wing | {"flag_to_source_arcsec": 0.45}, 97.0) == "present"
    assert m.classify_flag("present", wing | {"e_minus_t_offset_arcsec": 0.2}, 97.0) == "present"
    assert m.classify_flag("present", None, 97.0) == "present"
    assert m.classify_flag("absent", wing, 97.0) == "absent"
