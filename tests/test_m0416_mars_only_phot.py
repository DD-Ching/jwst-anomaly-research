import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import m0416_mars_extra as me  # noqa: E402
import m0416_mars_only_phot as mp  # noqa: E402


def test_match_pairs_agrees_with_match_observed():
    px, py, ox, oy = [0.0, 1.0, 10.0], [0, 0, 0], [0.2, 0.6], [0, 0]
    pair = mp.match_pairs(px, py, ox, oy)
    assert pair.tolist() == [0, 1]
    used = np.zeros(3, bool)
    used[pair[pair >= 0]] = True
    assert used.tolist() == me.match_observed(px, py, ox, oy).tolist()


def test_match_pairs_radius_and_empty():
    assert mp.match_pairs([5.0], [0.0], [0.0], [0.0]).tolist() == [-1]
    assert mp.match_pairs([], [], [0.0], [0.0]).tolist() == [-1]
    assert mp.match_pairs([1.0], [0.0], [], []).tolist() == []


def test_mars_only():
    # image 0 paired; image 1 is 2" from an extra prediction; image 2 is 3.5" from it
    m = mp.mars_only([0.0, 10.0, 11.5], [0, 0, 0], [0], [8.0], [0.0])
    assert m.tolist() == [False, False, True]
    assert mp.mars_only([0.0, 5.0], [0, 0], [-1], [], []).tolist() == [True, True]
    # exactly 3.0" from an extra prediction is not MARS-only (rule 2: farther than 3")
    assert mp.mars_only([3.0], [0.0], [], [0.0], [0.0]).tolist() == [False]


def test_reading():
    assert mp.reading("absent", True) == "lensing_residual_input"
    assert mp.reading("absent", False) == "no_reading"
    assert mp.reading("recovered", False) == "seen_where_only_mars_predicts"
    for c in ("confused", "ambiguous", "undetectable", "undetectable_local", "no_reference"):
        assert mp.reading(c, True) == "no_reading"
