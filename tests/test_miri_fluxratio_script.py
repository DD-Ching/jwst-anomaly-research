import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import miri_fluxratio as m  # noqa: E402


def test_leave_one_out_uses_other_usable_images_only():
    mag = np.array([25.0, 25.2, 26.0, 30.0, 24.0])
    system = np.array(["1", "1", "1", "1", "2"])
    usable = np.array([True, True, True, False, True])
    res = m.leave_one_out(mag, system, usable)
    assert np.isclose(res[0], 25.0 - 25.6)  # median of 25.2 and 26.0
    assert np.isclose(res[2], 26.0 - 25.1)
    assert np.isnan(res[3])  # unusable image is not tested
    assert np.isnan(res[4])  # a single usable image has no comparison


def test_pair_carries_full_difference():
    res = m.leave_one_out(np.array([25.0, 26.0]), np.array(["3", "3"]), np.array([True, True]))
    assert np.allclose(res, [-1.0, 1.0])
