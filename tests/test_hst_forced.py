from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import hst_forced as hf  # noqa: E402


def test_centroid_and_peaks():
    yy, xx = np.mgrid[0:60, 0:60].astype(float)
    img = np.random.default_rng(0).normal(0, 1, (60, 60))
    img += 200 * np.exp(-((xx - 20.3) ** 2 + (yy - 30.0) ** 2) / 2)
    img += 100 * np.exp(-((xx - 45.0) ** 2 + (yy - 12.6) ** 2) / 2)
    pk = hf.find_peaks(img, 0.05, 20.0, 5)
    assert len(pk) == 2
    assert abs(pk[0][0] - 20.3) < 0.2 and abs(pk[0][1] - 30.0) < 0.2  # brightest first
    assert abs(pk[1][1] - 12.6) < 0.2
    assert hf.centroid(img, 0, 0) == (0.0, 0.0)  # edge: no 3x3 box, returns the peak
