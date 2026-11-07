import sys
from pathlib import Path

import numpy as np
from astropy.wcs import WCS

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import miri_segoverlap as m  # noqa: E402


def test_adopt_f_takes_smallest_passing():
    assert m.adopt_f([0.1, 0.25, 0.5, 0.75], [0.3, 0.09, 0.05, None]) == 0.25
    assert m.adopt_f([0.1, 0.25], [0.10, 0.2]) == 0.1  # ceiling is inclusive
    assert m.adopt_f([0.1, 0.25], [0.2, None]) is None


def test_chance_fraction():
    assert m.chance_fraction(20, 1.0) == 0.05
    assert m.chance_fraction(0, 1.0) is None


def test_overlap_pairs_uses_smaller_area():
    a = np.zeros((4, 6), int)
    b = np.zeros((4, 6), int)
    a[0:2, 0:4] = 1  # 8 px
    b[0:2, 2:3] = 7  # 2 px, fully inside label 1
    b[3, 5] = 8  # no overlap
    i, j, o, n = m.overlap_pairs(a, b)
    assert list(zip(i.tolist(), j.tolist(), o.tolist(), strict=True)) == [(1, 7, 1.0)]
    assert n.tolist() == [2]
    # native band-b area larger than its visible part lowers the overlap
    _, _, o2, _ = m.overlap_pairs(a, b, np.array([0, 0, 0, 0, 0, 0, 0, 4, 1]))
    assert o2.tolist() == [0.5]


def test_greedy_one_to_one():
    i = np.array([1, 1, 2])
    j = np.array([5, 6, 5])
    s = np.array([0.9, 0.6, 0.8])
    assert m.greedy_one_to_one(i, j, s, 0.5) == [(1, 5)]
    assert m.greedy_one_to_one(i, j, s, 0.95) == []


def _wcs():
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [110.0, -73.0]
    w.wcs.crpix = [5.0, 5.0]
    w.wcs.cdelt = [-0.1 / 3600, 0.1 / 3600]
    return w


def test_resample_identity_and_shift():
    w = _wcs()
    seg = np.zeros((10, 10), int)
    seg[5, 5] = 3
    ra, dec = m.ref_sky(w, seg.shape)
    assert np.array_equal(m.resample_labels(seg, w, ra, dec), seg)
    # Moving the map 0.2" north puts the label two pixels up (dec increases with y here).
    sh = m.resample_labels(seg, w, ra, dec, (0.0, 0.2))
    assert sh[7, 5] == 3 and sh[5, 5] == 0


def test_band_shift_rotates_keeping_length():
    for k in range(3):
        dx, dy = m.band_shift(k, 8.0, 8.0)
        assert np.isclose(np.hypot(dx, dy), np.hypot(8.0, 8.0))
    assert np.allclose(m.band_shift(0, 8.0, -8.0), (8.0, -8.0))


def test_greedy_tiebreak_by_shared_pixels():
    i = np.array([1, 1])
    j = np.array([5, 9])
    s = np.array([1.0, 1.0])
    assert m.greedy_one_to_one(i, j, s, 0.5) == [(1, 5)]  # label order
    assert m.greedy_one_to_one(i, j, s, 0.5, np.array([3, 40])) == [(1, 9)]
