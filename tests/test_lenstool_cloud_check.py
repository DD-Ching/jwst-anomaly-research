import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import lenstool_cloud_check as m  # noqa: E402

FLAGS_OK = dict.fromkeys(m.FLAG_IMAGES, True)


def test_read_image_all(tmp_path):
    p = tmp_path / "image.all"
    p.write_text(
        "#REFERENCE 3 64.0382170 -24.0675010\n"
        "A.1 -1.0 2.0 0.3 0.1 10 1.0 25\n"
        "A.1 3.0 4.0 0.3 0.1 10 1.0 25\n"
        "B.2 5.0 6.0 0.3 0.1 10 2.0 25\n"
    )
    ref, rows = m.read_image_all(p)
    assert ref == (64.038217, -24.067501)
    assert rows["A.1"].shape == (2, 2) and rows["B.2"].tolist() == [[5.0, 6.0]]


def test_match_pairs_one_to_one_and_unmatched():
    ref = np.array([[0.0, 0.0], [1.0, 0.0], [5.0, 5.0]])
    port = np.array([[1.02, 0.0], [0.0, 0.05], [5.0, 5.5]])
    d = m.match_pairs(ref, port)
    assert np.allclose(d[:2], [0.05, 0.02]) and np.isnan(d[2])  # 0.5" pair is not a match
    assert np.isnan(m.match_pairs(ref, np.empty((0, 2)))).all()


def test_gate_r_needs_fraction_and_both_flags():
    good = np.r_[np.full(95, 0.05), np.full(5, np.nan)]
    assert m.gate_r(good, FLAGS_OK, FLAGS_OK)
    bad = np.r_[np.full(94, 0.05), np.full(6, 0.2)]
    assert not m.gate_r(bad, FLAGS_OK, FLAGS_OK)
    assert not m.gate_r(good, FLAGS_OK, {"D115b.1": True, "D115b.2": False})
    assert not m.gate_r(good, {}, FLAGS_OK)
    assert not m.gate_r(np.array([]), FLAGS_OK, FLAGS_OK)


def test_match_pairs_maximises_pairs_within_threshold():
    # both pairs lie within the threshold and are kept one-to-one
    ref = np.array([[0.0, 0.0], [0.19, 0.0]])
    port = np.array([[0.09, 0.0], [0.29, 0.0]])
    d = m.match_pairs(ref, port)
    assert np.allclose(d, [0.09, 0.10])


def test_near_flag():
    assert m.near_flag(np.array([[1.0, 1.0], [10.0, 10.0]]), 10.5, 10.2)
    assert not m.near_flag(np.array([[1.0, 1.0]]), 10.5, 10.2)
    assert not m.near_flag(np.empty((0, 2)), 0.0, 0.0)


def test_gate_r_fails_on_missing_images():
    good = np.full(100, 0.05)
    assert not m.gate_r(good, FLAGS_OK, FLAGS_OK, n_missing=1)
