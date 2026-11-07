import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import miri_fragments as m  # noqa: E402


def test_adjacent_labels_four_connected_only():
    segm = np.array([[1, 1, 0, 3], [1, 2, 0, 0], [0, 0, 4, 0]])
    adj = m.adjacent_labels(segm)
    assert adj[1] == {2}
    assert adj[2] == {1}
    assert 3 not in adj  # touches only background
    assert 4 not in adj  # 2 and 4 touch diagonally only


def test_bright_adjacent_needs_factor():
    adj = {1: {2}, 2: {1, 3}, 3: {2}}
    flag = m.bright_adjacent([1, 2, 3], [100.0, 10.0, 4.0], adj, factor=5.0)
    assert flag.tolist() == [False, True, False]  # 3 has neighbour 2 at only 2.5x


def test_mechanism_requires_all_tested_tertiles():
    rng = np.random.default_rng(0)
    n = 300
    flux = rng.uniform(1, 100, n)
    flag = rng.random(n) < 0.4
    ci = rng.normal(1.0, 0.05, n) + 0.3 * flag
    assert m.mechanism_test(ci, flux, flag)["supported"]
    assert not m.mechanism_test(rng.normal(1.0, 0.05, n), flux, flag)["supported"]
    # an effect only in the faintest tertile does not pass the "all tertiles" rule
    ci2 = rng.normal(1.0, 0.05, n) + 0.3 * flag * (flux < np.quantile(flux, 1 / 3))
    assert not m.mechanism_test(ci2, flux, flag)["supported"]
