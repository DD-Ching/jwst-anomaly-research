"""Offline tests for scripts/a2744_uncover.py (the WLS family source; no network)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import a2744_uncover as au  # noqa: E402


def test_wls_source_equal_weights_is_mean():
    beta = np.array([[1.0, 2.0], [3.0, -2.0], [-1.0, 0.0]])
    jac = np.broadcast_to(np.eye(2), (3, 2, 2))
    assert np.allclose(au.wls_source(beta, jac), beta.mean(0))


def test_wls_source_follows_the_most_magnified_image():
    # image 0 magnified 10x isotropically (A = 0.1 I): its source dominates by 100:1 in weight
    beta = np.array([[0.0, 0.0], [1.0, 1.0]])
    jac = np.array([0.1 * np.eye(2), np.eye(2)])
    got = au.wls_source(beta, jac)
    assert np.allclose(got, [1.0 / 101.0, 1.0 / 101.0])


def test_wls_source_minimises_image_plane_offsets():
    rng = np.random.default_rng(1)
    beta = rng.normal(size=(4, 2))
    jac = rng.normal(size=(4, 2, 2)) + 2 * np.eye(2)
    best = au.wls_source(beta, jac)
    minv = np.linalg.inv(jac)

    def cost(b):
        return sum(np.sum((minv[i] @ (b - beta[i])) ** 2) for i in range(4))

    for d in rng.normal(scale=0.05, size=(20, 2)):
        assert cost(best) <= cost(best + d)
