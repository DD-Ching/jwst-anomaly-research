"""Synthetic checks of scripts/parity_test.py: handedness recovered, symmetric source untestable."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

_P = Path(__file__).resolve().parents[1] / "scripts" / "parity_test.py"
spec = importlib.util.spec_from_file_location("parity_test", _P)
pt = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pt)

PIX = 0.031


def _grid(half):
    n = int(np.ceil(half / PIX))
    jj, ii = np.mgrid[-n : n + 1, -n : n + 1].astype(float)
    return -ii * PIX, jj * PIX  # x West, as the stamper


def _galaxy(x, y, asym=True):
    """A disc plus an off-axis knot and a bar: no mirror symmetry when ``asym``."""
    img = 50.0 * np.exp(-0.5 * (x**2 + (y / 0.6) ** 2) / 0.12**2)
    if asym:
        img += 30.0 * np.exp(-0.5 * ((x - 0.15) ** 2 + (y - 0.10) ** 2) / 0.04**2)
        img += 20.0 * np.exp(-0.5 * ((x + 0.05) ** 2 + (y + 0.18) ** 2) / 0.05**2)
    return img


def _pair(m_true, asym=True, noise=1.0, seed=1):
    rng = np.random.default_rng(seed)
    ax, ay = _grid(pt.TEMPLATE_HALF_ARCSEC)
    a = _galaxy(ax, ay, asym) + rng.normal(0, noise, ax.shape)
    tmpl = pt.Template(a, np.full(a.shape, noise), ax, ay)
    bx, by = _grid(pt.FIT_RADIUS_ARCSEC + 0.1)
    d = m_true @ np.vstack([bx.ravel(), by.ravel()])  # B offsets -> A offsets
    b = _galaxy(d[0], d[1], asym).reshape(bx.shape) + rng.normal(0, noise, bx.shape)
    return tmpl, *pt._fit_pixels(b, np.full(b.shape, noise), bx, by)


def test_decompose_roundtrip():
    for m in (np.array([[1.3, 0.4], [-0.2, 0.7]]), np.array([[0.5, 1.1], [0.9, -0.3]])):
        t1, t2, l1, l2 = pt.decompose(m)
        h = 1 if np.linalg.det(m) > 0 else -1
        got, _ = pt.affine(np.array([t1, t2, l1, l2, 0, 0]), h)
        assert np.allclose(got, m, atol=1e-9)


@pytest.mark.parametrize("mirror", [False, True])
def test_handedness_recovered(mirror):
    m_true = pt._rot(0.4) @ np.diag([1.4, -0.8 if mirror else 0.8]) @ pt._rot(-0.9)
    tmpl, bx, by, b, bvar = _pair(m_true)
    # the model predicts the true map: consistent, and the synthetics show power
    res = pt.pair_test(tmpl, bx, by, b, bvar, m_true)
    assert res["class"] == "consistent"
    assert res["h_pred"] == (-1 if mirror else 1)
    # a model predicting the opposite handedness: flagged as a parity residual
    wrong = m_true @ np.diag([1.0, -1.0])
    assert pt.pair_test(tmpl, bx, by, b, bvar, wrong)["class"] == "parity_residual"


def test_symmetric_source_untestable():
    m_true = np.diag([1.2, -0.9])  # mirror across an axis of the symmetric disc
    tmpl, bx, by, b, bvar = _pair(m_true, asym=False)
    assert pt.pair_test(tmpl, bx, by, b, bvar, m_true)["class"] == "untestable"


def test_segment_separates_neighbour_and_size_ratio():
    x, y = _grid(pt.TEMPLATE_HALF_ARCSEC)
    rng = np.random.default_rng(3)
    src = _galaxy(x, y)
    neighbour = 80.0 * np.exp(-0.5 * ((x - 0.9) ** 2 + (y + 0.4) ** 2) / 0.08**2)
    sci = src + neighbour + rng.normal(0, 1.0, x.shape)
    own, other = pt.segment(sci, np.ones_like(sci), x, y)
    assert own[np.hypot(x, y) < 0.05].all() and not own[np.hypot(x - 0.9, y + 0.4) < 0.05].any()
    assert other[np.hypot(x - 0.9, y + 0.4) < 0.05].all()
    assert not (own & other).any()
    # the extended galaxy is resolved; a PSF-sized point source is not
    assert pt.size_ratio(src, x, y, own) > pt.MIN_SIZE_RATIO
    s = pt.PSF_SIGMA_ARCSEC
    point = 200.0 * np.exp(-0.5 * (x**2 + y**2) / s**2)
    own_p, _ = pt.segment(point + rng.normal(0, 1.0, x.shape), np.ones_like(x), x, y)
    assert pt.size_ratio(point, x, y, own_p) < pt.MIN_SIZE_RATIO
