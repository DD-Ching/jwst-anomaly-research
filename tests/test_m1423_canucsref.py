"""Unit tests of the MACS1423 local-noise rule (plan rule 2)."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import m1423_canucsref as m  # noqa: E402


def test_local_class_thresholds():
    assert m.local_class("recovered", 4.9, 1.0) == "undetectable_local"
    assert m.local_class("recovered", 5.0, 1.0) == "recovered"
    assert m.local_class("absent", 9.9, 1.0) == "ambiguous"
    assert m.local_class("absent", 10.0, 1.0) == "absent"
    assert m.local_class("ambiguous", 6.0, 1.0) == "ambiguous"
    assert m.local_class("undetectable", 100.0, 1.0) == "undetectable"
    assert m.local_class("no_reference", float("nan"), float("nan")) == "no_reference"
    assert m.local_class("confused", 50.0, float("nan")) == "undetectable_local"


def test_local_noise_matches_white_noise():
    rng = np.random.default_rng(1)
    d = rng.normal(0.0, 1.0, (400, 400))
    npix = np.count_nonzero((np.mgrid[-9:10, -9:10] ** 2).sum(axis=0) <= (m.AP_ARCSEC / m.PIX) ** 2)
    s = m.local_noise(d, 200.0, 200.0, np.array([]), np.array([]), np.random.default_rng(2))
    assert s == pytest.approx(np.sqrt(npix), rel=0.2)


def test_local_noise_raises_when_all_excluded():
    d = np.zeros((400, 400))
    xs, ys = np.meshgrid(np.arange(0, 400, 5.0), np.arange(0, 400, 5.0))
    with pytest.raises(RuntimeError):
        m.local_noise(d, 200.0, 200.0, xs.ravel(), ys.ravel(), np.random.default_rng(0))
