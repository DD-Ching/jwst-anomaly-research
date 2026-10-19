"""Frozen rule and helpers of scripts/smacs_at2024jfz_psfmatch.py (offline)."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from smacs_at2024jfz_psfmatch import (  # noqa: E402
    control_passes,
    fit_mask,
    outcome_psfmatch,
    psf_matched_diff,
)

D = np.zeros(20)


def test_outcome_absent_present():
    assert outcome_psfmatch(D, 8.0, 100.0, 10.0, 200.0, True) == "absent"
    d = D.copy()
    d[3] = 120.0
    assert outcome_psfmatch(d, 8.0, 100.0, 10.0, 200.0, True) == "present"
    # above thr but below 5 sigma
    assert outcome_psfmatch(d, 30.0, 100.0, 10.0, 140.0, True) == "untestable"  # inj < 5 sigma
    assert outcome_psfmatch(d, 30.0, 100.0, 10.0, 400.0, True) == "absent"


def test_outcome_gates():
    assert outcome_psfmatch(D, 8.0, 100.0, 10.0, 200.0, False) == "untestable"
    assert outcome_psfmatch(D, 8.0, 100.0, -150.0, 200.0, True) == "untestable"
    assert outcome_psfmatch(D, 8.0, 100.0, np.nan, 200.0, True) == "untestable"
    assert outcome_psfmatch(D, 8.0, 100.0, 10.0, 90.0, True) == "untestable"
    assert outcome_psfmatch(D, 8.0, 100.0, 10.0, np.nan, True) == "untestable"
    d = D.copy()
    d[:5] = np.nan  # 25 % NaN
    assert outcome_psfmatch(d, 8.0, 100.0, 10.0, 200.0, True) == "untestable"
    assert outcome_psfmatch(np.array([]), 8.0, 100.0, 10.0, 200.0, True) == "untestable"


def test_control():
    assert control_passes(200.0, 8.0, 223.4)
    assert not control_passes(100.0, 8.0, 223.4)
    assert not control_passes(200.0, 50.0, 223.4)
    assert not control_passes(np.nan, 8.0, 223.4)


def test_fit_mask():
    m = fit_mask((21, 21), np.array([[10.0, 10.0]]), 3.0)
    assert not m[10, 10] and not m[10, 13] and m[10, 14] and m[0, 0]


def test_psf_matched_diff_removes_static_source():
    rng = np.random.default_rng(1)
    yy, xx = np.mgrid[0:80, 0:80]

    def gauss(x0, y0, s, f):
        g = np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * s * s))
        return f * g / g.sum()

    t = sum(gauss(x, y, 1.0, f) for x, y, f in [(20, 20, 50), (55, 30, 80), (40, 60, 30), (62, 62, 40)])
    e = sum(gauss(x, y, 1.8, f) for x, y, f in [(20, 20, 50), (55, 30, 80), (40, 60, 30), (62, 62, 40)])
    e = e + gauss(30, 45, 1.8, 20.0)  # a new source in E only
    noise = 0.01
    e = e + rng.normal(0, noise, e.shape)
    t = t + rng.normal(0, noise, t.shape)
    var = np.full(e.shape, noise**2)
    mask = fit_mask(e.shape, np.array([[30.0, 45.0]]), 6.0)
    d, k, _ = psf_matched_diff({"sci": e, "err": np.sqrt(var)}, t, var, mask)
    assert abs(k.sum() - 1.0) < 0.05
    # the static source at (55, 30) is subtracted to a few per cent; the new source stays
    assert abs(np.nansum(d[25:36, 50:61])) < 0.05 * 80
    assert np.nansum(d[40:51, 25:36]) > 0.8 * 20
