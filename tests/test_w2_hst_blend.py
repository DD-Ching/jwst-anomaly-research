"""Offline tests of the line-C HST rules (scripts/w2_hst_blend_chain.py); synthetic data only."""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from astropy.modeling import models
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import w2_hst_blend as cli  # noqa: E402
import w2_hst_blend_chain as ch  # noqa: E402
import w12_niq_hst as hst  # noqa: E402

P = ch.Params()
PIX = 0.04


def _classes():
    rows = [
        ("L1", "blended", True),
        ("L2", "too close", False),
        ("L3", "blended", False),
        ("L4", "too close", True),
        ("L5", "none", False),  # decided by LS: never selected
        ("L6", "blended", False),
    ]
    return Table(
        {
            "system_id": [r[0] for r in rows],
            "name": [f"N{r[0]}" for r in rows],
            "ra": np.arange(len(rows), dtype=float),
            "dec": np.zeros(len(rows)),
            "selection": ["quasar"] * len(rows),
            "status": [r[1] for r in rows],
            "ls_sep": np.full(len(rows), np.nan),
            "theta_e": np.full(len(rows), np.nan),
            "z_lens": np.full(len(rows), np.nan),
            "z_source": np.full(len(rows), np.nan),
            "lens_z_known": [r[2] for r in rows],
        }
    )


def _coverage():
    # (system, instrument, filter, exptime, programme, contains)
    rows = [
        ("L1", "ACS/WFC", "F814W", 2000.0, "100", True),
        ("L2", "WFC3/UVIS", "F814W", 1500.0, "100", True),
        ("L3", "ACS/WFC", "F814W", 3000.0, "200", True),
        ("L4", "ACS/WFC", "F814W", 5000.0, "300", True),  # control from another programme
        ("L5", "ACS/WFC", "F814W", 9000.0, "100", True),
        ("L6", "ACS/WFC", "F814W", 900.0, "100", True),  # too shallow
        ("L6", "WFC3/IR", "F814W", 5000.0, "100", True),  # wrong instrument
        ("L2", "ACS/WFC", "F814W", 9000.0, "100", False),  # does not contain the position
    ]
    return Table(
        {
            "system_id": [r[0] for r in rows],
            "instrument": [r[1] for r in rows],
            "filters": [r[2] for r in rows],
            "t_exptime": [r[3] for r in rows],
            "proposal_id": [r[4] for r in rows],
            "contains": [r[5] for r in rows],
        }
    )


def test_select_systems_rule():
    s = ch.select_systems(_classes(), _coverage(), P)
    tg = s[s["role"] == "target"]
    ct = s[s["role"] == "control"]
    assert list(tg["system_id"]) == ["L3", "L2"]  # deepest first; L6 too shallow; L5 decided
    assert list(tg["t814_s"]) == [3000.0, 1500.0]  # non-containing exposure not counted
    assert list(ct["system_id"]) == ["L1"]  # L4 shares no programme with a target


def test_select_systems_caps_and_ties():
    p = ch.Params(max_targets=1)
    c = _classes()
    cov = _coverage()
    cov["t_exptime"][1] = 3000.0  # L2 ties with L3: system_id order decides
    s = ch.select_systems(c, cov, p)
    assert list(s[s["role"] == "target"]["system_id"]) == ["L2"]


def test_classify_cases():
    c = ch.classify
    assert c(1, np.nan, np.nan, np.nan, 1.0, P).startswith("unresolved")
    assert c(2, 0.5, 10, np.nan, 1.0, P) == "too compact for the core mask"
    assert c(2, 1.5, 10, 2.2, 1.0, P) == "pair mismatch"
    assert c(4, 1.5, 10, 2.2, 1.0, P) == "light between images"  # quads are not pair-checked
    assert c(2, 1.5, 3.0, np.nan, 1.0, P) == "none between images"
    assert c(2, 1.5, -6.0, np.nan, 1.0, P).startswith("fit inconclusive")
    assert c(2, 1.5, 6.0, np.nan, np.nan, P) == "noise not estimable"
    assert c(2, 1.5, 5.0, 1.6, 1.0, P) == "light between images"


def test_control_gate():
    light, none = ch.DECIDED
    nine = [light] * 9 + [none]
    assert ch.control_gate(nine, P)["passed"]
    assert not ch.control_gate([light] * 8 + [none] * 2, P)["passed"]
    assert not ch.control_gate([light] * 4, P)["passed"]  # too few decided controls
    g = ch.control_gate([light] * 5 + ["unresolved (< 2 point sources)"] * 5, P)
    assert g["passed"] and g["n_decided"] == 5 and g["efficiency"] == 1.0


def _scene(lens_mag=None, reff=0.3, seed=1, quad=False):
    rng = np.random.default_rng(seed)
    n = 120
    yy, xx = np.indices((n, n))
    zp = 25.9
    sigma = 0.01
    img = rng.normal(0.0, sigma, (n, n))
    pts = [(45.0, 60.0, 50.0), (75.0, 60.0, 25.0)]
    if quad:
        pts = [(45.0, 60.0, 50.0), (75.0, 60.0, 30.0), (60.0, 45.0, 20.0), (60.0, 76.0, 15.0)]
    for x, y, a in pts:
        img += models.Moffat2D(a, x, y, gamma=1.6, alpha=2.8)(xx, yy)
    if lens_mag is not None:
        lx, ly = (60.0, 60.0)
        psf = ch.moffat_kernel(1.6, 2.8)
        g = ch.sersic_unit(img.shape, lx, ly, reff / PIX, 4, 0.2, psf=psf)
        img += 10 ** (-0.4 * (lens_mag - zp)) * g
    return img, sigma, zp, pts


def _noise_fn(sigma):
    return lambda npix: sigma


def test_find_point_sources_skips_extended_lens():
    img, sigma, _, pts = _scene(lens_mag=20.0, reff=0.4)
    pk = ch.find_point_sources(img, (60, 60), PIX, sigma, P)
    assert len(pk) == 2
    for x, y, _ in pts:
        assert min(np.hypot(x - a, y - b) for a, b in pk) <= 1.5


@pytest.mark.parametrize("quad", [False, True])
def test_measure_detects_lens_and_not_empty(quad):
    valid = np.ones((120, 120), bool)
    img, sigma, zp, _ = _scene(lens_mag=None, quad=quad)
    base = ch.measure(img, valid, (60, 60), PIX, P, sigma, _noise_fn(sigma))
    assert base["n_img"] == (4 if quad else 2)
    assert base["status"] == "none between images"
    img2, _, _, _ = _scene(lens_mag=21.0, quad=quad)
    m = ch.measure(img2, valid, (60, 60), PIX, P, sigma, _noise_fn(sigma))
    assert m["status"] == "light between images"
    assert ch.recovered(base, m, P)


def test_remove_halos_n_equals_d064_for_a_pair():
    rng = np.random.default_rng(3)
    res = rng.normal(0, 1, (60, 60))
    pos = [(20.3, 30.1), (41.7, 29.4)]
    fit = SimpleNamespace(
        x_0_1=SimpleNamespace(value=pos[0][0]), y_0_1=SimpleNamespace(value=pos[0][1]),
        x_0_2=SimpleNamespace(value=pos[1][0]), y_0_2=SimpleNamespace(value=pos[1][1]),
    )  # fmt: skip
    assert np.allclose(ch.remove_halos_n(res, pos), hst.remove_halos(res, fit))


def test_lens_position_sis_and_centroid():
    x, y = ch.lens_position([(0.0, 0.0), (10.0, 0.0)], [2.0, 1.0])
    assert x == pytest.approx(10 - 10 * 0.5 / 1.5) and y == 0.0
    assert ch.lens_position([(0, 0), (2, 0), (0, 2), (2, 2)], [1, 1, 1, 1]) == (1.0, 1.0)


def test_point_in_region():
    reg = "POLYGON 10.0 -1.0 10.1 -1.0 10.1 -0.9 10.0 -0.9"
    assert cli.point_in_region(10.05, -0.95, reg)
    assert not cli.point_in_region(10.15, -0.95, reg)
    wrap = "POLYGON 359.95 0.0 0.05 0.0 0.05 0.1 359.95 0.1"
    assert cli.point_in_region(0.01, 0.05, wrap)


def test_predicted_f814w_monotonic():
    a = ch.predicted_f814w(0.5, 0.5, 2.0, P)["m_typ"]
    b = ch.predicted_f814w(1.0, 0.5, 2.0, P)["m_typ"]
    assert b < a  # a larger Einstein radius needs a more luminous lens
    assert ch.predicted_f814w(1.0, 0.5, 2.0, P)["m_faint"] == pytest.approx(b + 2 * 0.88)
