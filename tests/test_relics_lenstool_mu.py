"""Frozen rules of scripts/relics_lenstool_mu.py (H_K mu-region plan), offline."""

import sys
from pathlib import Path

import numpy as np
from astropy.wcs import WCS

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import relics_lenstool_mu as rm  # noqa: E402


def _wcs(scale_arcsec: float) -> WCS:
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [10.0, -20.0]
    w.wcs.crpix = [51, 51]
    w.wcs.cdelt = [-scale_arcsec / 3600.0, scale_arcsec / 3600.0]
    return w


def test_jacobian_recovers_sis_convergence_and_critical_curve():
    # SIS: alpha = theta_E theta / r, kappa = theta_E / (2 r), mu infinite at r = theta_E.
    scale, theta_e = 0.5, 10.0
    yy, xx = (np.indices((101, 101)) - 50) * scale
    r = np.hypot(xx, yy)
    r[50, 50] = np.nan
    ax, ay = theta_e * xx / r, theta_e * yy / r
    dxx, dxy, dyx, dyy = rm.jacobian(ax, ay, _wcs(scale))
    kappa = 0.5 * (dxx + dyy)
    sel = (r > 3) & (r < 20) & np.isfinite(kappa)
    assert np.nanmedian(np.abs(kappa[sel] - theta_e / (2 * r[sel]))) < 0.01
    det = (1 - dxx) * (1 - dyy) - dxy * dyx
    ring = (np.abs(r - theta_e) < 0.3) & np.isfinite(det)
    far = (r > 18) & (r < 22) & np.isfinite(det)
    assert np.nanmedian(np.abs(det[ring])) < 0.05 < np.nanmedian(np.abs(det[far]))


def test_region_stats_selects_high_mu_and_minimum_area():
    a_ref = np.array([10.0, 10.0, 10.0, 10.0])
    a_lt = np.array([15.0, 15.0, 30.0, 30.0])
    mu = np.array([12.0, -50.0, 3.0, 9.9])
    ok = np.array([True, True, True, True])
    s = rm.region_stats(a_lt, a_ref, mu, ok, rm.MU_MIN)
    assert np.isclose(s["R"], 1.5)
    assert s["region_arcsec2"] == 2.0
    assert s["region_fraction"] == 0.5
    assert rm.region_stats(a_lt, a_ref, mu, ok, 100.0)["R"] is None


def _row(r, frac=0.1, cls=None):
    import relics_lenstool_batch as rb

    return {"R": r, "region_fraction": frac, "class": cls or rb.classify(r)}


FACTORS = {"macs0553": 0.665, "a3192": 0.692, "rxc2211": 1.017, "cl0152": 1.045}


def test_gate_v_needs_class_tolerance_and_subregion():
    good = {"macs0553": _row(1.6), "a3192": _row(1.3), "rxc2211": _row(1.0), "cl0152": _row(1.0)}
    assert rm.gate_v(good, FACTORS)["passed"]
    wide = dict(good, rxc2211=_row(1.0, frac=0.9))
    assert not rm.gate_v(wide, FACTORS)["passed"]
    unmeasured = dict(
        good, cl0152={"R": None, "region_fraction": None, "class": "unmeasured (gate T)"}
    )
    assert not rm.gate_v(unmeasured, FACTORS)["passed"]
