"""Offline tests for scripts/sunrise_zfit.py (synthetic SIS rows; no network)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import sunrise_zfit as zf  # noqa: E402


def test_zref_round_trip():
    for z in (1.5, 3.0, 6.2, 9.0):
        s = zf.PIXEL * zf.dratio(zf.Z_MAPS) / zf.dratio(z)
        assert zf.zref_for_scale(s) == pytest.approx(z, rel=1e-6)
    assert zf.zref_for_scale(zf.PIXEL) == pytest.approx(zf.Z_MAPS, rel=1e-6)


def _sis_rows(z_true: float, zref: float = 3.0, theta_e_ref: float = 8.0, beta=(1.0, 0.5)):
    """SIS images of a point source at z_true, stored like the map rows (raw = alpha/PIXEL)."""
    te = theta_e_ref * zf.dratio(z_true) / zf.dratio(zref)
    b = np.array(beta)
    u = b / np.hypot(*b)
    rows = []
    for k, sign in enumerate((1.0, -1.0)):
        th = b + sign * te * u
        r = np.hypot(*th)
        n = th / r
        # SIS at zref: alpha = theta_e_ref * n; psi_ij = theta_e_ref / r * (delta_ij - n_i n_j)
        psi = theta_e_ref / r * np.array([1 - n[0] ** 2, -n[0] * n[1], 1 - n[1] ** 2])
        rows.append(
            {
                "image": f"9{'ab'[k]}",
                "x_arcsec": th[0],
                "y_arcsec": th[1],
                "alpha_x": theta_e_ref * n[0] / zf.PIXEL,
                "alpha_y": theta_e_ref * n[1] / zf.PIXEL,
                "psi": list(psi / zf.PIXEL),
                "constraint": True,
            }
        )
    return rows


@pytest.mark.parametrize("z_true", [2.0, 6.2])
def test_fit_redshift_recovers_sis_source_redshift(z_true):
    fit = zf.fit_redshift(_sis_rows(z_true), 3.0)
    assert fit["z_best"] == pytest.approx(z_true, abs=0.03)
    assert fit["rms_beta_best"] < 1e-2
    assert not fit["at_grid_edge"]


def test_image_plane_rms_is_zero_at_the_true_scale():
    rows = _sis_rows(3.0)
    out = zf.image_plane_rms(rows, zf.PIXEL)
    assert out["rms_beta"] == pytest.approx(0.0, abs=1e-9)
    assert out["rms_theta"] == pytest.approx(0.0, abs=1e-9)


def test_interval_and_grid_edge():
    rows = _sis_rows(6.2)
    fit = zf.fit_redshift(rows, 3.0)
    lo, hi = zf.interval(fit, 0.05)
    assert lo < 6.2 < hi
    assert zf.interval(fit, -1.0) is None


def test_grid_edge_flag_when_the_true_redshift_is_beyond_the_grid():
    fit = zf.fit_redshift(_sis_rows(6.2), 3.0, grid=np.arange(1.0, 4.0, 0.02))
    assert fit["at_grid_edge"]
    assert fit["z_best"] == pytest.approx(3.98)


def test_unmatched_reference_set_raises_and_impossible_scale_is_nan():
    with pytest.raises(ValueError):
        zf.image_plane_rms(_sis_rows(3.0), zf.PIXEL, {"no-such-image"})
    assert np.isnan(zf.zref_for_scale(0.01))
