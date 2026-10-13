"""Offline tests for scripts/sunrise_relics.py (synthetic SIS rows and maps; no network)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import sunrise_relics as sr  # noqa: E402


def _sis_rows(z_true: float, knot: str = "2", theta_e: float = 10.0, beta=(1.0, 0.5)):
    """SIS images (theta_e at D_LS/D_S = 1) of a point source at z_true, stored like map rows."""
    te = theta_e * sr.dratio(z_true)
    b = np.array(beta)
    u = b / np.hypot(*b)
    rows = []
    for k, sign in enumerate((1.0, -1.0)):
        th = b + sign * te * u
        r = np.hypot(*th)
        n = th / r
        rows.append(
            {
                "image": f"{knot}{'ab'[k]}",
                "x_arcsec": th[0],
                "y_arcsec": th[1],
                "alpha_x": theta_e * n[0],
                "alpha_y": theta_e * n[1],
                "psi_unit": list(
                    theta_e / r * np.array([1 - n[0] ** 2, -n[0] * n[1], 1 - n[1] ** 2])
                ),
            }
        )
    return rows


@pytest.mark.parametrize("z_true", [2.0, 3.0, 6.2])
def test_zscan_minimum_at_true_redshift(z_true):
    rows = _sis_rows(z_true)
    s = sr.zscan(rows)
    assert sr.Z_GRID[int(np.argmin(s))] == pytest.approx(z_true, abs=0.03)
    assert sr.source_rms(rows, z_true) < 1e-9
    assert max(sr.image_plane_offsets(rows, z_true)) < 1e-9


def test_source_rms_pools_knots_about_their_own_means():
    a = _sis_rows(6.2, knot="1.1", beta=(1.0, 0.5))
    b = _sis_rows(6.2, knot="1.2", beta=(-2.0, 0.3))  # a different source position
    assert sr.source_rms(a + b, 6.2) < 1e-9


def test_image_plane_offset_inverts_the_local_jacobian():
    z = 3.0
    rows = _sis_rows(z)
    delta = np.array([0.04, -0.03])  # source-plane shift of image a, arcsec
    a = dict(rows[0])
    a["alpha_x"] -= delta[0] / sr.dratio(z)
    a["alpha_y"] -= delta[1] / sr.dratio(z)
    off = sr.image_plane_offsets([a, rows[1]], z)
    # each image sits delta/2 from the mean source position
    for k, q in enumerate(rows):
        expected = np.hypot(*np.linalg.solve(sr.jacobian(q, z), delta / 2))
        assert off[k] == pytest.approx(expected, rel=1e-9)


def test_local_lens_on_sampled_sis_map():
    step, theta_e = 0.05, 10.0
    j, i = np.mgrid[-3:4, -3:4]
    x0, y0 = 6.0, 4.0  # patch centre, arcsec (+x West along +i, +y North along +j)
    x, y = x0 + i * step, y0 + j * step
    r = np.hypot(x, y)
    lens = sr.local_lens(theta_e * x / r, theta_e * y / r, step, 3.0, 3.0, 1.0)
    rc = np.hypot(x0, y0)
    assert lens["kappa"] == pytest.approx(theta_e / (2 * rc), rel=1e-3)
    assert lens["mu"] == pytest.approx(1.0 / (1.0 - theta_e / rc), rel=1e-3)


def test_frame_check_rejects_rotated_grid():
    from astropy.wcs import WCS

    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cdelt = [-0.04 / 3600, 0.04 / 3600]
    assert sr.frame_check(w) == pytest.approx(0.04)
    w.wcs.pc = [[np.cos(0.3), -np.sin(0.3)], [np.sin(0.3), np.cos(0.3)]]
    with pytest.raises(RuntimeError):
        sr.frame_check(w)
