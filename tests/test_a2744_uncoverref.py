"""Unit tests of the Abell 2744 UNCOVER image-list copy (PLAN.md Inputs)."""

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

from jwst_anomaly import lensmodel

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import a2744_uncover as au  # noqa: E402
from a2744_uncoverref import (  # noqa: E402
    fold_deg,
    image_list,
    lenstool_id,
    segment_angle,
    theta_check,
    theta_sample,
    y_axis_pa_deg,
)


def test_lenstool_id_keeps_the_knot_level_source():
    for uid in ("1.11", "1.13", "78.21", "83.1", "83.3", "12.2"):
        assert lensmodel.image_family(lenstool_id(uid)) == au.uncover_system(uid)
    assert lenstool_id("78.22") == "78.2b"
    assert lenstool_id("83.3") == "83c"


def test_image_list_flags_spec_and_drops_none():
    raw = Table(
        {
            "ID": [1.11, 1.12, 83.1, 83.2, 90.1],
            "RA": [3.5, 3.6, 3.7, 3.8, 3.9],
            "Dec": [-30.4] * 5,
            "z": [2.6, 2.6, 4.0, 4.0, 0.0],
            "Redshift_type": ["spec", "spec", "phot", "phot", "None"],
        }
    )
    t = image_list(raw)
    assert list(t["image_id"]) == ["1.1a", "1.1b", "83a", "83b"]
    assert list(t["z_spec"]) == [True, True, False, False]
    assert np.all(np.isfinite(t["ra"]))


def test_ambiguous_ids_fail_loudly():
    with pytest.raises(ValueError):
        lenstool_id("34.1" + "0")
    raw = Table(
        {
            "ID": [1.11, 1.11],
            "RA": [3.5, 3.6],
            "Dec": [-30.4] * 2,
            "z": [2.6, 2.6],
            "Redshift_type": ["spec", "spec"],
        }
    )
    with pytest.raises(ValueError):
        image_list(raw)


def test_gate_k_wls_decides_with_the_unchanged_gate(monkeypatch):
    import a2744_uncoverref as ar

    sp = Table({"image_id": ["1a", "1b", "2a", "2b", "2c"], "system": ["1", "1", "2", "2", "2"]})
    monkeypatch.setattr(ar.lensmodel, "backtrace_images", lambda model, t, zml: t)

    def fake_imageplane(model, grid, bt, d=(0.2, 0.3, 0.1, 0.4, 0.5)):
        return Table(
            {
                "image_id": bt["image_id"],
                "system": bt["system"],
                "dtheta_arcsec": list(d),
                "n_predicted": [3, 3, 2, 2, 2],
                "shared_match": [False, False, False, True, True],
            }
        )

    monkeypatch.setattr(ar.au, "imageplane", fake_imageplane)
    _, k = ar.gate_k_wls(None, None, sp)
    assert k["gate_K"]["passes"] and k["gate_K"]["n_images"] == 5
    assert k["not_reproduced"] == {} and k["shared_match"] == 2
    assert k["extra_predicted"] == {"1": (2, 3)} and k["fewer_predicted"] == {"2": (3, 2)}

    # one image beyond 3.0'' (80 % matched still passes on fraction) but rms > 1.0'' fails
    monkeypatch.setattr(
        ar.au, "imageplane", lambda m, g, bt: fake_imageplane(m, g, bt, (0.2, 0.3, 0.1, 0.4, 5.0))
    )
    _, k = ar.gate_k_wls(None, None, sp)
    assert not k["gate_K"]["passes"] and k["gate_K"]["matched_fraction"] == 0.8
    assert k["not_reproduced"] == {"2c": 5.0}

    # an unmatched (NaN) image counts against the fraction
    monkeypatch.setattr(
        ar.au,
        "imageplane",
        lambda m, g, bt: fake_imageplane(m, g, bt, (0.2, np.nan, 0.1, np.nan, 0.5)),
    )
    _, k = ar.gate_k_wls(None, None, sp)
    assert not k["gate_K"]["passes"] and k["gate_K"]["matched_fraction"] == 0.6


def _ellipse_pixels(phi, a=12.0, b=4.0):
    yy, xx = np.mgrid[-30:31, -30:31]
    u = xx * np.cos(phi) + yy * np.sin(phi)
    v = -xx * np.sin(phi) + yy * np.cos(phi)
    return np.nonzero((u / a) ** 2 + (v / b) ** 2 <= 1.0)


@pytest.mark.parametrize("phi", [0.0, 0.4, np.pi / 2 - 0.1, -0.7])
def test_segment_angle_recovers_ccw_from_x(phi):
    ys, xs = _ellipse_pixels(phi)
    assert fold_deg(segment_angle(ys, xs) - phi) < 2.0  # pixelisation


def test_theta_check_adopts_the_matching_reading():
    rng = np.random.default_rng(1)
    theta = rng.uniform(-np.pi / 2, np.pi / 2, 200)
    noisy = theta + np.radians(rng.normal(0, 3, 200))
    r = theta_check(theta, noisy)
    assert r["H_a_ccw_from_x"]["passes"] and r["adopted"] == "H_a_ccw_from_x"
    r = theta_check(theta, np.pi / 2 - noisy)
    assert r["adopted"] == "H_b_cw_from_y"
    r = theta_check(theta, rng.uniform(-np.pi / 2, np.pi / 2, 200))
    assert r["adopted"] is None


def test_theta_check_needs_median_and_fraction():
    theta = np.zeros(10)
    # median 0 deg but only 60 % within 20 deg: rule 3's "and" fails
    phi = np.radians([0, 0, 0, 0, 0, 0, 45, 45, 45, 45])
    assert not theta_check(theta, phi)["H_a_ccw_from_x"]["passes"]
    phi = np.radians([0, 0, 0, 0, 0, 0, 0, 45, 45, 45])
    assert theta_check(theta, phi)["H_a_ccw_from_x"]["passes"]


def test_theta_sample_cuts_and_order():
    t = Table(
        {
            "id": [1, 2, 3, 4, 5, 6],
            "a_image": [5.0, 8.0, 8.0, 3.0, 9.0, 6.0],
            "b_image": [2.0, 2.0, 2.0, 1.0, 8.0, 2.0],
            "use_phot": [1, 1, 1, 1, 1, 0],
            "flag_kron": [0, 0, 0, 0, 0, 0],
            "flag_nearbcg": [0, 0, 0, 0, 0, 0],
        }
    )
    # 4: a < 4 px; 5: ell 0.11; 6: use_phot 0; 2 and 3 tie on a -> by id
    assert list(theta_sample(t)) == [1, 2, 0]


def test_y_axis_pa_of_a_north_up_and_rotated_grid():
    x, y = np.meshgrid(np.arange(0, 100, 10.0), np.arange(0, 100, 10.0))
    x, y = x.ravel(), y.ravel()
    s = 0.04 / 3600.0
    dec0 = -30.4
    ra = 3.58 - x * s / np.cos(np.radians(dec0))  # +x west, +y north
    dec = dec0 + y * s
    assert abs(y_axis_pa_deg(x, y, ra, dec)) < 1e-6
    c, n = np.cos(np.radians(30)), np.sin(np.radians(30))
    east, north = (-x * c + y * n) * s, (x * n + y * c) * s  # +y at PA 30 deg
    assert abs(y_axis_pa_deg(x, y, 3.58 + east / np.cos(np.radians(dec0)), dec0 + north) - 30) < 0.1
