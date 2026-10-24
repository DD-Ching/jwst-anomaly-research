"""Frozen stage 3d continuum criterion and F_c of scripts/m1931_stage3d.py (offline)."""

import sys
from pathlib import Path

import numpy as np
from astropy.wcs import WCS

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import m1931_stage3d as s3d  # noqa: E402

PIX = 0.04


def _wcs(n=400):
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [292.95, -26.57]
    w.wcs.crpix = [n / 2, n / 2]
    w.wcs.cdelt = [-PIX / 3600.0, PIX / 3600.0]
    return w


def _grid(n=12):
    ra = 292.95 + (np.arange(n) - n / 2) * 2.0 / 3600.0
    return ra, np.full(n, -26.57)


def test_continuum_class():
    assert s3d.continuum_class(5.0) == "continuum"
    assert s3d.continuum_class(4.99) == "faint"
    assert s3d.continuum_class(float("nan")) == s3d.OFF_SECTION


def test_continuum_snr_point_source_and_blank():
    n = 400
    rng = np.random.default_rng(1)
    d = rng.normal(0.0, 1.0, (n, n))
    w = _wcs(n)
    d[200, 200] += 200.0  # a bright pixel at the centre
    ra, dec = w.pixel_to_world_values([200.0, 260.0, -50.0], [200.0, 200.0, 200.0])
    snr = s3d.continuum_snr(d, w, ra, dec, np.array([]), np.array([]), PIX)
    assert snr[0] >= s3d.CONT_MIN_SNR
    assert abs(snr[1]) < s3d.CONT_MIN_SNR
    assert np.isnan(snr[2])  # off the section


def test_frame_check_fc_counts_only_continuum_images():
    ra, dec = _grid(12)
    snr = [10.0] * 10 + [1.0, float("nan")]
    f = s3d.frame_check_fc(ra, dec, snr, ra[:5], dec[:5])
    assert f["n_continuum"] == 10 and f["n_matched"] == 5
    assert f["continuum_classes"] == {"continuum": 10, "faint": 1, "off_section": 1}
    assert f["passes"]  # 5 >= 0.5 x 10, medians 0, no shifted match
    g = s3d.frame_check_fc(ra, dec, snr, ra[:4], dec[:4])
    assert not g["passes"]  # 4 < 5


def test_frame_check_fc_needs_min_n_and_chance():
    ra, dec = _grid(12)
    few = [10.0] * 9 + [1.0] * 3
    assert not s3d.frame_check_fc(ra, dec, few, ra, dec)["passes"]  # n_c 9 < 10
    snr = [10.0] * 12
    north = dec + 5.0 / 3600.0  # rows also 5" north: every shifted position matches at PA 0
    cra, cdec = np.concatenate([ra, ra]), np.concatenate([dec, north])
    f = s3d.frame_check_fc(ra, dec, snr, cra, cdec)
    assert f["passes_matcher_medians"] and f["chance_shifted_counts"][0] == 12
    assert f["chance_passes"]  # 12 >= 2 x 12 / 8


def test_frame_check_fc_no_continuum():
    ra, dec = _grid(12)
    f = s3d.frame_check_fc(ra, dec, [1.0] * 12, ra, dec)
    assert f["n_continuum"] == 0 and not f["passes"]


def test_matched_each():
    ra, dec = _grid(4)
    m = s3d.matched_each(ra, dec, ra[:2], dec[:2] + 0.4 / 3600.0)
    assert m.tolist() == [True, True, False, False]


def test_use_fc_replaces_rule5_and_checks_catalogue():
    import pytest

    stored = {"catalogue": {"n_rows_kept": 3}, "frame_F_c": {"passes": True}}
    f = {"passes": False}
    s3d.use_fc(stored)(f, None, {"ra": [1, 2, 3]})
    assert f["passes"] and not f["passes_stage3_rule"]
    with pytest.raises(ValueError):
        s3d.use_fc(stored)({"passes": False}, None, {"ra": [1, 2]})
