"""Tests of the Bullet Cluster extra-image script: constants, input pinning and frame check F
(PLAN.md rule 1)."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import bullet_extra as m  # noqa: E402
import m0417_extra as mx  # noqa: E402

RA0, DEC0 = 104.63, -55.95


def _grid(n=40):
    k = np.arange(n)
    return RA0 + 0.001 * (k % 8), DEC0 + 0.001 * (k // 8)


def _shift(ra, dec, dx, dy):
    """Move positions by (dx along RA cos dec, dy) arcsec (small-offset approximation)."""
    return ra + dx / 3600.0 / np.cos(np.radians(dec)), dec + dy / 3600.0


def test_constants_name_the_planned_field_model_and_tile():
    assert (m.FIELD, m.MODEL, m.OBS, m.TILE) == ("bullet", "bullet-canucs-maps", "o001", "t001")
    assert (m.F_MATCH_ARCSEC, m.F_MAX_MEDIAN_ARCSEC, m.F_MIN_N) == (0.5, 0.1, 30)
    assert (m.GATEK / m.GOLD_REL).resolve() == (
        Path("results/lens_consistency/bullet_inputs/bullet_gold_images.ecsv").resolve()
    )


def test_register_points_the_frame_tables_at_go4598_tile():
    m.register()
    assert mx.FIELDS["bullet"] == (m.GOLD_REL, m.MODEL, "o001", "t001")
    assert mx.S3.format("o001", "t001", "i2d.fits").endswith(
        "jw04598/L3/t/o001/jw04598-o001_t001_nircam_clear-f200w_i2d.fits"
    )


def test_refuses_to_run_until_the_canucs_files_are_pinned(tmp_path, monkeypatch):
    monkeypatch.setattr(m, "CAT_SHA", None)
    with pytest.raises(SystemExit, match="pin CAT_SHA"):
        m.check_inputs(tmp_path / "c", tmp_path / "i")


def test_refuses_files_with_other_content(tmp_path, monkeypatch):
    f = tmp_path / "f"
    f.write_bytes(b"x")
    monkeypatch.setattr(m, "CAT_SHA", "0" * 64)
    monkeypatch.setattr(m, "IMAGE_SHA", "0" * 64)
    with pytest.raises(SystemExit, match="not the CANUCS"):
        m.check_inputs(f, f)


def test_frame_check_passes_for_a_small_shift_and_reports_its_sign():
    cra, cdec = _grid()
    gra, gdec = _shift(cra, cdec, 0.05, -0.08)  # gold = catalogue + (0.05, -0.08)
    r = m.frame_check_f(gra, gdec, cra, cdec)
    assert r["passes"] and r["n_matched"] == 40
    assert r["median_dx_arcsec"] == pytest.approx(0.05, abs=2e-3)
    assert r["median_dy_arcsec"] == pytest.approx(-0.08, abs=2e-3)


def test_frame_check_fails_on_either_axis_beyond_0p1():
    cra, cdec = _grid()
    for dx, dy in ((0.15, 0.0), (0.0, -0.15)):
        r = m.frame_check_f(*_shift(cra, cdec, dx, dy), cra, cdec)
        assert r["n_matched"] == 40 and not r["passes"]


def test_frame_check_needs_30_matches_within_0p5():
    cra, cdec = _grid()
    gra, gdec = _shift(cra, cdec, 0.0, 0.0)
    gra = gra.copy()
    gra[:11] += 0.6 / 3600.0 / np.cos(np.radians(DEC0))  # 11 beyond 0.5": 29 left
    r = m.frame_check_f(gra, gdec, cra, cdec)
    assert r["n_matched"] == 29 and not r["passes"]


def test_frame_check_ignores_non_finite_catalogue_rows():
    cra, cdec = _grid()
    cra = np.r_[cra, np.nan]
    cdec = np.r_[cdec, DEC0]
    r = m.frame_check_f(*_grid(), cra, cdec)
    assert r["passes"] and r["median_dx_arcsec"] == pytest.approx(0.0, abs=1e-6)
