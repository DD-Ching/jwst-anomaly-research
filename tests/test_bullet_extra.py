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


def _bullet_header():
    from astropy.io import fits

    h = fits.Header()
    h.update(
        NAXIS=2, NAXIS1=200, NAXIS2=160, CRPIX1=100.5, CRPIX2=80.5, CRVAL1=RA0, CRVAL2=DEC0,
        CTYPE1="RA---TAN", CTYPE2="DEC--TAN", CUNIT1="deg", CUNIT2="deg",
        PC1_1=-1.0, PC1_2=0.0, PC2_1=0.0, PC2_2=1.0,
        CD1_1=-2.2222222222222e-05, CD1_2=0.0, CD2_1=0.0, CD2_2=2.22222222222222e-05,
    )  # fmt: skip
    return h


def test_header_wcs_drops_cd_and_sets_0p04_arcsec():
    from astropy.wcs.utils import proj_plane_pixel_scales

    w = m.header_wcs(_bullet_header())
    assert np.allclose(proj_plane_pixel_scales(w) * 3600.0, 0.04)
    # header CRPIX is 1-based: 0-based pixel (99.5, 79.5) is CRVAL
    ra, dec = w.pixel_to_world_values(99.5, 79.5)
    assert abs(ra - RA0) < 1e-9 and abs(dec - DEC0) < 1e-9
    # east is to -x
    ra1, _ = w.pixel_to_world_values(100.5, 79.5)
    assert ra1 < ra


@pytest.mark.parametrize("origin", [0, 1])
def test_wcs_catalogue_check_finds_the_catalogue_convention(origin):
    w = m.header_wcs(_bullet_header())
    rng = np.random.default_rng(1)
    x0, y0 = rng.uniform(5, 190, 300), rng.uniform(5, 150, 300)
    ra, dec = w.pixel_to_world_values(x0, y0)
    r = m.wcs_catalogue_check(w, x0 + origin, y0 + origin, ra, dec)
    assert r["catalogue_origin"] == origin and r["passes"]
    assert r[f"origin_{origin}"]["rms_arcsec"] < 1e-6
    assert not r[f"origin_{1 - origin}"]["passes"]  # one pixel = 0.04" rms > 0.01"
    # a 0.08" scale catalogue fails both readings
    assert not m.wcs_catalogue_check(w, 2 * x0, 2 * y0, ra, dec)["passes"]


def test_centroid_offsets_recover_a_subpixel_shift():
    yy, xx = np.mgrid[0:60, 0:60]
    xs, ys = np.array([20.3, 40.0]), np.array([30.0, 15.7])
    d = 5.0 + sum(
        100 * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * 1.2**2))
        for x, y in zip(xs, ys, strict=True)
    )
    off = m.centroid_offsets(d, xs - 0.2, ys + 0.1)  # positions given wrong by (-0.2, +0.1)
    assert np.allclose(off, [[0.2, -0.1], [0.2, -0.1]], atol=0.03)
    assert np.isnan(m.centroid_offsets(d, np.array([1.0]), np.array([1.0]))).all()  # off edge
