import io
import sys
import zipfile
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import macs1931_gatek as mg  # noqa: E402

SOURCEPOS = (
    " Image Pos_x Pos_y Avg_Pos_x Avg_Pos_y Bary_Pos_x Bary_Pos_y\n"
    " 1 1449.1339 1464.5245 1449.0907 1464.6595 1448.1954 1462.9011\n"
    " 100 1343.2266 1428.2135 1343.1387 1428.2612 1341.5558 1428.9364\n"
)


def _header(n=41):
    h = fits.Header()
    h.update(
        CTYPE1="RA---TAN",
        CTYPE2="DEC--TAN",
        CRVAL1=292.9570833,
        CRVAL2=-26.5761,
        CRPIX1=(n + 1) / 2,
        CRPIX2=(n + 1) / 2,
        CDELT1=-mg.PIXEL_ARCSEC / 3600,
        CDELT2=mg.PIXEL_ARCSEC / 3600,
    )
    return h


def test_k_dlds_check_reproduces_authors():
    r = mg.k_dlds_check()
    assert r["passes"]
    assert [row["ours"] for row in r["rows"]] == pytest.approx(
        [0.89043, 0.99885, 1.03257], abs=5e-5
    )


def test_parse_sourcepos_bary_columns():
    s = mg.parse_sourcepos(SOURCEPOS)
    assert s == {"1": (1448.1954, 1462.9011), "100": (1341.5558, 1428.9364)}
    with pytest.raises(ValueError):
        mg.parse_sourcepos("nothing here\n")


def test_unit_check_needs_exactly_one_reading():
    one = mg.unit_check({"A_arcsec": {"1": 0.1, "2": 0.3}, "P_pixels": {"1": 9.0, "2": 8.0}})
    assert one["passes"] and one["reading"] == "A_arcsec"
    both = mg.unit_check({"A_arcsec": {"1": 0.1}, "P_pixels": {"1": 0.2}})
    assert not both["passes"] and both["reading"] is None
    none = mg.unit_check({"A_arcsec": {"1": 1.0}, "P_pixels": {"1": 0.4, "2": float("nan")}})
    assert not none["passes"]  # a missing system (NaN) never accepts a reading


def test_build_model_scales_to_dlsds1():
    n = 41
    ax = np.full((n, n), 2.0)
    ay = np.zeros((n, n))
    m = mg.build_model(ax, ay, _header(n), mg.PIXEL_ARCSEC, "x", "t")
    a, _ = m.deflection_xy(np.array([0.0]), np.array([0.0]))
    assert a[0] == pytest.approx(2.0 * mg.PIXEL_ARCSEC / mg.dls_ds(mg.Z_NORM), rel=1e-6)


def test_authors_sources_reference_pixel_is_origin():
    n = 41
    m = mg.build_model(np.zeros((n, n)), np.zeros((n, n)), _header(n), 1.0, "x", "t")
    c = (n + 1) / 2
    src = mg.authors_sources(m, {"1": (c, c), "2": (c + 10, c)})
    assert src["1"] == pytest.approx((0.0, 0.0), abs=1e-6)
    assert src["2"][0] == pytest.approx(10 * mg.PIXEL_ARCSEC, rel=1e-3)  # +i is West = +x


def test_read_zip_map_single_member(tmp_path):
    buf = io.BytesIO()
    fits.PrimaryHDU(np.ones((3, 3), np.float32), header=_header(3)).writeto(buf)
    p = tmp_path / "Alpha_x.fits.zip"
    with zipfile.ZipFile(p, "w") as z:
        z.writestr("Alpha_x.fits", buf.getvalue())
    data, hdr = mg.read_zip_map(p)
    assert data.shape == (3, 3) and hdr["CTYPE1"] == "RA---TAN"
    with zipfile.ZipFile(p, "a") as z:
        z.writestr("other.fits", buf.getvalue())
    with pytest.raises(ValueError):
        mg.read_zip_map(p)
