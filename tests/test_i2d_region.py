"""Tests for scripts/i2d_region.py (synthetic mosaics)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("i2d_region", _DIR / "i2d_region.py")
ir = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ir)


def _wcs(n: int, scale_arcsec: float, ra: float = 3.58, dec: float = -30.39) -> WCS:
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [ra, dec]
    w.wcs.crpix = [n / 2 + 0.5, n / 2 + 0.5]
    w.wcs.cdelt = [-scale_arcsec / 3600, scale_arcsec / 3600]
    return w


def _mosaic(path: Path, n: int, scale: float, **kw: float) -> np.ndarray:
    rng = np.random.default_rng(1)
    sci = rng.normal(size=(n, n)).astype(np.float32)
    hdr = _wcs(n, scale, **kw).to_header()
    fits.HDUList(
        [
            fits.PrimaryHDU(header=fits.Header({"TELESCOP": "JWST"})),
            fits.ImageHDU(
                sci, header=hdr + fits.Header({"S_REGION": "POLYGON ICRS 0 0"}), name="SCI"
            ),
            fits.ImageHDU(sci * 0 + 0.5, header=hdr, name="ERR"),
        ]
    ).writeto(path)
    return sci


def test_cut_keeps_pixels_and_sky_positions(tmp_path: Path) -> None:
    big = tmp_path / "big.fits"
    sci = _mosaic(big, 400, 0.03)
    ref = tmp_path / "ref.fits"
    _mosaic(ref, 60, 0.06, ra=3.5805)  # 3.6" box, offset east of the big mosaic's centre
    out = tmp_path / "cut.fits"
    x0, x1, y0, y1 = ir.cut(str(big), ir.footprint_corners(str(ref)), out, pad=5)
    with fits.open(out) as h:
        assert h[0].header["CUTFROM"] == "big"
        assert "S_REGION" not in h["SCI"].header
        assert (h[0].header["CUTX0"], h[0].header["CUTY0"]) == (x0, y0)
        np.testing.assert_array_equal(h["SCI"].data, sci[y0:y1, x0:x1])
        assert np.all(h["ERR"].data == 0.5)
        w_cut = WCS(h["SCI"].header)
    w_big = WCS(fits.getheader(big, "SCI"))
    # the same sky position lands on the same pixel, shifted by the box origin
    ra, dec = w_big.all_pix2world([[x0 + 7.3, y0 + 11.6]], 0)[0]
    xy = w_cut.all_world2pix([[ra, dec]], 0)[0]
    np.testing.assert_allclose(xy, [7.3, 11.6], atol=1e-6)
    # the box covers the reference footprint (3.6" / 0.03" = 120 px, plus 2 x 5 px pad)
    assert 125 <= x1 - x0 <= 135 and 125 <= y1 - y0 <= 135


def test_cut_rejects_disjoint_reference(tmp_path: Path) -> None:
    big = tmp_path / "big.fits"
    _mosaic(big, 100, 0.03)
    ref = tmp_path / "ref.fits"
    _mosaic(ref, 50, 0.03, ra=10.0)
    with pytest.raises(ValueError, match="does not overlap"):
        ir.cut(str(big), ir.footprint_corners(str(ref)), tmp_path / "cut.fits")


@pytest.mark.network
def test_cut_remote_mosaic(tmp_path: Path) -> None:
    """Byte-range cut of a public JWST mosaic on S3 (jw02756 o003 F277W, 306 MB) to a 6" box."""
    src = (
        "s3://stpubdata/jwst/public/jw02756/L3/t/o003/jw02756-o003_t001_nircam_clear-f277w_i2d.fits"
    )
    ref = tmp_path / "ref.fits"
    _mosaic(ref, 100, 0.06, ra=3.5883, dec=-30.3972)  # 6" box on the Abell 2744 core
    out = tmp_path / "cut.fits"
    x0, x1, y0, y1 = ir.cut(src, ir.footprint_corners(str(ref)), out, pad=0)
    with fits.open(out) as h:
        assert h["SCI"].data.shape == (y1 - y0, x1 - x0)
        assert 90 <= x1 - x0 <= 150 and np.isfinite(h["SCI"].data).mean() > 0.9
        w = WCS(h["SCI"].header)
    x, y = w.all_world2pix([[3.5883, -30.3972]], 0)[0]
    assert 0 < x < x1 - x0 and 0 < y < y1 - y0
