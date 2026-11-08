import importlib.util
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

_P = Path(__file__).resolve().parents[1] / "scripts" / "arcs_shear_error_2157.py"
_spec = importlib.util.spec_from_file_location("arcs_shear_error_2157", _P)
m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m)


def _image(tmp_path, ny=60, nx=50):
    data = np.arange(ny * nx, dtype=np.float32).reshape(ny, nx) * 0.5
    hdr = fits.Header()
    hdr.update(
        CTYPE1="RA---TAN",
        CTYPE2="DEC--TAN",
        CRVAL1=47.24,
        CRVAL2=26.76,
        CRPIX1=20.0,
        CRPIX2=25.0,
        CD1_1=-1.666666667e-05,
        CD2_2=1.666666667e-05,
        CD1_2=0.0,
        CD2_1=0.0,
    )
    for k in range(60):  # pad the header past one 2880-byte block
        hdr[f"PAD{k:03d}"] = k
    path = tmp_path / "img.fits"
    fits.PrimaryHDU(data, header=hdr).writeto(path)
    raw = path.read_bytes()
    return data, raw, WCS(fits.getheader(path))


def _getter(raw):
    calls = []

    def get(start, end):
        calls.append((start, end))
        return raw[start : end + 1]

    get.calls = calls
    return get


def test_read_strip_matches_full_image_and_wcs(tmp_path):
    data, raw, wcs = _image(tmp_path)
    get = _getter(raw)
    strip, swcs = m.read_strip(get, 10, 31)
    assert np.array_equal(strip, data[10:31])
    assert len(get.calls) == 2  # header + one data range
    cut, cwcs = m.cut_columns(strip, swcs, 5, 40)
    assert np.array_equal(cut, data[10:31, 5:40])
    # the same sky position maps to the same value in both grids
    ra, dec = wcs.pixel_to_world_values(17.0, 22.0)
    i, j = cwcs.world_to_pixel_values(ra, dec)
    assert np.isclose(i, 12.0) and np.isclose(j, 12.0)
    assert cut[int(round(float(j))), int(round(float(i)))] == data[22, 17]
    assert np.isclose(m.pixel_scale_arcsec(cwcs), 0.06)


def test_read_strip_clips_and_rejects_short_reads(tmp_path):
    data, raw, _ = _image(tmp_path)
    strip, _ = m.read_strip(_getter(raw), -5, 999)
    assert np.array_equal(strip, data)

    def short(start, end):
        return raw[start : end + 1][:-4] if start > 0 else raw[start : end + 1]

    try:
        m.read_strip(short, 0, 10)
    except ValueError as e:
        assert "short read" in str(e)
    else:
        raise AssertionError("short read not detected")


def test_parse_header_offset(tmp_path):
    _, raw, _ = _image(tmp_path)
    hdr, off = m.parse_header(raw[: m.HEADER_BYTES])
    assert (
        off == len(raw) - -(-60 * 50 * 4 // 2880) * 2880 and off % 2880 == 0 and hdr["NAXIS1"] == 50
    )


def test_gates():
    assert m.v1_gate(1e-3) and not m.v1_gate(1.1e-3) and not m.v1_gate(float("nan"))
    assert m.v2_gate(["aligned", "anti"], ["aligned", "anti"], np.array([0.1, 0.5]))
    assert not m.v2_gate(["aligned", "anti"], ["aligned", "mixed"], np.array([0.1, 0.1]))
    assert not m.v2_gate(["aligned"], ["aligned"], np.array([0.6]))
    assert not m.v2_gate([], [], np.array([]))
