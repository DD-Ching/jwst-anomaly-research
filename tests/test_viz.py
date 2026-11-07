"""Tests for viz: stretch, single-cutout PNGs and contact sheets from synthetic cutouts."""

from __future__ import annotations

import astropy.units as u
import matplotlib
import matplotlib.image as mpimg
import numpy as np
import pytest
from astropy.io import fits
from astropy.table import MaskedColumn, Table, vstack

from jwst_anomaly import viz


def _cutout(path, *, band="F200W", uid="s1", flag="ok", nan_corner=False):
    yy, xx = np.mgrid[0:31, 0:31]
    data = 0.01 * np.random.default_rng(0).normal(size=(31, 31))
    data += np.exp(-((xx - 15) ** 2 + (yy - 15) ** 2) / 8.0)
    if nan_corner:
        data[:5, :5] = np.nan
    head = fits.Header({"SRCUID": uid, "BAND": band, "QUALFLAG": flag})
    fits.HDUList(
        [fits.PrimaryHDU(header=head), fits.ImageHDU(data.astype(np.float32), name="SCI")]
    ).writeto(path)
    return path


def _cutout_table(tmp_path, bands=("F200W",), n=3):
    rows = []
    for i in range(n):
        for band in bands:
            uid = f"src{i}"
            flag = "edge" if i == 1 else "ok"
            path = _cutout(tmp_path / f"{uid}_{band}.fits", band=band, uid=uid, flag=flag)
            rows.append((uid, band, path.as_posix(), 0.0, i == 1, flag))
    rows.append(("far", bands[-1], "", 1.0, True, "outside"))
    t = Table(
        rows=rows, names=("source_uid", "band", "path", "frac_nan", "on_edge", "quality_flag")
    )
    t.meta.update(provenance="observed", source="synthetic")
    return t


def test_normalize_keeps_nan_masked():
    data = np.array([[np.nan, 0.0], [1.0, 100.0]])
    scaled = viz.normalize(data)(np.ma.masked_invalid(data))
    assert scaled.mask[0, 0] and not scaled.mask[1].any()  # NaN -> colormap "bad" colour
    assert 0 < scaled[1, 0] < 1
    assert viz.normalize(np.full((3, 3), np.nan)) is None
    assert viz.normalize(np.zeros((3, 3))) is not None  # flat image does not divide by zero


def test_cutout_png_default_name_and_no_data_colour(tmp_path):
    path = _cutout(tmp_path / "a.fits", nan_corner=True)
    out = viz.cutout_png(path)
    assert out == tmp_path / "a.png"
    img = mpimg.imread(out)
    assert img.ndim == 3 and img.shape[0] > 100 and img.std() > 0
    blue = np.array(matplotlib.colors.to_rgb(viz.NO_DATA_COLOR))
    assert (np.abs(img[..., :3] - blue).max(axis=-1) < 0.02).sum() > 100  # NaN corner is blue


@pytest.mark.parametrize("bands", [("F200W",), ("F770W", "F200W")])
def test_contact_sheet(tmp_path, bands):
    table = _cutout_table(tmp_path, bands=bands)
    out = viz.contact_sheet(
        table, tmp_path / "sheet.png", ranks={"src2": 1, "src0": 2}, title="synthetic"
    )
    img = mpimg.imread(out)
    assert img.shape[0] > 200 and img.std() > 0


def test_ordering_and_band_sort(tmp_path):
    table = _cutout_table(tmp_path, bands=("F770W", "F200W"))
    assert viz._ordered_sources(table, {"src2": 1, "src0": 2}) == ["src2", "src0", "src1", "far"]
    assert viz._ordered_sources(table, None) == ["src0", "src1", "src2", "far"]
    assert sorted({"F770W", "F1000W", "F200W"}, key=viz._band_sort_key) == [
        "F200W",
        "F770W",
        "F1000W",
    ]
    table["rank"] = [3, 3, 1, 1, 2, 2, 9]
    assert viz._ordered_sources(table, None) == ["src1", "src2", "src0", "far"]


def test_duplicate_cells_prefer_rows_with_data_and_masked_ranks(tmp_path):
    table = _cutout_table(tmp_path)
    dup = table[:1].copy()
    dup["path"], dup["quality_flag"] = "", "outside"
    table = vstack([table, dup])  # same uid/band from a second image, after the real one
    assert viz._pick_cells(table)[("src0", "F200W")]["path"]
    table["rank"] = MaskedColumn([2, 1, 3, 9, 0], mask=[False, False, True, True, False])
    assert viz._ordered_sources(table, None) == ["src1", "src0", "src2", "far"]
    viz.contact_sheet(table, tmp_path / "sheet.png")


def test_contact_sheet_rejects_empty(tmp_path):
    empty = _cutout_table(tmp_path)[:0]
    with pytest.raises(ValueError, match="no cutouts"):
        viz.contact_sheet(empty, tmp_path / "x.png")


def _sky_wcs(theta_deg, *, flip=False, crpix=(16.0, 16.0), dec=-73.45, scale_arcsec=0.03):
    """TAN WCS rotated by ``theta_deg`` (east of north); ``flip`` makes east point right."""
    from astropy.wcs import WCS

    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [110.8, dec]
    w.wcs.crpix = list(crpix)
    s = scale_arcsec / 3600.0
    t = np.radians(theta_deg)
    rot = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    w.wcs.cd = rot @ np.diag([s if flip else -s, s])
    return w


def _offset_pixel(w, shape, east_arcsec, north_arcsec):
    """Input pixel of the sky point offset (east, north) from the image centre."""
    from astropy.coordinates import SkyCoord

    ny, nx = shape
    c = SkyCoord(w.pixel_to_world((nx - 1) / 2, (ny - 1) / 2))
    target = c.spherical_offsets_by(east_arcsec * u.arcsec, north_arcsec * u.arcsec)
    x, y = w.world_to_pixel(target)
    return int(round(float(y))), int(round(float(x)))


@pytest.mark.parametrize(
    ("theta", "flip", "crpix"),
    [
        (0.0, False, (16.0, 16.0)),
        (37.0, False, (16.0, 16.0)),
        (-121.0, False, (16.0, 16.0)),
        (90.0, True, (16.0, 16.0)),
        (210.0, True, (16.0, 16.0)),
        # Mosaic-like: cutout 12' from the tangent point at dec -73 (meridian convergence).
        (25.0, False, (-24000.0, 9000.0)),
    ],
)
def test_north_up_puts_north_up_and_east_left(theta, flip, crpix):
    w = _sky_wcs(theta, flip=flip, crpix=crpix)
    shape = (31, 31)
    k = 9  # pixels, i.e. 0.27"
    for (east, north), expect in (((0.0, k * 0.03), (15 + k, 15)), ((k * 0.03, 0.0), (15, 15 - k))):
        yy, xx = np.mgrid[0 : shape[0], 0 : shape[1]]
        py, px = _offset_pixel(w, shape, east, north)
        data = np.exp(-((xx - 15) ** 2 + (yy - 15) ** 2) / 4.0)  # target
        data += 5.0 * np.exp(-((xx - px) ** 2 + (yy - py) ** 2) / 4.0)  # neighbour
        out = viz.north_up(data, w)
        assert out.shape == shape
        assert np.unravel_index(
            np.nanargmax(np.where(np.hypot(xx - 15, yy - 15) < 3, out, -1)), shape
        ) == (15, 15)
        peak = np.unravel_index(np.nanargmax(out), shape)
        # within one pixel: the neighbour was placed on the nearest input pixel
        assert max(abs(int(peak[0]) - expect[0]), abs(int(peak[1]) - expect[1])) <= 1, peak


def test_north_up_identity_and_rotated_corners():
    data = np.random.default_rng(1).normal(size=(21, 21))
    # Sub-1e-6 differences come from the spherical offsets (the TAN plane is not flat).
    np.testing.assert_allclose(
        viz.north_up(data, _sky_wcs(0.0, crpix=(11.0, 11.0))), data, atol=1e-6
    )
    np.testing.assert_allclose(
        viz.north_up(data, _sky_wcs(90.0, crpix=(11.0, 11.0))), np.rot90(data, 1), atol=1e-6
    )
    out = viz.north_up(data, _sky_wcs(45.0, crpix=(11.0, 11.0)))
    assert np.isnan(out[0, 0]) and np.isnan(out[-1, -1])  # uncovered corners are no data
    assert np.isfinite(out[10, :]).all()  # the central row is fully covered


def test_north_up_keeps_a_single_hot_pixel_and_nan_holes():
    data = np.zeros((31, 31))
    data[20, 7] = 5.0  # a one-pixel artefact must stay visible after rotation
    data[5, 25] = np.nan
    w = _sky_wcs(210.0, flip=True)
    out = viz.north_up(data, w)
    assert np.nanmax(out) > 1.0
    assert np.isfinite(out[15, 15])
    clean = data.copy()
    clean[5, 25] = 0.0
    # the planted hole adds NaN pixels beyond the uncovered corners
    assert np.isnan(out).sum() > np.isnan(viz.north_up(clean, w)).sum()


def test_north_up_keeps_interior_nan_visible():
    data = np.zeros((31, 31))
    data[12, 14] = np.nan
    out = viz.north_up(data, _sky_wcs(37.0))
    assert np.isnan(out[8:23, 8:23]).any()  # not filled with interpolated sky


def test_north_up_uses_icrs_north_for_galactic_wcs():
    from astropy.coordinates import SkyCoord
    from astropy.wcs import WCS

    w = WCS(naxis=2)
    w.wcs.ctype = ["GLON-TAN", "GLAT-TAN"]
    w.wcs.crval = [30.0, 10.0]
    w.wcs.crpix = [16.0, 16.0]
    w.wcs.cdelt = [-0.03 / 3600, 0.03 / 3600]
    shape = (31, 31)
    c = SkyCoord(w.pixel_to_world(15, 15)).icrs
    x, y = w.world_to_pixel(c.spherical_offsets_by(0 * u.arcsec, 0.27 * u.arcsec))
    yy, xx = np.mgrid[0:31, 0:31]
    data = np.exp(-((xx - float(x)) ** 2 + (yy - float(y)) ** 2) / 4.0)
    peak = np.unravel_index(np.nanargmax(viz.north_up(data, w)), shape)
    assert abs(int(peak[0]) - 24) <= 1 and abs(int(peak[1]) - 15) <= 1, peak


def test_north_up_accepts_header_and_rejects_missing_wcs():
    w = _sky_wcs(30.0)
    data = np.zeros((31, 31))
    np.testing.assert_array_equal(
        viz.north_up(data, w.to_header()), viz.north_up(data, w), strict=False
    )
    with pytest.raises(ValueError, match="celestial"):
        viz.north_up(data, fits.Header())


def _sky_cutout(path, theta, *, uid, band="F200W", amp=1.0):
    w = _sky_wcs(theta)
    yy, xx = np.mgrid[0:31, 0:31]
    data = amp * np.exp(-((xx - 15) ** 2 + (yy - 15) ** 2) / 8.0)
    data[_offset_pixel(w, (31, 31), 0.0, 0.27)] = 3.0 * amp  # a "neighbour" due north
    head = fits.Header({"SRCUID": uid, "BAND": band, "QUALFLAG": "ok"})
    fits.HDUList(
        [
            fits.PrimaryHDU(header=head),
            fits.ImageHDU(data.astype(np.float32), header=w.to_header(), name="SCI"),
        ]
    ).writeto(path)
    return path


def test_contact_sheet_north_up_and_shared_scales(tmp_path, monkeypatch):
    rows = []
    for i, (theta, amp) in enumerate(((0.0, 1.0), (60.0, 10.0), (200.0, 0.1))):
        p = _sky_cutout(tmp_path / f"s{i}.fits", theta, uid=f"s{i}", amp=amp)
        rows.append((f"s{i}", "F200W", p.as_posix(), 0.0, False, "ok"))
    nowcs = _cutout(tmp_path / "nowcs.fits", uid="nw")
    rows.append(("nw", "F200W", nowcs.as_posix(), 0.0, False, "ok"))
    t = Table(
        rows=rows, names=("source_uid", "band", "path", "frac_nan", "on_edge", "quality_flag")
    )

    drawn = []
    real = viz._draw_panel

    def spy(ax, data, lines, flag, norm, wrap):
        drawn.append((data, list(lines), norm))
        real(ax, data, lines, flag, norm, wrap)

    monkeypatch.setattr(viz, "_draw_panel", spy)
    viz.contact_sheet(t, tmp_path / "sheet.png", north_up=True, scale="all")
    assert len(drawn) == 4
    for data, _, _ in drawn[:3]:  # every rotated panel shows the neighbour straight above
        peak = np.unravel_index(np.nanargmax(data), data.shape)
        assert abs(int(peak[0]) - 24) <= 1 and abs(int(peak[1]) - 15) <= 1, peak
    assert "no WCS" in drawn[3][1][-1]  # unrotatable panel is labeled, not silently wrong
    assert all(norm is drawn[0][2] for _, _, norm in drawn)  # one shared interval

    drawn.clear()
    viz.contact_sheet(t, tmp_path / "sheet_panel.png")
    assert len({id(norm) for _, _, norm in drawn}) == 4  # default: one per panel
    with pytest.raises(ValueError, match="scale"):
        viz.contact_sheet(t, tmp_path / "x.png", scale="row")


def test_contact_sheet_band_scale_shares_within_band(tmp_path, monkeypatch):
    t = _cutout_table(tmp_path, bands=("F150W", "F444W"))
    norms = []
    real = viz._draw_panel

    def spy(ax, data, lines, flag, norm, wrap):
        norms.append((lines[1].split()[0], norm))
        real(ax, data, lines, flag, norm, wrap)

    monkeypatch.setattr(viz, "_draw_panel", spy)
    viz.contact_sheet(t, tmp_path / "band.png", scale="band")
    by_band = {}
    for band, norm in norms:
        if norm is not None:
            by_band.setdefault(band.rstrip(":"), set()).add(id(norm))
    assert sorted(by_band) == ["F150W", "F444W"]
    assert all(len(ids) == 1 for ids in by_band.values())
    assert by_band["F150W"] != by_band["F444W"]
