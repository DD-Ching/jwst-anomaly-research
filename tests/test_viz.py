"""Tests for viz: stretch, single-cutout PNGs and contact sheets from synthetic cutouts."""

from __future__ import annotations

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
