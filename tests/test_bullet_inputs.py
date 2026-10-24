"""Offline checks of the Bullet Cluster input-check parsers and the committed image list."""

import sys
from pathlib import Path

import pytest
from astropy.io import fits
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import bullet_inputs as bi  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / bi.OUT

GOLD = """K1b.1 104.63017980829942 -55.94343598701029 0.3 0.3 0 3.24 25
K1b.2 104.63167480860355 -55.942185987037675 0.3 0.3 0 3.24 25
K12.1 104.6 -55.95 0.3 0.3 0 0.897 25
"""
FULL = """#ID RA(ICRS) DEC(ICRS) specz quality
K1b.1 104.63017980829942 -55.94343598701029 3.24 gold
K1b.2 104.63167480860355 -55.942185987037675 3.24 gold
K12.1 104.6 -55.95 0.897 gold
K40.1 104.61 -55.96 2.0 bronze
"""


def test_split_id():
    assert bi.split_id("K1b.2") == ("K1b", "K1")
    assert bi.split_id("K12.1") == ("K12", "K12")
    with pytest.raises(ValueError):
        bi.split_id("K1b")


def test_parse_lenstool_images_and_counts():
    t = bi.parse_lenstool_images(GOLD)
    assert list(t["image_id"]) == ["K1b.1", "K1b.2", "K12.1"]
    assert t["z_spec"].all()
    c = bi.counts(t)
    assert (c["images"], c["systems"], c["galaxies"], c["distinct_z"]) == (3, 2, 2, 2)
    with pytest.raises(ValueError):
        bi.parse_lenstool_images("K1b.1 104.6 -55.9 0.3 0.3 0 3.24\n")
    with pytest.raises(ValueError):
        bi.parse_lenstool_images("K1b.1 104.6 -55.9 0.3 0.3 0 0 25\n")


def test_parse_multim_full_and_compare():
    full = bi.parse_multim_full(FULL)
    assert list(full["grade"]) == ["gold", "gold", "gold", "bronze"]
    r = bi.gold_vs_full(bi.parse_lenstool_images(GOLD), full)
    assert (r["n_shared"], r["max_sep_arcsec"], r["max_abs_dz"]) == (3, 0.0, 0.0)
    with pytest.raises(ValueError):
        bi.parse_multim_full("K1b.1 104.6 -55.9 3.24 platinum\n")


def _blob(cols, unit="nJy"):
    t = fits.BinTableHDU.from_columns(
        [fits.Column(name=c, format="E", unit=unit if c.startswith("FLUX") else "") for c in cols]
    )
    return (fits.PrimaryHDU().header.tostring() + t.header.tostring()).encode()


def test_header_check():
    assert bi.header_check(_blob(bi.MAPPING_COLS))["passes"]
    r = bi.header_check(_blob(bi.MAPPING_COLS[1:]))
    assert r["mapping_columns_missing"] == ["RA"] and not r["passes"]
    assert not bi.header_check(_blob(bi.MAPPING_COLS, unit="uJy"))["passes"]


def test_committed_image_list():
    t = Table.read(OUT / bi.IMAGES)
    c = bi.counts(t)
    assert (c["images"], c["systems"], c["galaxies"]) == (135, 49, 27)
    assert t.meta["provenance"] == "observed"
