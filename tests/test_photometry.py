"""Tests for photometry (D-013): DJA loader, positional join, checksum-verified fetch."""

from __future__ import annotations

import hashlib

import numpy as np
import pytest
from astropy.table import Table

from jwst_anomaly import photometry, schema


def _dja(path, *, unit="uJy"):
    t = Table(
        {
            "id": [1, 2, 3, 4],
            "ra": [10.0, 10.001, 10.002, 10.003],
            "dec": [-5.0, -5.0, -5.0, -5.0],
            "f200w_flux_aper_1": [1.0, 10.0, -1.0, 5.0],
            "f200w_fluxerr_aper_1": [0.1, 0.5, 0.1, 0.5],
            "f200w_flag_aper_1": [0, 0x20, 0, 0x10],  # HASMASKED kept, TRUNC dropped
            "f277w_flux_aper_1": [2.0, 10.0, 1.0, 5.0],
            "f277w_fluxerr_aper_1": [0.2, 1.0, 0.1, 0.5],
            "f277w_flag_aper_1": [0, 0, 0, 0],
        }
    )
    for c in t.colnames:
        if "flux" in c:
            t[c].unit = unit
    t.meta["ASEC_1"] = 0.5
    t.write(path, overwrite=True)
    return path


def _sources(ras):
    t = Table(
        {
            "source_uid": [f"s{i}" for i in range(len(ras))],
            "ra": ras,
            "dec": [-5.0] * len(ras),
            "ref_band": ["F200W"] * len(ras),
            "n_bands": [2] * len(ras),
            "f200w_detected": [True] * len(ras),
            "f277w_detected": [True] * len(ras),
            "f444w_detected": [True] * len(ras),
        }
    )
    t.meta.update(provenance="derived", source="test sources")
    return t


def test_load_dja_catalog_magnitudes_and_flags(tmp_path):
    cat = photometry.load_dja_catalog(_dja(tmp_path / "x_phot.fits"), ["f200w", "f277w", "f444w"])
    assert cat.meta["aperture_diameter_arcsec"] == 0.5
    assert cat.meta["provenance"] == "observed"
    assert "f444w_mag" not in cat.colnames  # absent from the catalog
    np.testing.assert_allclose(cat["f200w_mag"][:2], [23.9, 23.9 - 2.5])
    assert np.isnan(cat["f200w_mag"][2])  # non-positive flux
    assert np.isnan(cat["f200w_mag"][3])  # APER_TRUNC
    assert cat["f200w_mag_err"][0] == pytest.approx(2.5 / np.log(10) * 0.1)


def test_load_dja_catalog_rejects_unexpected_units(tmp_path):
    with pytest.raises(ValueError, match="expected uJy"):
        photometry.load_dja_catalog(_dja(tmp_path / "x.fits", unit="Jy"), ["f200w"])


def test_join_matched_photometry(tmp_path):
    cat = photometry.load_dja_catalog(_dja(tmp_path / "x_phot.fits"), ["f200w", "f277w"])
    src = _sources([10.0, 10.001 + 0.05 / 3600 / np.cos(np.deg2rad(5)), 10.1])
    out = photometry.join_matched_photometry(src, cat, "dja05", radius_arcsec=0.2)
    schema.validate(out, schema.SOURCE_COLUMNS)
    assert out["f200w_dja05_abmag"][0] == pytest.approx(23.9)
    assert out["f277w_dja05_abmag"][1] == pytest.approx(23.9 - 2.5)
    assert np.isnan(out["f200w_dja05_abmag"][2])  # nothing within 0.2"
    assert out["dja05_match_sep_arcsec"][1] == pytest.approx(0.05, abs=0.01)
    assert "f444w_dja05_abmag" not in out.colnames
    assert out.meta["matched_photometry"]["n_matched"] == 2
    assert "f200w_dja05_abmag" not in src.colnames  # input untouched


def test_join_rejects_bad_label_and_disjoint_bands(tmp_path):
    cat = photometry.load_dja_catalog(_dja(tmp_path / "x_phot.fits"), ["f200w"])
    with pytest.raises(ValueError, match="label"):
        photometry.join_matched_photometry(_sources([10.0]), cat, "dja_05")
    other = _sources([10.0])
    other.remove_columns(["f200w_detected", "f277w_detected"])
    with pytest.raises(ValueError, match="no common bands"):
        photometry.join_matched_photometry(other, cat, "dja05")


def test_fetch_catalog_verifies_checksum(tmp_path, monkeypatch):
    payload = b"fits bytes"
    good = hashlib.sha256(payload).hexdigest()

    class Resp:
        def __init__(self):
            self.done = False

        def read(self, n=-1):
            if self.done:
                return b""
            self.done = True
            return payload

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    calls = []

    def fake_urlopen(url, timeout=0):
        calls.append(url)
        return Resp()

    monkeypatch.setattr(photometry.urllib.request, "urlopen", fake_urlopen)
    path = photometry.fetch_catalog("https://e.org/a/x_phot.fits", good, cache_dir=tmp_path)
    assert path.read_bytes() == payload and path.name == "x_phot.fits"
    photometry.fetch_catalog("https://e.org/a/x_phot.fits", good, cache_dir=tmp_path)
    assert len(calls) == 1  # cached
    with pytest.raises(ValueError, match="sha256"):
        photometry.fetch_catalog("https://e.org/a/x_phot.fits", "0" * 64, cache_dir=tmp_path)
