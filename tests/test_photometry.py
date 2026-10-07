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
    cat = photometry.load_dja_catalog(_dja(tmp_path / "x_phot.fits"), ["f200w", "f277w"])
    assert cat.meta["aperture_diameter_arcsec"] == 0.5
    assert cat.meta["provenance"] == "derived"
    np.testing.assert_allclose(cat["f200w_mag"][:2], [23.9, 23.9 - 2.5])
    assert np.isnan(cat["f200w_mag"][2])  # non-positive flux
    assert np.isnan(cat["f200w_mag"][3])  # APER_TRUNC
    assert cat["f200w_mag_err"][0] == pytest.approx(2.5 * np.log10(1.1))  # pipeline convention


def test_load_dja_catalog_rejects_unexpected_units(tmp_path):
    with pytest.raises(ValueError, match="expected uJy"):
        photometry.load_dja_catalog(_dja(tmp_path / "x.fits", unit="Jy"), ["f200w"])


def test_join_matched_photometry(tmp_path):
    cat = photometry.load_dja_catalog(_dja(tmp_path / "x_phot.fits"), ["f200w", "f277w"])
    src = _sources([10.0, 10.001 + 0.05 / 3600 / np.cos(np.deg2rad(5)), 10.1])
    src.remove_columns(["f444w_detected"])
    out = photometry.join_matched_photometry(src, cat, "dja05", radius_arcsec=0.2)
    schema.validate(out, schema.SOURCE_COLUMNS)
    assert out["f200w_dja05_abmag"][0] == pytest.approx(23.9)
    assert out["f277w_dja05_abmag"][1] == pytest.approx(23.9 - 2.5)
    assert np.isnan(out["f200w_dja05_abmag"][2])  # nothing within 0.2"
    assert out["dja05_match_sep_arcsec"][1] == pytest.approx(0.05, abs=0.01)
    assert out.meta["matched_photometry"]["n_matched"] == 2
    assert "f200w_dja05_abmag" not in src.colnames  # input untouched


def test_join_rejects_bad_label_and_disjoint_bands(tmp_path):
    cat = photometry.load_dja_catalog(_dja(tmp_path / "x_phot.fits"), ["f200w"])
    with pytest.raises(ValueError, match="alphanumeric"):
        photometry.join_matched_photometry(_sources([10.0]), cat, "dja_05")
    other = _sources([10.0])
    other.remove_columns(["f200w_detected", "f277w_detected"])
    with pytest.raises(ValueError, match="no photometry for sample bands"):
        other["f150w_detected"] = [True]
        photometry.join_matched_photometry(other, cat, "dja05")


class _Resp:
    def __init__(self, payload, length=None):
        self.payload, self.done = payload, False
        self.headers = {} if length is None else {"Content-Length": str(length)}

    def read(self, n=-1):
        if self.done:
            return b""
        self.done = True
        return self.payload

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _serve(monkeypatch, payloads):
    calls = []

    def fake_urlopen(url, timeout=0):
        calls.append(url)
        return _Resp(payloads[min(len(calls), len(payloads)) - 1])

    monkeypatch.setattr(photometry.urllib.request, "urlopen", fake_urlopen)
    return calls


def test_fetch_catalog_verifies_before_caching(tmp_path, monkeypatch):
    payload = b"fits bytes"
    good = hashlib.sha256(payload).hexdigest()
    calls = _serve(monkeypatch, [b"fits by", payload])  # first download truncated
    with pytest.raises(ValueError, match="truncated or changed"):
        photometry.fetch_catalog("https://e.org/a/x_phot.fits", good, cache_dir=tmp_path)
    assert list(tmp_path.iterdir()) == []  # nothing cached, no .part left
    path = photometry.fetch_catalog("https://e.org/a/x_phot.fits", good, cache_dir=tmp_path)
    assert path.read_bytes() == payload and path.name == f"{good[:12]}_x_phot.fits"
    photometry.fetch_catalog("https://e.org/a/x_phot.fits", good, cache_dir=tmp_path)
    assert len(calls) == 2  # third call served from the cache


def test_fetch_catalog_replaces_a_corrupt_cached_file(tmp_path, monkeypatch):
    payload = b"fits bytes"
    good = hashlib.sha256(payload).hexdigest()
    (tmp_path / f"{good[:12]}_x.fits").write_bytes(b"corrupt")
    calls = _serve(monkeypatch, [payload])
    path = photometry.fetch_catalog("https://e.org/x.fits", good, cache_dir=tmp_path)
    assert path.read_bytes() == payload and len(calls) == 1


def test_fetch_catalog_size_limit(tmp_path, monkeypatch):
    monkeypatch.setattr(
        photometry.urllib.request, "urlopen", lambda url, timeout=0: _Resp(b"x" * 10, length=10)
    )
    with pytest.raises(ValueError, match="max_bytes"):
        photometry.fetch_catalog("https://e.org/x.fits", "0" * 64, cache_dir=tmp_path, max_bytes=5)


def test_load_dja_catalog_requires_bands_and_is_derived(tmp_path):
    path = _dja(tmp_path / "x_phot.fits")
    assert photometry.load_dja_catalog(path, ["f200w"]).meta["provenance"] == "derived"
    with pytest.raises(ValueError, match="f444w"):
        photometry.load_dja_catalog(path, ["f200w", "f444w"])


def test_join_is_one_to_one(tmp_path):
    cat = photometry.load_dja_catalog(_dja(tmp_path / "x_phot.fits"), ["f200w", "f277w"])
    near = 0.05 / 3600 / np.cos(np.deg2rad(5.0))
    src = _sources([10.0, 10.0 + near])  # two fragments around DJA object 1
    src.remove_columns(["f444w_detected"])
    out = photometry.join_matched_photometry(src, cat, "dja05", radius_arcsec=0.2)
    assert np.isfinite(out["f200w_dja05_abmag"][0]) and np.isnan(out["f200w_dja05_abmag"][1])
    assert out.meta["matched_photometry"]["n_contested"] >= 1


def test_join_rejects_reserved_or_clashing_labels(tmp_path):
    cat = photometry.load_dja_catalog(_dja(tmp_path / "x_phot.fits"), ["f200w", "f277w"])
    src = _sources([10.0])
    src.remove_columns(["f444w_detected"])
    with pytest.raises(ValueError, match="alphanumeric"):
        photometry.join_matched_photometry(src, cat, "aper50")
    src["f200w_dja05_abmag"] = [20.0]
    with pytest.raises(ValueError, match="overwrite"):
        photometry.join_matched_photometry(src, cat, "dja05")


def test_joined_columns_keeps_only_the_join(tmp_path):
    cat = photometry.load_dja_catalog(_dja(tmp_path / "x_phot.fits"), ["f200w", "f277w"])
    src = _sources([10.0])
    src.remove_columns(["f444w_detected"])
    out = photometry.joined_columns(photometry.join_matched_photometry(src, cat, "dja05"), "dja05")
    assert out.colnames[0] == "source_uid"
    assert set(out.colnames[1:]) == {
        "dja05_match_sep_arcsec",
        "f200w_dja05_abmag",
        "f200w_dja05_abmag_err",
        "f277w_dja05_abmag",
        "f277w_dja05_abmag_err",
    }
    assert out.meta["provenance"] == "derived"


@pytest.mark.network
def test_live_dja_smacs0723_header_matches_loader_assumptions():
    """Read only the FITS headers of the pinned DJA file (HTTP range) and check the format."""
    import io
    import urllib.request

    from astropy.io import fits

    url = "https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/smacs0723-grizli-v7.4-fix_phot.fits"
    req = urllib.request.Request(url, headers={"Range": "bytes=0-399999"})
    with urllib.request.urlopen(req, timeout=60) as resp:  # noqa: S310
        head = resp.read()
    with fits.open(io.BytesIO(head), lazy_load_hdus=True, ignore_missing_end=True) as hdul:
        hdr = hdul[1].header
    names = {hdr[k] for k in hdr if k.startswith("TTYPE")}
    units = {hdr[f"TTYPE{i}"]: hdr.get(f"TUNIT{i}") for i in range(1, hdr["TFIELDS"] + 1)}
    for band in ("f090w", "f150w", "f200w", "f277w", "f356w", "f444w"):
        assert {f"{band}_flux_aper_1", f"{band}_fluxerr_aper_1", f"{band}_flag_aper_1"} <= names
        assert units[f"{band}_flux_aper_1"] == "uJy"
    assert abs(float(hdr["ASEC_1"]) - 0.5) < 1e-3
