"""Tests for cached, checksummed acquisition. Offline tests fake the MAST download endpoint."""

from __future__ import annotations

import hashlib
import importlib.util
import io
import logging
import os
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import MaskedColumn, Table

from jwst_anomaly import acquire, paths, query, schema

REPO = Path(__file__).resolve().parents[1]
CAT = "jw02736-o001_t001_nircam_clear-f200w"
MIRI = "jw02736-o002_t001_miri_f770w"


def ecsv_bytes(versions: dict) -> bytes:
    t = Table({"label": [1, 2], "xcentroid": [1.5, 2.5]})
    t.meta["version"] = versions
    buf = io.StringIO()
    t.write(buf, format="ascii.ecsv")
    return buf.getvalue().encode()


def fits_bytes(**cards) -> bytes:
    hdu = fits.PrimaryHDU(data=np.zeros((2, 2), dtype=np.float32))
    hdu.header.update(cards)
    buf = io.BytesIO()
    hdu.writeto(buf)
    return buf.getvalue()


class FakeServer:
    """Stands in for MAST: serves registered files and records every download request."""

    def __init__(self):
        self.files: dict[str, bytes] = {}
        self.calls: list[str] = []
        self.truncate: set[str] = set()
        self.missing: set[str] = set()

    def add(self, obs_id: str, fname: str, data: bytes, sub: str = "CAT") -> dict:
        uri = f"mast:JWST/product/{fname}"
        self.files[uri] = data
        return {
            "obsID": "1",
            "obs_id": obs_id,
            "productFilename": fname,
            "productSubGroupDescription": sub,
            "dataURI": uri,
            "size": len(data),
            "calib_level": 3,
            "prvversion": "2.0.1",
        }

    def download_file(self, uri, *, local_path=None, cache=True, verbose=True, **kwargs):
        self.calls.append(uri)
        assert cache is False and Path(local_path).suffix == ".part"
        if uri in self.missing:
            msg = "HTTPError: 404 Client Error: Not Found for url: https://mast.stsci.edu/..."
            return "ERROR", msg, "https://mast.stsci.edu/..."
        data = self.files[uri]
        Path(local_path).write_bytes(data[:-1] if uri in self.truncate else data)
        return "COMPLETE", None, None

    def served_size(self, uri):
        """Content-Length of the full file, as the real download service declares it."""
        return len(self.files[uri]) if uri in self.files else None


@pytest.fixture(autouse=True)
def no_mast_token(monkeypatch):
    """Never send a developer's real MAST token from these tests."""
    monkeypatch.delenv(query.TOKEN_ENV, raising=False)


@pytest.fixture
def server(monkeypatch):
    s = FakeServer()
    monkeypatch.setattr(query.Observations, "download_file", s.download_file)
    monkeypatch.setattr(acquire, "_served_size", s.served_size)
    return s


@pytest.fixture
def two_products(server):
    rows = [
        server.add(
            CAT,
            f"{CAT}_cat.ecsv",
            ecsv_bytes({"jwst": "2.0.1", "photutils": "2.3.0", "astropy": "7.2.0"}),
        ),
        server.add(
            MIRI, f"{MIRI}_i2d.fits", fits_bytes(CAL_VER="2.0.1", CRDS_CTX="jwst_1234.pmap"), "I2D"
        ),
    ]
    return Table(rows=[list(r.values()) for r in rows], names=list(rows[0]))


def test_fetch_downloads_verifies_and_writes_manifest(tmp_path, server, two_products):
    manifest = tmp_path / "m.ecsv"
    out = acquire.fetch_products(two_products, data_root=tmp_path, manifest_path=manifest)

    schema.validate(out, schema.MANIFEST_COLUMNS)
    assert out.meta["provenance"] == "observed"
    assert list(out["status"]) == ["downloaded", "downloaded"]
    assert list(out["dataURI"]) == sorted(two_products["dataURI"])
    for row in out:
        local = tmp_path / row["local_path"]
        assert (
            row["local_path"]
            == f"cache/mast/{row['productFilename'].rsplit('_', 1)[0]}/{row['productFilename']}"
        )
        assert local.read_bytes() == server.files[row["dataURI"]]
        assert row["sha256"] == hashlib.sha256(local.read_bytes()).hexdigest()
        assert row["size"] == local.stat().st_size
        assert row["retrieved_utc"].endswith("Z")
    versions = dict(zip(out["productFilename"], out["pipeline_version"], strict=True))
    assert versions[f"{CAT}_cat.ecsv"] == "jwst=2.0.1;photutils=2.3.0;astropy=7.2.0"
    assert versions[f"{MIRI}_i2d.fits"] == "jwst=2.0.1;crds=jwst_1234.pmap"

    text = manifest.read_bytes()
    assert b"\r\n" not in text
    saved = acquire.read_manifest(manifest)
    assert saved.colnames == list(schema.MANIFEST_COLUMNS)
    assert list(saved["sha256"]) == list(out["sha256"])
    assert not list(tmp_path.rglob("*.part"))


def test_second_fetch_is_a_noop(tmp_path, server, two_products):
    manifest = tmp_path / "m.ecsv"
    acquire.fetch_products(two_products, data_root=tmp_path, manifest_path=manifest)
    before = (manifest.read_bytes(), manifest.stat().st_mtime_ns)
    server.calls.clear()

    out = acquire.fetch_products(two_products, data_root=tmp_path, manifest_path=manifest)

    assert server.calls == []
    assert list(out["status"]) == ["cached", "cached"]
    assert (manifest.read_bytes(), manifest.stat().st_mtime_ns) == before


def test_corrupted_local_file_is_redownloaded(tmp_path, server, two_products):
    manifest = tmp_path / "m.ecsv"
    out = acquire.fetch_products(two_products, data_root=tmp_path, manifest_path=manifest)
    local = tmp_path / out["local_path"][0]
    good = local.read_bytes()
    local.write_bytes(b"X" * len(good))  # same size, different content
    server.calls.clear()

    again = acquire.fetch_products(two_products, data_root=tmp_path, manifest_path=manifest)

    assert server.calls == [out["dataURI"][0]]
    assert list(again["status"]) == ["downloaded", "cached"]
    assert local.read_bytes() == good


def test_size_mismatch_keeps_nothing_and_records_the_rest(tmp_path, server, two_products):
    manifest = tmp_path / "m.ecsv"
    bad = two_products["dataURI"][0]
    server.truncate.add(bad)

    with pytest.raises(acquire.DownloadError, match="MAST reports"):
        acquire.fetch_products(two_products, data_root=tmp_path, manifest_path=manifest)

    assert not (tmp_path / str(acquire.local_relpath(two_products[0]))).exists()
    assert not list(tmp_path.rglob("*.part"))
    assert list(acquire.read_manifest(manifest)["dataURI"]) == [two_products["dataURI"][1]]


def test_stale_mast_size_is_accepted_when_the_service_serves_it(
    tmp_path, server, two_products, caplog
):
    """A reprocessed product: MAST still lists the old size, the service declares the new one."""
    manifest = tmp_path / "m.ecsv"
    true_size = int(two_products["size"][0])
    two_products["size"][0] = true_size - 7  # stale listing
    with caplog.at_level(logging.WARNING, logger="jwst_anomaly.acquire"):
        out = acquire.fetch_products(two_products, data_root=tmp_path, manifest_path=manifest)
    assert list(out["status"]) == ["downloaded", "downloaded"]
    assert out["size"][0] == true_size
    assert "listing is stale" in caplog.text

    server.calls.clear()  # the verified copy is kept, not re-downloaded on every run
    again = acquire.fetch_products(two_products, data_root=tmp_path, manifest_path=manifest)
    assert server.calls == []
    assert list(again["status"]) == ["cached", "cached"]

    local = tmp_path / out["local_path"][0]
    local.write_bytes(b"Z" * true_size)  # verified size, wrong content: re-downloaded
    again = acquire.fetch_products(two_products, data_root=tmp_path, manifest_path=manifest)
    assert server.calls == [out["dataURI"][0]]
    assert local.read_bytes() == server.files[out["dataURI"][0]]


def test_stale_mast_size_and_short_transfer_still_fail(tmp_path, server, two_products):
    bad = two_products["dataURI"][0]
    two_products["size"][0] += 3  # stale listing, and the transfer is also incomplete
    server.truncate.add(bad)
    with pytest.raises(acquire.DownloadError, match="download service"):
        acquire.fetch_products(two_products, data_root=tmp_path, manifest_path=tmp_path / "m.ecsv")
    assert not (tmp_path / str(acquire.local_relpath(two_products[0]))).exists()


def test_served_size_reads_content_length(monkeypatch):
    seen = {}

    class Response:
        headers = {"Content-Length": "3405952"}

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    def fake_urlopen(request, timeout):
        seen.update(url=request.full_url, method=request.get_method())
        return Response()

    monkeypatch.setattr(acquire.urllib.request, "urlopen", fake_urlopen)
    assert acquire._served_size(f"mast:JWST/product/{CAT}_cat.ecsv") == 3405952
    assert seen["method"] == "HEAD"
    assert seen["url"].endswith(f"/api/v0.1/Download/file?uri=mast:JWST/product/{CAT}_cat.ecsv")

    def failing_urlopen(request, timeout):
        raise acquire.urllib.error.URLError("offline")

    monkeypatch.setattr(acquire.urllib.request, "urlopen", failing_urlopen)
    assert acquire._served_size("mast:JWST/product/x_cat.ecsv") is None


def test_404_mentions_exclusive_access(tmp_path, server, two_products):
    server.missing.add(two_products["dataURI"][0])
    with pytest.raises(acquire.DownloadError, match="EXCLUSIVE_ACCESS"):
        acquire.fetch_products(two_products, data_root=tmp_path, manifest_path=tmp_path / "m.ecsv")


def test_oversized_product_is_not_downloaded(tmp_path, server, two_products):
    with pytest.raises(acquire.DownloadError, match="max_size_bytes"):
        acquire.fetch_products(
            two_products, data_root=tmp_path, manifest_path=tmp_path / "m.ecsv", max_size_bytes=10
        )
    assert server.calls == []


def test_unknown_size_is_refused_not_adopted(tmp_path, server, two_products):
    manifest = tmp_path / "m.ecsv"
    acquire.fetch_products(two_products, data_root=tmp_path, manifest_path=manifest)
    manifest.unlink()  # file stays on disk, but without a manifest row it is unverifiable
    two_products["size"] = MaskedColumn(two_products["size"], mask=[True, False])
    server.calls.clear()

    with pytest.raises(acquire.DownloadError, match="no size"):
        acquire.fetch_products(two_products, data_root=tmp_path, manifest_path=manifest)

    assert server.calls == []
    assert list(acquire.read_manifest(manifest)["dataURI"]) == [two_products["dataURI"][1]]


def test_cached_file_with_newer_mast_version_warns(tmp_path, server, two_products, caplog):
    manifest = tmp_path / "m.ecsv"
    acquire.fetch_products(two_products, data_root=tmp_path, manifest_path=manifest)
    two_products["prvversion"] = ["3.0.0", "2.0.1"]  # MAST reprocessed the first one
    with caplog.at_level(logging.WARNING, logger="jwst_anomaly.acquire"):
        out = acquire.fetch_products(two_products, data_root=tmp_path, manifest_path=manifest)
    assert list(out["status"]) == ["cached", "cached"]
    assert "prvversion=3.0.0" in caplog.text and f"{CAT}_cat.ecsv" in caplog.text


def test_file_on_disk_without_manifest_row_is_adopted(tmp_path, server, two_products):
    manifest = tmp_path / "m.ecsv"
    acquire.fetch_products(two_products, data_root=tmp_path, manifest_path=manifest)
    manifest.unlink()
    server.calls.clear()

    out = acquire.fetch_products(two_products, data_root=tmp_path, manifest_path=manifest)

    assert server.calls == []
    assert list(out["status"]) == ["adopted", "adopted"]
    assert len(acquire.read_manifest(manifest)) == 2


def test_manifest_keeps_rows_of_other_products(tmp_path, server, two_products):
    manifest = tmp_path / "m.ecsv"
    acquire.fetch_products(two_products[:1], data_root=tmp_path, manifest_path=manifest)
    out = acquire.fetch_products(two_products[1:], data_root=tmp_path, manifest_path=manifest)
    assert len(out) == 1
    assert list(acquire.read_manifest(manifest)["dataURI"]) == sorted(two_products["dataURI"])


def test_defaults_use_data_root_and_manifests_dir(tmp_path, monkeypatch, server, two_products):
    monkeypatch.setenv(paths.DATA_ENV, str(tmp_path / "data"))
    monkeypatch.setattr(paths, "manifests_dir", lambda: tmp_path / "manifests")
    out = acquire.fetch_products(two_products)
    assert (tmp_path / "data" / out["local_path"][0]).is_file()
    assert (tmp_path / "manifests" / acquire.DEFAULT_MANIFEST_NAME).is_file()


@pytest.mark.parametrize("fname", ["../evil_cat.ecsv", "a/b_cat.ecsv", "c:evil.fits", ""])
def test_unsafe_names_are_rejected_before_any_download(tmp_path, server, two_products, fname):
    two_products["productFilename"] = two_products["productFilename"].astype("U80")
    two_products["productFilename"][1] = fname
    with pytest.raises(ValueError, match="unsafe"):
        acquire.fetch_products(two_products, data_root=tmp_path, manifest_path=tmp_path / "m.ecsv")
    assert server.calls == []


def test_read_pipeline_version(tmp_path):
    ecsv = tmp_path / "x_cat.ecsv"
    ecsv.write_bytes(ecsv_bytes({"jwst": "3.0.0", "numpy": "2.4.4"}))
    assert acquire.read_pipeline_version(ecsv) == "jwst=3.0.0"
    fit = tmp_path / "x.part"
    fit.write_bytes(fits_bytes(CAL_VER="1.2.3"))
    assert acquire.read_pipeline_version(fit, "x_i2d.fits") == "jwst=1.2.3"
    assert acquire.read_pipeline_version(fit, "x_asn.json") == ""
    broken = tmp_path / "broken_cat.ecsv"
    broken.write_text("not an ecsv file\n")
    assert acquire.read_pipeline_version(broken) == ""


def test_verify_manifest_reports_problems(tmp_path, server, two_products):
    manifest = tmp_path / "m.ecsv"
    out = acquire.fetch_products(two_products, data_root=tmp_path, manifest_path=manifest)
    assert set(acquire.verify_manifest(manifest, tmp_path)["check"]) == {"ok"}

    first, second = (tmp_path / p for p in out["local_path"])
    first.write_bytes(b"Y" * first.stat().st_size)
    second.unlink()
    checks = acquire.verify_manifest(manifest, tmp_path)
    assert list(checks["check"]) == ["sha256_mismatch", "missing"]
    assert checks.meta["provenance"] == "derived"
    only = acquire.verify_manifest(manifest, tmp_path, uris=[out["dataURI"][1]])
    assert list(only["check"]) == ["missing"]


def test_read_manifest_keeps_provenance(tmp_path, server, two_products):
    manifest = tmp_path / "m.ecsv"
    acquire.fetch_products(two_products, data_root=tmp_path, manifest_path=manifest)
    saved = acquire.read_manifest(manifest)
    assert schema.validate(saved, schema.MANIFEST_COLUMNS) is saved
    assert saved.meta["source"] == acquire.MANIFEST_SOURCE


def test_fetch_script_end_to_end_offline(tmp_path, monkeypatch, server, capsys):
    """The fetch script against a fake MAST: catalogs only, then a second no-op run."""
    obs = Table(
        {
            "obsid": ["1", "2"],
            "obs_id": [CAT, MIRI],
            "proposal_id": ["2736", "2736"],
            "instrument_name": ["NIRCAM/IMAGE", "MIRI/IMAGE"],
            "filters": ["F200W", "F770W"],
            "calib_level": [3, 3],
            "t_exptime": [1.0, 1.0],
            "s_ra": [110.8, 110.8],
            "s_dec": [-73.5, -73.5],
            "dataRights": ["PUBLIC", "PUBLIC"],
            "target_name": ["SMACS", "SMACS"],
        }
    )
    rows = [
        server.add(CAT, f"{CAT}_cat.ecsv", ecsv_bytes({"jwst": "2.0.1"})),
        server.add(CAT, f"{CAT}_i2d.fits", fits_bytes(CAL_VER="2.0.1"), "I2D"),
        server.add(MIRI, f"{MIRI}_cat.ecsv", ecsv_bytes({"jwst": "2.0.1"})),
    ]
    for r, obsid in zip(rows, ["1", "1", "2"], strict=True):
        r.update(obsID=obsid, dataRights="PUBLIC")
    products = Table(rows=[list(r.values()) for r in rows], names=list(rows[0]))
    criteria = []
    monkeypatch.setattr(
        query.Observations, "query_criteria", lambda **c: criteria.append(c) or obs.copy()
    )
    monkeypatch.setattr(query.Observations, "get_product_list", lambda ids: products.copy())
    monkeypatch.setattr(
        query, "mast_relative_path", lambda uris, verbose=True: [f"jwst/public/{u}" for u in uris]
    )
    config = tmp_path / "sample.yaml"
    config.write_text(
        "archive: {collection: JWST, calib_level: 3, data_rights: PUBLIC, "
        "product_subgroups: [CAT, I2D]}\n"
        f"samples:\n  - id: nircam\n    obs_ids: [{CAT}]\n  - id: miri\n    obs_ids: [{MIRI}]\n"
    )
    script = _load_script()
    argv = [
        "--config",
        str(config),
        "--catalogs-only",
        "--data-root",
        str(tmp_path / "d"),
        "--manifest-dir",
        str(tmp_path / "manifests"),
    ]

    assert script.main(argv) == 0
    first = capsys.readouterr().out
    assert "summary: downloaded 2, adopted 0, cached 0" in first
    assert "verified 2/2 files" in first
    assert criteria[0]["obs_id"] == sorted([CAT, MIRI])
    manifest = tmp_path / "manifests" / "sample.ecsv"
    index = Table.read(tmp_path / "manifests" / "sample_products.ecsv", format="ascii.ecsv")
    assert len(acquire.read_manifest(manifest)) == 2
    assert sorted(index["productSubGroupDescription"]) == ["CAT", "CAT", "I2D"]  # i2d listed too
    assert set(index["sample_id"]) == {"nircam", "miri"}
    assert all(u.startswith("s3://stpubdata/") for u in index["cloud_uri"])
    snapshot = {p.name: p.read_bytes() for p in (tmp_path / "manifests").iterdir()}

    server.calls.clear()
    assert script.main(argv) == 0
    assert "summary: downloaded 0, adopted 0, cached 2" in capsys.readouterr().out
    assert server.calls == []
    assert {p.name: p.read_bytes() for p in (tmp_path / "manifests").iterdir()} == snapshot

    assert script.main([*argv, "--sample", "nope"]) == 2
    verify = ["--config", str(config), "--verify-only", "--data-root", str(tmp_path / "d")]
    assert script.main([*verify, "--manifest-dir", str(tmp_path / "manifests")]) == 0
    assert script.main([*verify, "--manifest-dir", str(tmp_path / "nowhere")]) == 1

    # MAST stops listing the NIRCam i2d: re-listing that observation drops it from the index,
    # while rows of observations not re-listed this time (MIRI) are kept.
    products = products[products["productSubGroupDescription"] != "I2D"]  # read by the fake
    assert script.main([*argv, "--sample", "nircam"]) == 0
    index = Table.read(tmp_path / "manifests" / "sample_products.ecsv", format="ascii.ecsv")
    assert sorted(index["productFilename"]) == [f"{CAT}_cat.ecsv", f"{MIRI}_cat.ecsv"]


def _load_script():
    path = REPO / "scripts" / "fetch_reference_sample.py"
    spec = importlib.util.spec_from_file_location("fetch_reference_sample", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# --- live MAST ---------------------------------------------------------------------------


@pytest.mark.network
def test_live_fetch_small_catalog_is_idempotent(tmp_path, monkeypatch):
    obs = query.query_observations(obs_id="jw02736-o002_t001_miri_f1500w")
    prods = query.list_products(obs, subgroups="CAT")
    assert len(prods) == 1 and prods["size"][0] < 1_000_000  # ~120 kB on 2026-10-07
    manifest = tmp_path / "m.ecsv"

    out = acquire.fetch_products(prods, data_root=tmp_path, manifest_path=manifest)
    row = out[0]
    assert row["status"] == "downloaded"
    assert row["size"] == prods["size"][0]
    assert row["sha256"] == acquire.sha256_file(tmp_path / row["local_path"])
    assert row["pipeline_version"].startswith("jwst=")

    calls = []
    real = query.Observations.download_file
    monkeypatch.setattr(
        query.Observations, "download_file", lambda *a, **k: calls.append(a) or real(*a, **k)
    )
    mtime = os.stat(manifest).st_mtime_ns
    again = acquire.fetch_products(prods, data_root=tmp_path, manifest_path=manifest)
    assert calls == [] and again["status"][0] == "cached"
    assert os.stat(manifest).st_mtime_ns == mtime
