"""Tests for the MAST query layer. Offline tests fake astroquery; live ones are marked network."""

from __future__ import annotations

import logging

import numpy as np
import pytest
from astropy.table import MaskedColumn, Table

from jwst_anomaly import query, schema

NIRCAM_F200W = "jw02736-o001_t001_nircam_clear-f200w"


@pytest.fixture(autouse=True)
def no_mast_token(monkeypatch):
    """Never send a developer's real MAST token from these tests."""
    monkeypatch.delenv(query.TOKEN_ENV, raising=False)


def fake_observations() -> Table:
    """Three observation rows shaped like MAST CAOM results (obsid as strings, masked cols)."""
    t = Table(
        {
            "obsid": ["3", "1", "2"],
            "obs_id": ["jw00001-o001_t001_nircam_clear-f150w", NIRCAM_F200W, "jw00002-o001_x"],
            "proposal_id": ["1", "2736", "2"],
            "instrument_name": ["NIRCAM/IMAGE"] * 3,
            "filters": ["F150W", "F200W", "F200W"],
            "calib_level": [3, 3, 3],
            "t_exptime": [1.0, 2.0, 3.0],
            "s_ra": [1.0, 2.0, 3.0],
            "s_dec": [-1.0, -2.0, -3.0],
            "target_name": ["a", "b", "c"],
        }
    )
    t["dataRights"] = MaskedColumn(["PUBLIC", "PUBLIC", "EXCLUSIVE_ACCESS"])
    return t


def product(obsid, obs_id, fname, sub, calib=3, rights="PUBLIC", size=100, parent=None):
    return {
        "obsID": obsid,
        "obs_id": obs_id,
        "parent_obsid": parent or obsid,
        "productFilename": fname,
        "productSubGroupDescription": sub,
        "dataURI": f"mast:JWST/product/{fname}",
        "size": size,
        "calib_level": calib,
        "dataRights": rights,
        "prvversion": "2.0.1",
    }


def fake_products() -> Table:
    f150 = "jw00001-o001_t001_nircam_clear-f150w"
    rows = [
        product("1", NIRCAM_F200W, f"{NIRCAM_F200W}_cat.ecsv", "CAT", size=3_444_394),
        product("1", NIRCAM_F200W, f"{NIRCAM_F200W}_i2d.fits", "I2D", size=1_762_162_560),
        product("1", NIRCAM_F200W, f"{NIRCAM_F200W}_cat.ecsv", "CAT", size=3_444_394),  # dup
        product("1", NIRCAM_F200W, f"{NIRCAM_F200W}_segm.fits", "SEGM"),
        product("1", NIRCAM_F200W, "jw02736-o001_20260802t150051_image3_00002_asn.json", "ASN"),
        # level-2 member listed with the level-3 observation: own obsID, parent = level-3
        product(
            "99",
            "jw02736001001_02105_00004_nrca4",
            "jw02736001001_02105_00004_nrca4_i2d.fits",
            "I2D",
            calib=2,
            parent="1",
        ),
        # defensive: level-3 row under the right obsID but named after another observation
        product("1", NIRCAM_F200W, "jw09999-o001_t001_nircam_clear-f200w_cat.ecsv", "CAT"),
        # F150W2 shares the "clear-f150w" prefix; only obs_id + "_" may match
        product("3", f150, f"{f150}2_cat.ecsv", "CAT"),
        product("3", f150, f"{f150}_cat.ecsv", "CAT"),
        product("3", f150, f"{f150}_i2d.fits", "I2D", rights="EXCLUSIVE_ACCESS"),
    ]
    t = Table(rows=[list(r.values()) for r in rows], names=list(rows[0]))
    t["productSubGroupDescription"] = MaskedColumn(t["productSubGroupDescription"])
    return t


@pytest.fixture
def fake_mast(monkeypatch):
    """Replace the network methods of the shared Observations client; record the calls."""
    calls: dict[str, list] = {"query_criteria": [], "get_product_list": [], "path_lookup": []}

    def query_criteria(**criteria):
        calls["query_criteria"].append(criteria)
        return fake_observations()

    def get_product_list(obsids):
        calls["get_product_list"].append(obsids)
        return fake_products()

    def path_lookup(uris, verbose=True):
        calls["path_lookup"].append(list(uris))
        return [None if "segm" in u else f"jwst/public/{u.rsplit('/', 1)[1]}" for u in uris]

    monkeypatch.setattr(query.Observations, "query_criteria", query_criteria)
    monkeypatch.setattr(query.Observations, "get_product_list", get_product_list)
    monkeypatch.setattr(query, "mast_relative_path", path_lookup)
    monkeypatch.delenv(query.TOKEN_ENV, raising=False)
    return calls


def test_query_observations_defaults_to_public_jwst(fake_mast):
    obs = query.query_observations(proposal_id="2736", calib_level=3)

    sent = fake_mast["query_criteria"][0]
    assert sent["obs_collection"] == "JWST"
    assert sent["dataRights"] == "PUBLIC"
    assert sent["proposal_id"] == "2736"
    # the EXCLUSIVE_ACCESS row the fake "server" returned anyway is dropped
    assert set(obs["dataRights"]) == {"PUBLIC"}
    assert list(obs["obs_id"]) == sorted(obs["obs_id"])
    assert obs.meta["provenance"] == schema.Provenance.OBSERVED.value
    assert "proposal_id='2736'" in obs.meta["source"]
    schema.validate(obs, schema.OBSERVATION_COLUMNS)


def test_query_observations_public_only_false_keeps_exclusive(fake_mast):
    obs = query.query_observations(public_only=False, proposal_id="2736")
    assert "dataRights" not in fake_mast["query_criteria"][0]
    assert "EXCLUSIVE_ACCESS" in set(obs["dataRights"])


@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({}, "constraint"),
        ({"obs_collection": "JWST", "dataRights": "PUBLIC"}, "constraint"),  # whole collection
        ({"proposal_id": "1", "dataRights": "EXCLUSIVE_ACCESS"}, "contradicts"),
    ],
)
def test_query_observations_rejects_bad_criteria(fake_mast, kwargs, match):
    with pytest.raises(ValueError, match=match):
        query.query_observations(**kwargs)
    assert fake_mast["query_criteria"] == []


def test_query_observations_empty_result_keeps_contract(monkeypatch):
    monkeypatch.setattr(query.Observations, "query_criteria", lambda **c: Table())
    obs = query.query_observations(proposal_id="0")
    assert len(obs) == 0
    schema.validate(obs, schema.OBSERVATION_COLUMNS)


def test_list_products_keeps_only_level3_products_of_the_observations(fake_mast):
    obs = query.query_observations(proposal_id="2736")
    prods = query.list_products(obs)

    assert fake_mast["get_product_list"] == [["1", "3"]]  # one batched call, PUBLIC obs only
    assert list(prods["productFilename"]) == [
        "jw00001-o001_t001_nircam_clear-f150w_cat.ecsv",
        f"{NIRCAM_F200W}_cat.ecsv",
        f"{NIRCAM_F200W}_i2d.fits",
    ]
    assert list(prods["dataURI"]) == sorted(prods["dataURI"])
    assert set(prods["calib_level"]) == {3}
    assert "cloud_uri" not in prods.colnames
    assert prods.meta["provenance"] == "observed"
    schema.validate(prods, schema.PRODUCT_COLUMNS)


def test_list_products_subgroups_and_rights(fake_mast):
    obs = query.query_observations(public_only=False, proposal_id="2736")
    cats = query.list_products(obs, subgroups="CAT")
    assert set(cats["productSubGroupDescription"]) == {"CAT"}
    everything = query.list_products(obs, subgroups=("CAT", "I2D"), public_only=False)
    assert "EXCLUSIVE_ACCESS" in set(everything["dataRights"])


def test_list_products_cloud_uris(fake_mast):
    obs = query.query_observations(public_only=False, proposal_id="2736")
    prods = query.list_products(
        obs, subgroups=("CAT", "I2D", "SEGM"), public_only=False, cloud_uris=True
    )
    uri = dict(zip(prods["productFilename"], prods["cloud_uri"], strict=True))
    assert uri[f"{NIRCAM_F200W}_cat.ecsv"] == f"s3://stpubdata/jwst/public/{NIRCAM_F200W}_cat.ecsv"
    assert uri[f"{NIRCAM_F200W}_segm.fits"] == ""  # MAST had no path
    assert uri["jw00001-o001_t001_nircam_clear-f150w_i2d.fits"] == ""  # not in the public bucket
    assert len(fake_mast["path_lookup"]) == 1  # batched


def test_resolve_cloud_uris_rejects_misaligned_lookup(monkeypatch):
    monkeypatch.setattr(query, "mast_relative_path", lambda uris, verbose=True: ["only/one"])
    with pytest.raises(RuntimeError, match="2 URIs"):
        query.resolve_cloud_uris(["mast:JWST/product/a", "mast:JWST/product/b"])


def test_str_values_handles_masked_numbers_and_strings():
    ints = MaskedColumn([87602476, 2], mask=[False, True])
    strs = MaskedColumn(["CAT", "x"], mask=[False, True])
    assert list(query.str_values(ints)) == ["87602476", ""]
    assert list(query.str_values(strs)) == ["CAT", ""]
    assert list(query.str_values(np.array([1, 2]))) == ["1", "2"]


def test_list_products_accepts_integer_obsids(fake_mast):
    obs = query.query_observations(proposal_id="2736")
    obs["obsid"] = MaskedColumn(np.asarray(obs["obsid"], dtype=int))  # older astroquery: ints
    assert len(query.list_products(obs)) == 3


def test_list_products_without_observations_makes_no_request(fake_mast):
    empty = fake_observations()[:0]
    prods = query.list_products(empty, cloud_uris=True)
    assert len(prods) == 0
    assert fake_mast["get_product_list"] == [] and fake_mast["path_lookup"] == []
    schema.validate(prods, schema.PRODUCT_COLUMNS)


def test_mast_client_logs_in_once_with_env_token(monkeypatch, caplog, capsys):
    secret = "not-a-real-token-123"
    logins = []
    monkeypatch.setattr(query, "_logged_in", False)
    monkeypatch.setattr(query.Observations, "login", lambda token=None: logins.append(token))
    monkeypatch.setenv(query.TOKEN_ENV, secret)
    with caplog.at_level(logging.DEBUG):
        query.mast_client()
        query.mast_client()
    assert logins == [secret]
    out = capsys.readouterr()
    assert secret not in caplog.text + out.out + out.err


def test_mast_client_without_token_does_not_log_in(monkeypatch):
    monkeypatch.setattr(query, "_logged_in", False)
    monkeypatch.setattr(query.Observations, "login", lambda **k: pytest.fail("login called"))
    monkeypatch.delenv(query.TOKEN_ENV, raising=False)
    assert query.mast_client() is query.Observations


# --- live MAST ---------------------------------------------------------------------------


@pytest.mark.network
def test_live_program_2736_level3_images():
    obs = query.query_observations(proposal_id="2736", calib_level=3, dataproduct_type="image")
    assert len(obs) == 12  # 6 NIRCam + 4 MIRI + 2 NIRISS (checked 2026-10-07)
    assert set(obs["dataRights"]) == {"PUBLIC"}
    assert NIRCAM_F200W in set(obs["obs_id"])


@pytest.mark.network
def test_live_list_products_excludes_level2_members():
    obs = query.query_observations(obs_id=NIRCAM_F200W)
    prods = query.list_products(obs, cloud_uris=True)
    # MAST lists 1,000+ rows for this observation, including 72 level-2 per-detector i2d files
    assert sorted(prods["productSubGroupDescription"]) == ["CAT", "I2D"]
    assert all(np.char.startswith(np.asarray(prods["productFilename"], dtype=str), NIRCAM_F200W))
    assert set(prods["calib_level"]) == {3}
    key = "s3://stpubdata/jwst/public/jw02736/L3/t/o001/"
    assert all(u.startswith(key) for u in prods["cloud_uri"])

    s3fs = pytest.importorskip("s3fs")
    fs = s3fs.S3FileSystem(anon=True)
    cat = prods[prods["productSubGroupDescription"] == "CAT"][0]
    assert fs.info(cat["cloud_uri"].removeprefix("s3://"))["size"] == cat["size"]
