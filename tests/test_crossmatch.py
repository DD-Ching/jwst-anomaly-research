"""Tests for jwst_anomaly.crossmatch.

Offline tests replay raw service responses recorded from the live services for five real
positions in the SMACS J0723.3-7327 field (tests/data/xmatch_*.ecsv; provenance in each file's
header). Re-record with ``python tests/test_crossmatch.py --record`` (needs network).
"""

from __future__ import annotations

import sys
from datetime import UTC, datetime
from pathlib import Path

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.table import Table

from jwst_anomaly import crossmatch as xm
from jwst_anomaly import schema

DATA = Path(__file__).parent / "data"
TARGETS_FILE = DATA / "xmatch_targets.ecsv"
RADIUS = 1.0
GAIA_STAR = "gaia_dr3_5262900027979791104"  # G=14.6, saturated in F200W: Gaia position
GAIA_STAR_CAT = "f200w_2596"  # F200W point source matched to Gaia DR3 5262899954964194816
CLUSTER = "simbad_smacs0723_cluster"
BCG = "f200w_1205_bcg"
LENSED = "f200w_1793"  # F200W counterpart of SIMBAD LeG [CLL2024] SMACS J0723.3-7327 15.2


def raw_file(service: str, backend: str) -> Path:
    return DATA / f"xmatch_raw_{service}_{backend}.ecsv"


@pytest.fixture
def targets() -> Table:
    return Table.read(TARGETS_FILE, format="ascii.ecsv")


def recorded(service: str, backend: str):
    def fetch(coords, radius_arcsec, opts):
        assert radius_arcsec == RADIUS
        return Table.read(raw_file(service, backend), format="ascii.ecsv")

    return fetch


def failing(service: str, backend: str, calls: list | None = None):
    def fetch(coords, radius_arcsec, opts):
        if calls is not None:
            calls.append((service, backend))
        raise ConnectionError(f"{service}/{backend} unreachable (simulated)")

    return fetch


@pytest.fixture
def replay(monkeypatch):
    """Point every backend at its recorded response; return a setter for overrides."""
    monkeypatch.setattr(xm.time, "sleep", lambda s: None)
    for key in xm._FETCHERS:
        monkeypatch.setitem(xm._FETCHERS, key, recorded(*key))

    def override(key, fetch):
        monkeypatch.setitem(xm._FETCHERS, key, fetch)

    return override


def by_uid(table: Table) -> dict:
    return {row["source_uid"]: row for row in table}


def check_smacs0723_summary(summary: Table) -> None:
    """Expectations that hold for the recorded and the live responses."""
    schema.validate(summary, schema.XMATCH_COLUMNS)
    assert summary.meta["provenance"] == schema.Provenance.OBSERVED.value
    rows = by_uid(summary)
    for uid in (GAIA_STAR, GAIA_STAR_CAT):
        assert rows[uid]["is_star"], uid
        assert rows[uid]["n_gaia"] >= 1
    assert rows[CLUSTER]["best_match_id"] == "SMACS J0723.3-7327"
    assert rows[CLUSTER]["n_simbad"] >= 1 and rows[CLUSTER]["n_ned"] >= 1
    assert rows[BCG]["is_known_object"] and not rows[BCG]["is_star"]
    assert rows[BCG]["n_simbad"] + rows[BCG]["n_ned"] >= 1
    assert rows[LENSED]["is_lens_related"]
    assert "simbad:LeG" in rows[LENSED]["lens_types"]
    for row in summary:
        if row["is_known_object"]:
            assert 0 <= row["best_match_sep_arcsec"] <= RADIUS
            assert row["best_match_service"] in xm.SERVICES


# --- Offline: recorded responses ----------------------------------------------------------


def test_recorded_primary_backends(targets, replay):
    summary = xm.crossmatch(targets, RADIUS, cache=False)
    check_smacs0723_summary(summary)
    assert list(summary["source_uid"]) == list(targets["source_uid"])
    for col in xm.SUMMARY_EXTRA_COLUMNS:
        assert col in summary.colnames
    meta = summary.meta
    assert meta["services_ok"] == ["simbad", "ned", "gaia"] and meta["services_failed"] == []
    assert meta["backend"] == {"simbad": "xmatch", "ned": "tap", "gaia": "xmatch"}
    assert "vizier:I/355/gaiadr3" in meta["catalog"]["gaia"]
    assert set(meta["query_utc"]) == {"simbad", "ned", "gaia"}
    assert meta["from_cache"] == {"simbad": False, "ned": False, "gaia": False}
    assert "ASSUMPTION" in meta["star_criteria"]
    assert "5 targets" in meta["source"]
    rows = by_uid(summary)
    assert rows[GAIA_STAR]["best_match_id"] == "Gaia DR3 5262900027979791104"
    assert rows[GAIA_STAR]["best_match_type"] == "astrometric_star"


def test_recorded_fallback_backends_agree(targets, replay):
    calls: list = []
    for svc, backend in (("simbad", "xmatch"), ("ned", "tap"), ("gaia", "xmatch")):
        replay((svc, backend), failing(svc, backend, calls))
    summary = xm.crossmatch(targets, RADIUS, cache=False, retries=0)
    check_smacs0723_summary(summary)
    assert summary.meta["backend"] == {"simbad": "tap", "ned": "cone", "gaia": "tap"}
    assert set(summary.meta["service_errors"]) == {"simbad", "ned", "gaia"}
    assert "ConnectionError" in summary.meta["service_errors"]["gaia"]
    assert len(calls) == 3


def test_query_matches_long_format(targets, replay):
    matches = xm.query_matches(targets, RADIUS, cache=False)
    schema.validate(matches, xm.MATCH_COLUMNS)
    assert set(matches["service"]) <= set(xm.SERVICES)
    assert np.all(matches["sep_arcsec"] <= RADIUS)
    gaia = matches[matches["service"] == "gaia"]
    assert np.all(np.isfinite(gaia["parallax_over_error"]))
    star = gaia[gaia["source_uid"] == GAIA_STAR_CAT][0]
    assert star["match_id"] == "Gaia DR3 5262899954964194816"
    assert star["parallax_over_error"] > xm.GAIA_MIN_PARALLAX_SNR
    lensed = matches[(matches["source_uid"] == LENSED) & (matches["service"] == "simbad")][0]
    assert lensed["match_type"] == "LeG" and "LeI" in lensed["all_types"]
    summary = xm.summarize_matches(targets, matches)
    assert list(summary["n_matches"]) == [
        int(np.sum(matches["source_uid"] == uid)) for uid in targets["source_uid"]
    ]


def test_failed_service_is_recorded_and_others_continue(targets, replay):
    replay(("ned", "tap"), failing("ned", "tap"))
    replay(("ned", "cone"), failing("ned", "cone"))
    summary = xm.crossmatch(targets, RADIUS, cache=False, retries=0)
    assert summary.meta["services_ok"] == ["simbad", "gaia"]
    assert summary.meta["services_failed"] == ["ned"]
    assert "tap: ConnectionError" in summary.meta["service_errors"]["ned"]
    assert "cone: ConnectionError" in summary.meta["service_errors"]["ned"]
    assert np.all(summary["n_ned"] == -1)
    assert by_uid(summary)[GAIA_STAR]["is_star"]
    assert "ned" not in summary.meta["source"].split("[")[0]


def test_all_services_failing_raises(targets, replay):
    for key in xm._FETCHERS:
        replay(key, failing(*key))
    with pytest.raises(xm.CrossmatchError, match="all external services failed"):
        xm.crossmatch(targets, RADIUS, cache=False, retries=0)


def test_retries_with_exponential_backoff(targets, replay, monkeypatch):
    delays: list[float] = []
    monkeypatch.setattr(xm.time, "sleep", delays.append)
    attempts = {"n": 0}

    def flaky(coords, radius_arcsec, opts):
        attempts["n"] += 1
        if attempts["n"] < 3:
            raise TimeoutError("simulated timeout")
        return recorded("gaia", "xmatch")(coords, radius_arcsec, opts)

    replay(("gaia", "xmatch"), flaky)
    summary = xm.crossmatch(targets, RADIUS, ["gaia"], cache=False, retries=2, backoff_s=0.5)
    assert attempts["n"] == 3 and delays == [0.5, 1.0]
    assert summary.meta["backend"] == {"gaia": "xmatch"}
    assert summary.meta["service_errors"] == {}
    assert summary.meta["services_ok"] == ["gaia"]
    assert np.all(summary["n_simbad"] == -1) and np.all(summary["n_ned"] == -1)


def test_cache_roundtrip_refresh_and_expiry(targets, replay, tmp_path):
    first = xm.crossmatch(targets, RADIUS, cache_dir=tmp_path)
    assert len(list(tmp_path.glob("*.ecsv"))) == 3
    for key in xm._FETCHERS:
        replay(key, failing(*key))
    cached = xm.crossmatch(targets, RADIUS, cache_dir=tmp_path)
    assert cached.meta["from_cache"] == {"simbad": True, "ned": True, "gaia": True}
    assert cached.meta["query_utc"] == first.meta["query_utc"]
    for col in schema.XMATCH_COLUMNS:
        assert list(cached[col]) == list(first[col]), col
    with pytest.raises(xm.CrossmatchError):
        xm.crossmatch(targets, RADIUS, cache_dir=tmp_path, refresh=True, retries=0)
    with pytest.raises(xm.CrossmatchError):
        xm.crossmatch(targets, RADIUS, cache_dir=tmp_path, cache_max_age_days=0, retries=0)
    moved = targets.copy()
    moved["dec"][0] += 1.0 / 3600
    with pytest.raises(xm.CrossmatchError):  # different positions -> different cache key
        xm.crossmatch(moved, RADIUS, cache_dir=tmp_path, retries=0)


def test_empty_targets_make_no_queries(replay):
    for key in xm._FETCHERS:
        replay(key, failing(*key))
    empty = Table({"source_uid": np.array([], dtype=str), "ra": [], "dec": []})
    empty.meta["source"] = "nothing"
    summary = xm.crossmatch(empty, cache=False)
    assert len(summary) == 0
    schema.validate(summary, schema.XMATCH_COLUMNS)


def test_quantity_columns_and_extra_input_columns(targets, replay):
    qt = targets.copy()
    qt["ra"] = (np.asarray(targets["ra"]) * u.deg).to(u.rad)
    qt["dec"] = (np.asarray(targets["dec"]) * u.deg).to(u.rad)
    a = xm.crossmatch(targets, RADIUS, cache=False)
    b = xm.crossmatch(qt, RADIUS, cache=False)
    assert list(a["best_match_id"]) == list(b["best_match_id"])
    assert np.allclose(a["best_match_sep_arcsec"], b["best_match_sep_arcsec"], equal_nan=True)


@pytest.mark.parametrize(
    "mutate, kwargs, match",
    [
        (lambda t: t.remove_column("dec"), {}, "missing required columns"),
        (lambda t: t["ra"].__setitem__(0, np.nan), {}, "non-finite"),
        (lambda t: t["source_uid"].__setitem__(1, t["source_uid"][0]), {}, "unique"),
        (lambda t: None, {"services": ["simbad", "vizier"]}, "services must be"),
        (lambda t: None, {"services": []}, "services must be"),
        (lambda t: None, {"backends": {"gaia": "cone"}}, "unknown backend"),
        (lambda t: None, {"radius_arcsec": 0.0}, "positive"),
        (lambda t: None, {"radius_arcsec": float("nan")}, "positive"),
        (lambda t: None, {"retries": -1}, "retries >= 0"),
        (lambda t: None, {"timeout_s": 0}, "timeout_s > 0"),
    ],
)
def test_input_validation(targets, replay, mutate, kwargs, match):
    t = targets.copy()
    mutate(t)
    with pytest.raises(ValueError, match=match):
        xm.crossmatch(t, **{"radius_arcsec": RADIUS, "cache": False, **kwargs})


# --- Offline: classification rules ---------------------------------------------------------


def test_gaia_astrometric_snr_uses_correlation():
    plx, pm = xm.gaia_astrometric_snr(
        np.array([2.0, np.nan]),
        np.array([0.5, np.nan]),
        np.array([3.0, np.nan]),
        np.array([1.0, np.nan]),
        np.array([4.0, np.nan]),
        np.array([1.0, np.nan]),
        np.array([np.nan, np.nan]),
    )
    assert plx[0] == pytest.approx(4.0) and pm[0] == pytest.approx(5.0)
    assert np.isnan(plx[1]) and np.isnan(pm[1])
    _, pm_corr = xm.gaia_astrometric_snr(*[np.array([v]) for v in (0, 1, 3, 1, 4, 1, 0.5)])
    expected = np.sqrt((9 - 2 * 0.5 * 12 + 16) / (1 - 0.25))
    assert pm_corr[0] == pytest.approx(expected)


def objects(**cols) -> Table:
    n = len(next(iter(cols.values())))
    t = Table()
    for col in xm._STR_COLS:
        t[col] = np.array(cols.get(col, [""] * n), dtype=str)
    for col in xm._FLOAT_COLS:
        t[col] = np.array(cols.get(col, [np.nan] * n), dtype=float)
    return t


def test_gaia_star_thresholds():
    objs = objects(
        match_id=["plx", "pm", "neither", "two_param"],
        parallax=[1.0, 0.1, 0.1, np.nan],
        parallax_error=[0.2, 0.1, 0.1, np.nan],  # S/N 5, 1, 1
        pmra=[0.0, 5.0, 1.0, np.nan],
        pmra_error=[1.0, 1.0, 1.0, np.nan],  # pm S/N 0, 5, 1
        pmdec=[0.0, 0.0, 0.0, np.nan],
        pmdec_error=[1.0, 1.0, 1.0, np.nan],
    )
    out = xm._classify("gaia", objs)
    assert list(out["is_star"]) == [True, True, False, False]
    assert list(out["match_type"]) == [
        "astrometric_star",
        "astrometric_star",
        "unclassified",
        "unclassified",
    ]
    assert not out["is_lens_related"].any()


def test_simbad_and_ned_classification():
    simbad = xm._classify(
        "simbad",
        objects(
            match_type=["PM*", "G", "SN*", "LS?", "**?", "BD*"],
            all_types=["PM*|*", "G|G?|LeG", "SN*", "LS?|G", "**?", "BD*|LM*"],
        ),
    )
    assert list(simbad["is_star"]) == [True, False, False, False, False, True]
    assert list(simbad["lens_types"]) == ["", "LeG", "", "LS?", "", ""]
    ned = xm._classify("ned", objects(match_type=["*", "G_Lens", "Q_Lens", "G"]))
    assert not ned["is_star"].any()
    assert list(ned["is_lens_related"]) == [False, True, True, False]


def test_derive_simbad_otype_sets_rule():
    otypedef = Table(
        rows=[
            ("*", "*", 0),
            ("PM*", "* > PM*", 0),
            ("BD?", "* > LM* > BD*", 1),
            ("PN", "* > Ev* > PN", 0),
            ("HH", "* > Y*O > out > HH", 0),
            ("G", "G", 0),
            ("LeG", "grv > gLS > LeI > LeG", 0),
            ("Le?", "grv > gLS > gLe", 1),
            ("Lev", "grv > Lev", 0),
            ("BH", "grv > BH", 0),
        ],
        names=("otype", "path", "is_candidate"),
    )
    stellar, lens = xm.derive_simbad_otype_sets(otypedef)
    assert stellar == {"*", "PM*"}
    assert lens == {"LeG", "Le?", "Lev"}
    assert stellar <= xm.SIMBAD_STAR_OTYPES and lens <= xm.SIMBAD_LENS_OTYPES


def test_adql_cones_and_chunking():
    coords = SkyCoord(np.linspace(10, 11, 250), np.full(250, -73.0), unit="deg")
    chunks = xm._chunks(coords, xm._ADQL_CHUNK)
    assert [len(c) for c in chunks] == [100, 100, 50]
    adql = xm._adql_cones(coords[:2], 1.8)
    assert adql.count("CONTAINS(") == 2 and adql.count(" OR ") == 1
    assert "CIRCLE('ICRS', 10.000000000, -73.000000000, 0.000500000)" in adql


def match_rows(uid, rows, services_ok=("simbad", "ned", "gaia")) -> Table:
    """Long-format matches for one target from (service, id, type, sep, is_star) tuples."""
    service, match_id, match_type, sep, star = zip(*rows, strict=True)
    matches = Table(
        {
            "source_uid": [uid] * len(rows),
            "service": list(service),
            "match_id": list(match_id),
            "match_type": list(match_type),
            "all_types": [""] * len(rows),
            "lens_types": [""] * len(rows),
            "sep_arcsec": list(sep),
            "is_star": list(star),
        }
    )
    for col in xm.MATCH_COLUMNS:
        if col not in matches.colnames:
            matches[col] = False if col.startswith("is_") else np.nan
    matches.meta.update(
        provenance="observed",
        source="unit test",
        services_requested=["simbad", "ned", "gaia"],
        services_ok=list(services_ok),
    )
    return matches


def test_best_match_prefers_smaller_separation_then_service_order(targets):
    t = targets[:1]
    rows = [("gaia", "g", "unclassified", 0.3, False), ("ned", "n", "G", 0.3, False)]
    rows.append(("simbad", "s", "G", 0.5, False))
    row = xm.summarize_matches(t, match_rows(t["source_uid"][0], rows))[0]
    assert row["best_match_id"] == "n" and row["best_match_service"] == "ned"
    assert row["n_matches"] == 3 and row["n_gaia"] == 1


@pytest.mark.parametrize(
    "rows, expected",
    [
        # A galaxy with a Gaia star 0.9" away is not a star: the nearest SIMBAD/Gaia match is G.
        ([("simbad", "s", "G", 0.1, False), ("gaia", "g", "astrometric_star", 0.9, True)], False),
        # NED cannot classify stars, so a nearer NED IR source does not hide the Gaia star.
        ([("ned", "n", "IrS", 0.05, False), ("gaia", "g", "astrometric_star", 0.08, True)], True),
        ([("ned", "n", "*", 0.05, False)], False),
    ],
)
def test_is_star_follows_nearest_simbad_or_gaia_match(targets, rows, expected):
    t = targets[:1]
    row = xm.summarize_matches(t, match_rows(t["source_uid"][0], rows))[0]
    assert row["is_star"] is np.bool_(expected)


def test_summarize_requires_query_meta(targets):
    matches = match_rows(targets["source_uid"][0], [("ned", "n", "G", 0.1, False)])
    del matches.meta["services_ok"]
    with pytest.raises(ValueError, match="services_ok"):
        xm.summarize_matches(targets, matches)


def test_fallback_cache_does_not_mask_healthy_primary(targets, replay, tmp_path):
    replay(("ned", "tap"), failing("ned", "tap"))
    first = xm.crossmatch(targets, RADIUS, ["ned"], cache_dir=tmp_path, retries=0)
    assert first.meta["backend"] == {"ned": "cone"}
    replay(("ned", "tap"), recorded("ned", "tap"))
    second = xm.crossmatch(targets, RADIUS, ["ned"], cache_dir=tmp_path)
    assert second.meta["backend"] == {"ned": "tap"}
    assert second.meta["from_cache"] == {"ned": False}
    replay(("ned", "tap"), failing("ned", "tap"))  # primary down again: its cache entry is used
    third = xm.crossmatch(targets, RADIUS, ["ned"], cache_dir=tmp_path, retries=0)
    assert third.meta["backend"] == {"ned": "tap"} and third.meta["from_cache"] == {"ned": True}


def test_corrupt_cache_entry_is_a_miss(targets, replay, tmp_path):
    xm.crossmatch(targets, RADIUS, ["gaia"], cache_dir=tmp_path)
    (entry,) = tmp_path.glob("gaia_xmatch_*.ecsv")
    entry.write_text("not an ecsv file")
    again = xm.crossmatch(targets, RADIUS, ["gaia"], cache_dir=tmp_path)
    assert again.meta["from_cache"] == {"gaia": False}


def ned_tap_empty() -> Table:
    return Table(
        names=list(xm.RAW_COLUMNS[("ned", "tap")].values()), dtype=(str,) * 2 + (float,) * 3
    )


def test_retries_are_per_request_not_per_service(replay, monkeypatch):
    delays: list[float] = []
    monkeypatch.setattr(xm.time, "sleep", delays.append)
    n = xm._ADQL_CHUNK + 20
    many = Table({"source_uid": [f"t{i}" for i in range(n)], "ra": np.full(n, 10.0)})
    many["dec"] = np.linspace(-30, -29, n)
    calls: list[int] = []

    def flaky_second_chunk(coords, radius_arcsec, opts):
        calls.append(len(coords))
        if len(calls) == 2:
            raise ConnectionError("simulated drop on chunk 2")
        return ned_tap_empty()

    replay(("ned", "tap"), flaky_second_chunk)
    out = xm.crossmatch(many, RADIUS, ["ned"], cache=False, retries=1, backoff_s=0.25)
    assert calls == [xm._ADQL_CHUNK, 20, 20] and delays == [0.25]
    assert out.meta["backend"] == {"ned": "tap"} and not out["is_known_object"].any()


def test_permanent_errors_are_not_retried(replay, monkeypatch):
    delays: list[float] = []
    monkeypatch.setattr(xm.time, "sleep", delays.append)
    n = xm._MAX_TARGETS[("ned", "cone")] + 1
    many = Table({"source_uid": [f"t{i}" for i in range(n)], "ra": np.full(n, 10.0)})
    many["dec"] = np.linspace(-30, -29, n)
    tap_calls: list = []

    def truncated(coords, radius_arcsec, opts):
        tap_calls.append(1)
        raise xm._PermanentError("NED TAP result truncated")

    cone_calls: list = []
    replay(("ned", "tap"), truncated)
    replay(("ned", "cone"), failing("ned", "cone", cone_calls))
    with pytest.raises(xm.CrossmatchError, match="capped at 100 targets"):
        xm.crossmatch(many, RADIUS, ["ned"], cache=False, retries=3)
    assert len(tap_calls) == 1 and cone_calls == [] and delays == []


def test_ned_cone_empty_result_has_correct_dtypes(monkeypatch):
    from astroquery.ipac.ned import Ned

    monkeypatch.setattr(Ned, "query_region", lambda *a, **k: Table())
    coords = SkyCoord([10.0], [-30.0], unit="deg")
    opts = xm._Options(timeout_s=5.0, retries=0, backoff_s=0.0)
    empty = xm._fetch_ned_cone(coords, RADIUS, opts)
    assert len(empty) == 0
    assert empty["Type"].dtype.kind == "U" and empty["DEC"].dtype.kind == "f"
    assert len(xm._normalize("ned", "cone", empty)) == 0


def test_backend_names_are_case_insensitive(targets, replay):
    out = xm.crossmatch(targets, RADIUS, ["Gaia"], cache=False, backends={"GAIA": "TAP"})
    assert out.meta["backend"] == {"gaia": "tap"}


# --- Live services ----------------------------------------------------------------------------


@pytest.mark.network
def test_live_smacs0723(targets):
    summary = xm.crossmatch(targets, RADIUS, cache=False)
    assert summary.meta["services_failed"] == [], summary.meta["service_errors"]
    check_smacs0723_summary(summary)


@pytest.mark.network
def test_live_fallback_backends(targets):
    backends = {"simbad": "tap", "ned": "cone", "gaia": "tap"}
    summary = xm.crossmatch(targets, RADIUS, cache=False, backends=backends)
    assert summary.meta["backend"] == backends, summary.meta["service_errors"]
    check_smacs0723_summary(summary)


@pytest.mark.network
def test_live_simbad_otype_sets_are_current():
    from astroquery.simbad import Simbad

    otypedef = Simbad.query_tap("SELECT otype, path, is_candidate FROM otypedef")
    stellar, lens = xm.derive_simbad_otype_sets(otypedef)
    assert stellar == xm.SIMBAD_STAR_OTYPES
    assert lens == xm.SIMBAD_LENS_OTYPES


# --- Fixture recording (manual) -------------------------------------------------------------


def record_fixtures() -> None:
    """Re-record tests/data/xmatch_raw_*.ecsv from the live services (network)."""
    from importlib.metadata import version

    targets = Table.read(TARGETS_FILE, format="ascii.ecsv")
    coords = SkyCoord(targets["ra"], targets["dec"], unit="deg")
    opts = xm._Options(timeout_s=60.0, retries=2, backoff_s=2.0)
    for (service, backend), fetch in xm._FETCHERS.items():
        raw = fetch(coords, RADIUS, opts)
        keep = [c for c in ("xm_idx", "angDist") if c in raw.colnames]
        raw = raw[keep + list(xm.RAW_COLUMNS[(service, backend)].values())]
        raw.meta.clear()
        raw.meta.update(
            provenance=schema.Provenance.OBSERVED.value,
            source=xm.CATALOGS[(service, backend)],
            note=(
                f"raw {service}/{backend} response for {TARGETS_FILE.name} at radius "
                f"{RADIUS} arcsec, trimmed to the columns crossmatch.py reads"
            ),
            recorded_utc=datetime.now(UTC).isoformat(timespec="seconds"),
            software={p: version(p) for p in ("astroquery", "pyvo", "astropy")},
        )
        raw.write(raw_file(service, backend), format="ascii.ecsv", overwrite=True)
        print(f"{service}/{backend}: {len(raw)} rows -> {raw_file(service, backend).name}")


if __name__ == "__main__":
    if "--record" in sys.argv:
        record_fixtures()
