"""Tests for classify.classify_sources (D-012)."""

from __future__ import annotations

import pytest
from astropy.table import Table

from jwst_anomaly import classify, schema


def _sources(uids):
    t = Table(
        {
            "source_uid": uids,
            "ra": [10.0 + i * 1e-3 for i in range(len(uids))],
            "dec": [-5.0] * len(uids),
            "ref_band": ["F200W"] * len(uids),
            "n_bands": [3] * len(uids),
        }
    )
    t.meta.update(provenance="derived", source="test")
    return t


def _matches(rows, failed=(), radius=0.5):
    t = Table(
        rows=rows,
        names=("source_uid", "service", "match_id", "match_type", "sep_arcsec", "is_star"),
        dtype=(str, str, str, str, float, bool),
    )
    t.meta.update(
        provenance="observed",
        source="test",
        services_requested=["gaia", "simbad"],
        services_failed=list(failed),
        radius_arcsec=radius,
    )
    return t


def test_star_rules():
    src = _sources(["astro", "simbad_star", "gaia_close", "gaia_far", "galaxy", "none"])
    m = _matches(
        [
            ("astro", "gaia", "G1", "astrometric_star", 0.05, True),
            ("simbad_star", "simbad", "S1", "*", 0.2, True),
            ("gaia_close", "gaia", "G2", "unclassified", 0.2, False),
            ("gaia_far", "gaia", "G3", "unclassified", 0.45, False),
            ("galaxy", "simbad", "S2", "G", 0.1, False),
        ]
    )
    out = classify.classify_sources(src, m)
    pop = dict(zip(out["source_uid"], out["population"], strict=True))
    basis = dict(zip(out["source_uid"], out["star_basis"], strict=True))
    assert pop == {
        "astro": "star",
        "simbad_star": "star",
        "gaia_close": "star",
        "gaia_far": "other",
        "galaxy": "other",
        "none": "other",
    }
    assert basis["astro"] == "simbad_or_gaia_astrometry"
    assert basis["gaia_close"] == "gaia_position"
    assert out.meta["provenance"] == "derived"
    assert out.meta["thresholds"]["provenance"] == "assumption"
    schema.validate(out, schema.CLASSIFY_COLUMNS)


def test_nearest_match_decides():
    src = _sources(["galaxy_near_star", "star_behind_galaxy_entry"])
    m = _matches(
        [
            ("galaxy_near_star", "simbad", "S", "G", 0.05, False),
            ("galaxy_near_star", "gaia", "G1", "astrometric_star", 0.45, True),
            ("star_behind_galaxy_entry", "gaia", "G2", "astrometric_star", 0.05, True),
            ("star_behind_galaxy_entry", "simbad", "S2", "G", 0.4, False),
        ]
    )
    out = classify.classify_sources(src, m)
    assert list(out["population"]) == ["other", "star"]


def test_ned_failure_does_not_block_classification():
    src = _sources(["a"])
    m = _matches([("a", "gaia", "G", "astrometric_star", 0.1, True)], failed=["ned"])
    assert classify.classify_sources(src, m)["population"][0] == "star"


def test_gaia_radius_cannot_exceed_query_radius():
    with pytest.raises(ValueError, match="exceeds the query radius"):
        classify.classify_sources(_sources(["a"]), _matches([]), gaia_radius_arcsec=1.0)


def test_no_matches_means_no_stars():
    out = classify.classify_sources(_sources(["a", "b"]), _matches([]))
    assert list(out["population"]) == ["other", "other"]


def test_failed_service_raises():
    with pytest.raises(ValueError, match="services failed"):
        classify.classify_sources(_sources(["a"]), _matches([], failed=["gaia"]))
