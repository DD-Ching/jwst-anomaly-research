"""Tests for classify.classify_sources (D-012)."""

from __future__ import annotations

import numpy as np
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


def _locus_table(n_ref=12):
    rng = np.random.default_rng(1)
    rows = []
    for i in range(n_ref):  # catalogued stars: r50 ~2.9 px, ordinary colours
        rows.append((f"ref{i}", 2.9 + rng.normal(0, 0.05), 21.0, -0.7, -0.4, True))
    rows += [
        ("faint_star", 2.8, 23.0, -0.65, -0.35, False),  # point-like, stellar colours -> star
        ("compact_star", 2.0, 23.0, -0.65, -0.35, False),  # below r50_psf: still a star (D-016)
        ("noise_like", 1.0, 23.0, -0.65, -0.35, False),  # under the r50 floor -> stays
        ("brown_dwarf", 2.9, 22.0, 1.5, 0.8, False),  # point-like, unusual colours -> stays
        ("galaxy", 6.0, 22.0, -0.7, -0.4, False),  # extended -> stays
        ("too_faint", 2.9, 26.0, -0.7, -0.4, False),  # below mag_max -> stays
    ]
    t = Table(
        rows=[(u, r, m, 20.0 + c1, 20.0, 20.0 + c2, 20.0) for u, r, m, c1, c2, _ in rows],
        names=(
            "source_uid",
            "dja05_r50_pix",
            "dja05_mag_auto",
            "f150w_dja05_abmag",
            "f444w_dja05_abmag",
            "f200w_dja05_abmag",
            "f356w_dja05_abmag",
        ),
    )
    known = np.array([k for *_, k in rows])
    return t, known


def test_stellar_locus_membership():
    t, known = _locus_table()
    member, info = classify.stellar_locus(t, "dja05", known)
    got = dict(zip(t["source_uid"], member, strict=True))
    assert got["faint_star"] and got["compact_star"] and not got["brown_dwarf"]
    assert not got["galaxy"] and not got["too_faint"] and not got["noise_like"]
    assert info["n_calibration_stars"] == 12 and info["r50_psf_pix"] == pytest.approx(2.9, abs=0.1)
    assert info["provenance"] == "assumption"


def test_stellar_locus_needs_calibration_stars():
    t, known = _locus_table(n_ref=3)
    with pytest.raises(ValueError, match="catalogued stars"):
        classify.stellar_locus(t, "dja05", known)


def test_star_basis_is_wide_enough_for_later_bases():
    out = classify.classify_sources(_sources(["a"]), _matches([]))
    out["star_basis"][0] = "stellar_locus"
    assert out["star_basis"][0] == "stellar_locus"


def test_stellar_locus_calibration_needs_finite_colours():
    t, known = _locus_table(n_ref=12)
    t["f444w_dja05_abmag"][:5] = np.nan  # 7 calibration stars left with every colour
    with pytest.raises(ValueError, match="found 7"):
        classify.stellar_locus(t, "dja05", known)


def test_stellar_locus_rejects_unknown_keys_and_reads_masked_columns():
    t, known = _locus_table()
    with pytest.raises(ValueError, match="r50_tol"):
        classify.stellar_locus(t, "dja05", known, r50_tol=0.1)
    masked = Table(t, masked=True)
    k = list(t["source_uid"]).index("faint_star")
    masked["dja05_r50_pix"].mask[k] = True  # faint_star loses its size
    member, _ = classify.stellar_locus(masked, "dja05", known)
    assert not member[k]


def _populations(t, known):
    pops = Table(
        {
            "source_uid": list(t["source_uid"]),
            "population": ["star" if k else "other" for k in known],
            "star_basis": np.array(["gaia_position" if k else "" for k in known], dtype="U32"),
        }
    )
    pops.meta.update(provenance="derived", source="catalogue classification")
    return pops


def test_apply_stellar_locus_adds_members_and_records_the_catalog():
    t, known = _locus_table()
    t.meta["matched_photometry"] = {"catalog": "DJA catalog x_phot.fits"}
    pops = _populations(t, known)
    out = classify.apply_stellar_locus(pops, t, "dja05", {"enabled": True})
    basis = dict(zip(out["source_uid"], out["star_basis"], strict=True))
    assert basis["faint_star"] == "stellar_locus" and basis["ref0"] == "gaia_position"
    k = list(t["source_uid"]).index("faint_star")
    assert basis["brown_dwarf"] == "" and pops["star_basis"][k] == ""  # input untouched
    assert out.meta["stellar_locus"]["n_added"] == 2  # faint_star and compact_star
    assert "DJA catalog x_phot.fits" in out.meta["source"]
    assert (
        classify.apply_stellar_locus(pops, t, "dja05", True).meta["stellar_locus"]["n_added"] == 2
    )
    with pytest.raises(ValueError, match="mapping"):
        classify.apply_stellar_locus(pops, t, "dja05", ["colours"])
    narrow = pops.copy()
    narrow["star_basis"] = np.asarray(narrow["star_basis"]).astype("U13")
    out = classify.apply_stellar_locus(narrow, t, "dja05", True)
    assert "stellar_locus" in list(out["star_basis"])  # not truncated


def test_stellar_locus_two_sided_with_a_high_floor():
    t, known = _locus_table()
    member, info = classify.stellar_locus(t, "dja05", known, r50_floor=0.8)  # the D-015 rule
    got = dict(zip(t["source_uid"], member, strict=True))
    assert got["faint_star"] and not got["compact_star"]
    assert info["thresholds"]["r50_floor"] == 0.8
    with pytest.raises(ValueError, match="r50_floor"):
        classify.stellar_locus(t, "dja05", known, r50_floor=2.0)
