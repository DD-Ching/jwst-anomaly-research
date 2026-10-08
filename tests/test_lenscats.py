"""Offline tests for the published-lens-catalogue adapter and the no-visible-deflector screen."""

import math

import numpy as np
import pytest
from astropy.cosmology import Planck18
from astropy.table import Table

from jwst_anomaly import lenscats, schema, signatures

LENSCAT = """name,RA [deg],DEC [deg],zlens,type,grading,ref
J0011-0845,2.83435,-8.7643,-,galaxy,confident,lemon
SDSSJ0100+0000,15.0,0.0,0.31 | 0.3102,galaxy,confident,slacs
DESI-dark,30.0,10.0,not measured,galaxy,probable,storfer
CL0001,40.0,-5.0,0.4,cluster,confident,clash
AGEL 002527+101107,6.364166667,0.185305556,L,galaxy,confident,agel
NSCS J000001+051909,0.00504,5.31906,not measured,galaxy,probable,Lopes (2004) doi:10.1086/423038
"""
SUGOHI = """\
010000+000000,15.0002,0.0001,0.310,-99.000,-99.00,-99.00,1.20,19.50,-99.00,GG,V,A,SuGOHI1
020000+100000,30.0003,10.0002,-99.000,2.100,-99.00,-99.00,0.80,-99.00,-99.00,GQ,CHITAH,B,SuGOHI9
"""
EUCLID = (
    "subset,id_str,object_id,tile_index,right_ascension,declination,segmentation_area,"
    "flux_vis_1fwhm_aper,expert_score,grade,expert_total_votes,notes,references\n"
    """discovery_engine,1_A,1,1,60.0,-50.0,100.0,2.0,2.5,A,10,,
discovery_engine,2_B,2,1,61.0,-50.0,100.0,2.0,0.5,C,10,possible cluster,
"""
)
EUCLID_MASS = """,id_str,einstein_radius_effective_median_pdf
0,1_A,0.9
"""
EUCLID_MAG = """,id_str,vis_lens_magnitude_ab_median_pdf
0,1_A,20.5
"""


@pytest.fixture
def tables(tmp_path):
    files = {}
    for name, text in [
        ("lenscat.csv", LENSCAT),
        ("sugohi.csv", SUGOHI),
        ("euclid.csv", EUCLID),
        ("mass.csv", EUCLID_MASS),
        ("mag.csv", EUCLID_MAG),
    ]:
        (tmp_path / name).write_text(text)
        files[name] = tmp_path / name
    return [
        lenscats.read_lenscat(files["lenscat.csv"]),
        lenscats.read_euclid_q1(files["euclid.csv"], files["mass.csv"], files["mag.csv"]),
        lenscats.read_sugohi(files["sugohi.csv"]),
    ]


def test_readers_normalise(tables):
    lc, eu, su = tables
    assert list(lc.colnames) == list(lenscats.ENTRY_COLUMNS)
    assert lc.meta["provenance"] == schema.Provenance.OBSERVED.value
    assert list(lc["grade"]) == ["A", "A", "B", "A", "A", "B"]
    assert list(lc["scale"]) == ["galaxy", "galaxy", "galaxy", "cluster", "galaxy", "cluster"]
    assert lc["z_lens"][1] == pytest.approx(0.31)
    assert list(lc["lens_z_known"]) == [
        False,
        True,
        False,
        True,
        True,
        False,
    ]  # "L" counts as known
    assert su["scale"][1] == "galaxy" and su["system_type"][1] == "GQ"
    assert np.isnan(su["z_lens"][1]) and su["z_source"][1] == pytest.approx(2.1)
    assert su["theta_e"][0] == pytest.approx(1.2) and su["lens_mag"][0] == pytest.approx(19.5)
    assert eu["theta_e"][0] == pytest.approx(0.9) and np.isnan(eu["theta_e"][1])
    assert eu["lens_mag"][0] == pytest.approx(20.5) and eu["scale"][1] == "group"


def test_merge_dedups_and_keeps_provenance(tables):
    s = lenscats.merge(tables, radius_arcsec=3.0)
    assert len(s) == 6 + 2 + 2 - 2  # two SuGOHI entries within 3'' of lenscat entries
    assert s.meta["provenance"] == schema.Provenance.DERIVED.value
    assert s.meta["n_entries"] == 10
    row = s[np.array(["sugohi" in c for c in s["catalogues"]]) & (np.abs(s["ra"] - 15) < 0.01)][0]
    assert row["n_entries"] == 2
    assert set(row["entries"].split(";")) == {"lenscat:SDSSJ0100+0000", "sugohi:010000+000000"}
    assert row["theta_e"] == pytest.approx(1.2) and row["lens_mag_known"]
    dark = s[np.abs(s["ra"] - 30) < 0.01][0]
    assert not dark["lens_z_known"] and dark["z_source"] == pytest.approx(2.1)
    assert dark["grade"] == "B"  # both entries graded B


def test_merge_radius_separates(tables):
    s = lenscats.merge(tables, radius_arcsec=0.1)
    assert len(s) == 10


def test_select_no_lens_info(tables):
    s = lenscats.merge(tables)
    sel = lenscats.select_no_lens_info(s)
    picked = set(s["name"][sel])
    # Lemon quasar (zlens "-"), the dark-ish SuGOHI GQ / DESI pair; not the cluster, not "L",
    # not the Euclid group, not lenses with a lens magnitude
    assert "J0011-0845" in picked
    assert any(n in picked for n in ("DESI-dark", "020000+100000"))
    assert "CL0001" not in picked and "NSCS J000001+051909" not in picked
    assert "AGEL 002527+101107" not in picked
    assert "1_A" not in picked and "2_B" not in picked


def test_survey_adapter(tables):
    surv = lenscats.PublishedLensSurvey(tables, area=123.0)
    assert isinstance(surv, signatures.CatalogueSurvey)
    assert len(surv.catalogue()) == 8 and surv.area_deg2() == 123.0


def test_name_position_offset():
    # consistent designations: no excess beyond the truncation cell
    assert lenscats.name_position_offset("SDSSJ0029-0055", 7.4, -0.92) < 1
    assert lenscats.name_position_offset("010000+000000", 15.0002, 0.0001) == pytest.approx(
        0, abs=1
    )
    # AGEL entry with the declination degrees dropped: ~10 deg off
    off = lenscats.name_position_offset("AGEL 002527+101107", 6.364166667, 0.185305556)
    assert off > 3.5e4
    assert np.isnan(lenscats.name_position_offset("HE0435-1223", 69.56, -12.29))
    assert np.isnan(lenscats.name_position_offset("DESI-dark", 30.0, 10.0))


def test_classify_deflectors():
    systems = Table(
        {"system_id": ["a", "b", "c"], "ra": [10.0, 20.0, 30.0], "dec": [0.0, 0.0, 0.0]}
    )
    d = 1 / 3600
    sources = Table(
        {
            "ra": [10.0 + 0.5 * d, 20.0 + 0.8 * d, 20.0 - 0.8 * d, 20.0 + 4 * d],
            "dec": [0.0, 0.0, 0.0, 0.0],
            "type": ["DEV", "PSF", "PSF", "EXP"],
            "mag_z": [19.0, 21.0, 21.5, 22.0],
            "galdepth_z": [400.0, 400.0, 400.0, 400.0],
        }
    )
    out = lenscats.classify_deflectors(systems, sources, radius_arcsec=1.5, cover_radius_arcsec=5)
    a, b, c = out
    assert a["visible_deflector"] and a["sep_ext"] == pytest.approx(0.5, abs=0.01)
    assert a["mag_ext"] == 19.0 and a["covered"]
    assert not b["visible_deflector"] and b["n_psf"] == 2 and b["covered"]  # quasar images only
    assert not c["covered"] and not c["visible_deflector"]
    assert a["depth_z"] == pytest.approx(22.5 - 2.5 * math.log10(5 / 20))


def test_sis_sigma_round_trip():
    zl, zs = 0.5, 2.0
    ratio = (
        Planck18.angular_diameter_distance(zl, zs) / Planck18.angular_diameter_distance(zs)
    ).value
    theta = math.degrees(4 * math.pi * (200 / lenscats.C_KM_S) ** 2 * ratio) * 3600
    assert lenscats.sis_sigma(theta, zl, zs) == pytest.approx(200, rel=1e-6)
    assert np.isnan(lenscats.sis_sigma(1.0, 2.0, 1.5))


def test_fit_fj_recovers_and_required_mag():
    rng = np.random.default_rng(1)
    n = 400
    zl = rng.uniform(0.2, 0.8, n)
    zs = np.full(n, 2.0)
    th = rng.uniform(0.6, 2.0, n)
    truth = lenscats.FJCalibration(a=19.0, k=1.5, slope=-10.0, rms=0.0, n=0, band="z")
    mag = truth.mag(lenscats.sis_sigma(th, zl, zs), zl) + rng.normal(0, 0.3, n)
    cal = lenscats.fit_fj(th, zl, zs, mag, band="z")
    assert cal.a == pytest.approx(19.0, abs=0.06) and cal.k == pytest.approx(1.5, abs=0.2)
    assert cal.rms == pytest.approx(0.3, abs=0.05)
    m_small, z_small = lenscats.required_lens_mag([0.5], [2.0], cal)
    m_big, _ = lenscats.required_lens_mag([2.0], [2.0], cal)
    assert m_small[0] > m_big[0]  # smaller Einstein radius -> fainter lens allowed
    assert 0.1 <= z_small[0] < 1.9
    fainter, _ = lenscats.required_lens_mag([0.5], [2.0], cal, n_sigma_faint=3.0)
    assert fainter[0] == pytest.approx(m_small[0] + cal.rms, abs=1e-9)


def test_galdepth_to_mag():
    assert lenscats.galdepth_to_mag(25.0) == pytest.approx(22.5 - 2.5 * math.log10(1.0))


def test_name_offset_prefixes_and_tenths():
    # SPT designations are J2000 without a "J": a 1.7 deg declination error is caught
    assert lenscats.name_position_offset("SPT0418-47", 64.66529, -46.13536) > 3000
    assert lenscats.name_position_offset("SPT0418-47", 64.665, -47.86) == 0
    # HHMM+DDd designation (declination in tenths of a degree)
    assert lenscats.name_position_offset("0302+006", 45.6, 0.62) == 0
    assert lenscats.name_position_offset("0302+006", 0.62875, 0.100583) > 1e5


def test_position_quantum():
    assert lenscats.position_quantum_arcsec(12.345678, 1.234567) == 0  # 0.1'' precision
    # whole seconds of RA printed to 5 decimals (00:42:32, -17:13:45)
    q = lenscats.position_quantum_arcsec(10.63333, -17.22917)
    assert q == pytest.approx(15 * math.cos(math.radians(17.22917)), rel=1e-6)
    assert lenscats.position_quantum_arcsec(84.57, -49.48556) == pytest.approx(
        36 * math.cos(math.radians(49.48556))
    )
    assert lenscats.position_quantum_arcsec(175.8, 0.48) == pytest.approx(
        360 * math.cos(math.radians(0.48))
    )


def test_merge_keeps_unique_refs(tables):
    s = lenscats.merge(tables)
    row = s[np.abs(s["ra"] - 15) < 0.01][0]
    assert row["refs"] == "slacs | SuGOHI1"


def test_name_offset_catches_sign_errors():
    # declination sign flipped relative to the designation
    assert lenscats.name_position_offset("J0100+0030", 15.0, -0.51) > 3000
    assert lenscats.name_position_offset("J0100-0030", 15.0, -0.51) < 1


def test_position_quantum_needs_both_axes():
    # fine steps on one axis only (whole RA second, precise declination) -> not rounded
    assert lenscats.position_quantum_arcsec(10.63333, -17.229123) == 0
    assert lenscats.position_quantum_arcsec(84.57123, -49.481) == 0
    # a coarse 0.01 deg step on one axis is enough
    q = lenscats.position_quantum_arcsec(84.57, -49.48731)
    assert q == pytest.approx(36 * math.cos(math.radians(49.48731)))


def test_brick_coverage_independent_of_sources():
    bricks = Table(
        {
            "ra1": [0.0, 0.25, 359.75],
            "ra2": [0.25, 0.5, 360.0],
            "dec1": [-0.125, -0.125, -0.125],
            "dec2": [0.125, 0.125, 0.125],
            "nexp_r": [2, 0, 3],
            "nexp_z": [2, 4, 3],
            "galdepth_z": [23.4, 23.0, 23.1],
        }
    )
    systems = Table({"ra": [0.1, 0.3, 359.9, 10.0], "dec": [0.0, 0.0, 0.1, 0.0]})
    out = lenscats.brick_coverage(systems, bricks)
    assert list(out["covered"]) == [True, False, True, False]  # no r data; outside footprint
    assert out["depth_z"][0] == pytest.approx(23.4) and np.isnan(out["depth_z"][3])
