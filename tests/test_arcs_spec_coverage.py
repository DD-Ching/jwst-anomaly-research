import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import arcs_spec_coverage as sc  # noqa: E402


def test_polygon_contains_and_excludes():
    s = "POLYGON ICRS 10.0 -20.0 10.01 -20.0 10.01 -19.99 10.0 -19.99"
    assert sc.region_contains(s, 10.005, -19.995)
    assert not sc.region_contains(s, 10.02, -19.995)


def test_polygon_across_ra_zero_and_circle_and_multi():
    s = "POLYGON 359.995 0.0 0.005 0.0 0.005 0.01 359.995 0.01"
    assert sc.region_contains(s, 0.0, 0.005)
    assert sc.region_contains("CIRCLE ICRS 50.0 10.0 0.01", 50.005, 10.0)
    assert not sc.region_contains("CIRCLE ICRS 50.0 10.0 0.01", 50.02, 10.0)
    multi = "POLYGON 1 1 1.01 1 1.01 1.01 1 1.01 POLYGON 2 2 2.01 2 2.01 2.01 2 2.01"
    assert sc.region_contains(multi, 2.005, 2.005)
    assert not sc.region_contains("", 1.0, 1.0)


def test_mode_groups():
    assert sc.mode_group("ESO", "MUSE", "", "cube") == "ifu_slitless"
    assert sc.mode_group("ESO", "ALMA", "", "cube") == "other"
    assert sc.mode_group("ESO", "VIMOS", "", "spectrum") == "slit"
    assert sc.mode_group("ESO", "XSHOOTER", "", "spectrum") == "slit"
    assert sc.mode_group("MAST", "NIRISS/WFSS", "F150W;GR150R", "spectrum") == "ifu_slitless"
    assert sc.mode_group("MAST", "WFC3/IR", "G141", "spectrum") == "ifu_slitless"
    assert sc.mode_group("MAST", "NIRSPEC/MSA", "F170LP;G235M", "spectrum") == "slit"
    assert sc.mode_group("MAST", "NIRSPEC/IFU", "F170LP;G235H", "cube") == "ifu_slitless"


def test_classify_and_summary():
    assert sc.classify(["slit", "ifu_slitless"]) == "spectrum_available"
    assert sc.classify(["slit", "other"]) == "msa_footprint_only"
    assert sc.classify(["other"]) == "none"
    assert sc.classify([]) == "none"
    rows = [
        {
            "archive": "ESO",
            "instrument": "MUSE",
            "dataproduct": "cube",
            "contains": True,
            "public": False,
        },
        {
            "archive": "MAST",
            "instrument": "NIRSPEC/MSA",
            "mode": "x",
            "dataproduct": "spectrum",
            "contains": True,
            "public": True,
        },
        {
            "archive": "MAST",
            "instrument": "NIRISS/WFSS",
            "mode": "GR150C",
            "dataproduct": "spectrum",
            "contains": False,
            "public": True,
        },
    ]
    s = sc.summarize(rows)
    assert s["class"] == "msa_footprint_only"  # proprietary MUSE and non-covering WFSS do not count
    assert s["n_public_covering"] == 1


def test_gate_q():
    pos = {"by_instrument": {"ESO:MUSE:ifu_slitless": 3, "MAST:NIRISS/WFSS:ifu_slitless": 2}}
    assert sc.gate_q(pos, {"class": "none"})["passed"]
    assert not sc.gate_q(pos, {"class": "msa_footprint_only"})["passed"]
    assert not sc.gate_q({"by_instrument": {"ESO:MUSE:ifu_slitless": 1}}, {"class": "none"})[
        "passed"
    ]


def test_released_iso():
    assert sc.released_iso("2019-11-17T05:40:56Z", "2026-11-08T19:00:00")
    assert not sc.released_iso("2027-01-01T00:00:00Z", "2026-11-08T19:00:00")
    assert not sc.released_iso("--", "2026-11-08T19:00:00")
    assert not sc.released_iso("", "2026-11-08T19:00:00")
