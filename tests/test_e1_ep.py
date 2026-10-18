from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
_spec = importlib.util.spec_from_file_location("e1_ep", ROOT / "scripts" / "e1_ep.py")
ep = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ep)

BODY_FIRST = (
    "We report on the detection of an X-ray transient by the Wide-field X-ray Telescope (WXT) on "
    "board the Einstein Probe (EP) mission, designated EP260601a. The transient triggered EP-WXT "
    "(ID: 01709266307) at 2026-06-01T19:13:22 UTC. The WXT position of the source is R.A. = "
    "257.649 deg, Dec = -1.647 deg (J2000) with an uncertainty of 3 arcmin in radius (90% C.L.)."
)
BODY_ONSET = (
    "The EP-WXT observation (ID: 11900946826) started at 2026-10-15 17:29:32 (UTC), with an "
    "exposure time of 1200 s. The WXT light curve of the transient began at 2026-10-15T17:37:52 "
    "UTC lasted for 700 s. The WXT position of the source is R.A. = 21.920 deg, DEC = 66.601 deg "
    "(J2000) with an uncertainty of 2.695 arcmin in radius. EP-FXT gives R.A. = 21.93, "
    "DEC = 66.60 with an uncertainty of 10 arcsec."
)


def test_parse_first_timestamp_and_negative_dec():
    r = ep.parse_circular("EP260601a: EP-WXT detection of an X-ray transient", BODY_FIRST)
    assert r["name"] == "EP260601a" and r["time_rule"] == "first"
    assert abs(r["mjd"] - (61192 + (19 * 3600 + 13 * 60 + 22) / 86400)) * 86400 < 0.01
    assert (r["ra"], r["dec"]) == (257.649, -1.647)
    assert np.isclose(r["r90_deg"], 3 / 60) and r["r90_unit"] == "arcmin"


def test_parse_onset_wins_and_arcmin_radius_preferred():
    r = ep.parse_circular(
        "EP261015c: Einstein Probe detection of a fast X-ray transient", BODY_ONSET
    )
    assert r["time_rule"] == "onset"
    assert abs(r["mjd"] - (61328 + (17 * 3600 + 37 * 60 + 52) / 86400)) * 86400 < 0.01
    assert np.isclose(r["r90_deg"], 2.695 / 60)


def test_parse_rejects_name_date_mismatch_and_no_wxt():
    assert ep.parse_circular("EP260501a: EP-WXT detection", BODY_FIRST) is None
    assert ep.parse_circular("EP260601a: detection", BODY_FIRST.replace("WXT", "XRT")) is None


def test_wxt_must_precede_position():
    body = (
        "FXT position R.A. = 257.649 deg, Dec = 1.6 deg with an uncertainty of 10 arcsec at "
        "2026-06-01T19:13:22. WXT also."
    )
    assert ep.parse_circular("EP260601a: EP detection", body) is None


def test_retraction_only_of_the_trigger():
    assert ep.is_retraction(
        "The EP-WXT trigger 01709133021 (designated as EP250324a) is not a real source"
    )
    assert ep.is_retraction("EP260707a is not a real source\nH. Zhou (PMO)")
    assert not ep.is_retraction("The optical candidate is not a real source.")
    assert not ep.is_retraction("The optical candidate of EP250903a is not a real source.")
    assert not ep.is_retraction("EP250903a: the afterglow candidate is not a real source.")


def test_ep_events_lowest_circular_and_retraction():
    circ = [
        {
            "circularId": 20,
            "subject": "EP260601a: follow-up",
            "body": BODY_FIRST.replace("257.649", "1.0"),
        },
        {"circularId": 10, "subject": "EP260601a: EP-WXT detection", "body": BODY_FIRST},
        {"circularId": 30, "subject": "EP261015c: EP-WXT detection", "body": BODY_ONSET},
        {
            "circularId": 31,
            "subject": "EP261015c is not a real source",
            "body": "The trigger is not a real source.",
        },
    ]
    t = ep.ep_events(circ)
    assert list(t["name"]) == ["EP260601a"]
    assert t["circular"][0] == 10 and t["ra"][0] == 257.649
    assert np.isclose(t["sigma"][0], 3 / 60 / 2.146)


def test_mask_family_has_18_cells():
    m = ep.test_mask(ep.P)
    assert int(m.sum()) == 18
    epep = ep.CHANNELS.index(("EP", "EP"))
    assert m[epep, :5, :2].all() and not m[epep, :, 2].any()
    for ch in ep.CONTROL:
        c = ep.CHANNELS.index(ch)
        assert not m[c, :3].any() and m[c, 3:5, :2].all()


def test_control_sum_and_rule():
    counts = np.zeros((len(ep.CHANNELS), 8, 3))
    for ch in ep.CONTROL:
        c = ep.CHANNELS.index(ch)
        counts[c, :3, 0] = [1, 2, 3]
        counts[c, 3, 0] = 100  # 1 h-1 d: not in the control
    counts[ep.CHANNELS.index(("EP", "EP")), 0, 0] = 50
    assert ep.control_sum(counts) == 12
    ens = np.stack([counts, 2 * counts, 3 * counts])  # an ensemble keeps its scramble axis
    assert ep.control_sum(ens).tolist() == [12, 24, 36]
    null = np.random.default_rng(0).poisson(0.1, 10000)
    ok, p = ep.control_passes(12, null)
    assert ok and p < 1e-10
    assert not ep.control_passes(4, null)[0]  # too few pairs, however significant


def test_outcome_rule():
    assert ep.outcome(np.array([0.5, 1e-3]), 18) == "null"
    assert ep.outcome(np.array([0.5, 1.3e-3 / 18]), 18) == "excess"
    assert ep.outcome(np.array([0.5, 1.4e-3 / 18]), 18) == "null"
