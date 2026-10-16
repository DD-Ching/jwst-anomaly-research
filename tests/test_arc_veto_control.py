import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import arc_veto_control as avc  # noqa: E402


def test_has_negative_lobe():
    assert avc.has_negative_lobe({"F150W": {"lobe_peaks": [5.0, -3.0]}})
    assert not avc.has_negative_lobe({"F150W": {"lobe_peaks": [5.0, -2.9]}, "F444W": {}})


def test_decide_outcomes():
    assert avc.decide([True] * 8 + [False] * 2, [True] * 3 + [False] * 7)["outcome"] == "H1"
    # vetoed fraction high but no excess over passed -> mixed
    assert avc.decide([True] * 8 + [False] * 2, [True] * 6 + [False] * 4)["outcome"] == "mixed"
    assert avc.decide([True] * 3 + [False] * 7, [False] * 5)["outcome"] == "H2"
    assert avc.decide([True] * 5 + [False] * 5, [False] * 5)["outcome"] == "mixed"
    assert avc.decide([], [True])["outcome"] == "undetermined"


def test_injected_hits_and_compare(tmp_path):
    run = tmp_path / "inj0"
    run.mkdir()
    (run / "injected.csv").write_text(
        "name,ra,dec,mag,hit,near_bright_core\n"
        "i0,1.0,2.0,27.5,h1,False\ni1,1.1,2.0,27.5,h2,False\n"
        "i2,1.2,2.0,27.5,,False\ni3,1.3,2.0,27.5,h3,True\n"
    )
    morph = {
        "targets": [
            {"id": "h1", "verdict": "rosette", "bands": {"F150W": {"s_hit": 8.0}}},
            {"id": "h2", "verdict": "pass", "bands": {"F150W": {"s_hit": -6.0}}},
        ]
    }
    (run / "morph.json").write_text(json.dumps(morph))
    rows = avc.injected_hits(run)
    assert [r["name"] for r in rows] == ["inj0:i0", "inj0:i1"]
    ctrl = {
        "targets": [
            {"id": "inj0:i0", "verdict": "pass", "bands": {"F150W": {"lobe_peaks": [-4.0]}}},
            {"id": "inj0:i1", "verdict": "pass", "bands": {"F150W": {"lobe_peaks": [4.0]}}},
        ]
    }
    cp = tmp_path / "ctrl.json"
    cp.write_text(json.dumps(ctrl))
    out = avc.compare([run], cp, None, None)
    s = out["strata"]["27.5/all"]
    assert (s["n_vetoed"], s["n_passed"], s["frac_vetoed"], s["frac_passed"]) == (1, 1, 1.0, 0.0)
    assert s["outcome"] == "H1"
    assert s["post_hoc_neg_ratio_vetoed"] == [0.5]
    assert s["post_hoc_neg_ratio_passed"] == [0.0]


def test_negative_ratio():
    assert avc.negative_ratio({"a": {"lobe_peaks": [-3.0, -6.0, 9.0]}}, 12.0) == 0.5
    assert avc.negative_ratio({"a": {"lobe_peaks": [-2.0]}}, 12.0) == 0.0
