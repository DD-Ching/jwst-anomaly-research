import sys
from datetime import date
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import sneos_lead_epochs as m  # noqa: E402


def test_band_lya_cut():
    assert m.LYA_A * (1 + m.Z_SN) == pytest.approx(7457.7, abs=0.1)
    assert not m.rule_band("FR505N")
    assert not m.rule_band("F606W")
    assert m.rule_band("F814W")
    assert m.rule_band("F110W")
    with pytest.raises(KeyError):
        m.rule_band("F999X")


def test_phase_window_with_allowance():
    # 1000 d lead: t' range [t + 800, t + 1200] d
    assert m.rule_phase(date(2022, 1, 1), -1000.0)  # t' 2024-03-11 .. 2025-04-15 overlaps
    assert not m.rule_phase(date(2019, 1, 1), -1000.0)  # t' up to 2022-04 < window start
    assert not m.rule_phase(date(2024, 3, 11), -1000.0)  # t' from 2026-05-20 > window end
    assert m.rule_phase(date(2024, 3, 11), -1000.0, m.WINDOW_EXT)
    # edge: hi exactly at window start counts as overlap
    assert m.rule_phase(date(2024, 3, 11), 0.0)


def test_reference_and_classify():
    assert m.rule_reference("F110W") and not m.rule_reference("F2100W")
    r = m.classify(date(2022, 1, 1), "F110W", -1000.0)
    assert r["pass"]
    r = m.classify(date(2022, 1, 1), "F2100W", -1000.0)
    assert r["band"] and r["phase"] and not r["reference"] and not r["pass"]
    assert m.outcome([r]) == "U"
    assert m.outcome([r, m.classify(date(2022, 1, 1), "F110W", -1000.0)]) == "P"


def test_parse_filters_and_unknown_filter():
    assert m.parse_filters("F110W;CLEAR") == ["F110W"]
    assert m.parse_filters("CLEAR;F444W") == ["F444W"]
    assert m.parse_filters("F444W;F405N") == ["F444W", "F405N"]
    assert m.parse_filters("detection") == []
    r = m.classify(date(2022, 1, 1), "F405N", -1000.0)
    assert not r["band_known"] and not r["band"] and not r["pass"]
