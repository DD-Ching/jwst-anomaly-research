"""Tests for scripts/fade_calib.py and arc_diff pre-matching injection (synthetic, offline)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("fade_calib", _DIR / "fade_calib.py")
fc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fc)


def test_class_mags_and_cycle():
    assert fc.class_mags("vanish", 26.0) == (26.0, None)
    assert fc.class_mags("appear", 26.0) == (None, 26.0)
    assert fc.class_mags("partial", 26.0) == (26.0, 26.75)
    rows = fc.sample_rows([(1.0, 2.0)] * 12, [26.0, 27.0, 27.5])
    cells = {(r["cls"], r["mag"]) for r in rows}
    assert len(cells) == 12  # 4 classes x 3 mags each appear once in 12 rows


def test_sign_ok():
    assert fc.sign_ok("vanish", "-", {})
    assert fc.sign_ok("vanish", "mixed", {"F150W": -8.0, "F444W": 3.0})
    assert not fc.sign_ok("vanish", "mixed", {"F150W": -3.0, "F444W": 8.0})
    assert not fc.sign_ok("vanish", "+", {})
    assert fc.sign_ok("appear", "+", {})


def test_summarize_and_pass_rule():
    rows = []
    for i in range(10):
        ok = i < 8
        rows.append(
            {
                "cls": "vanish",
                "mag": "26.00",
                "hit": ok,
                "unflagged": ok,
                "morph_pass": ok,
                "sign_ok": ok,
                "review": ok,
                "review_sign": ok,
            }
        )
        rows.append(
            {
                "cls": "unchanged",
                "mag": "26.00",
                "hit": i == 0,
                "unflagged": False,
                "morph_pass": False,
                "sign_ok": False,
                "review": False,
                "review_sign": False,
            }
        )
    s = fc.summarize(rows)
    assert s["vanish"]["26.00"]["review_sign_frac"] == 0.8
    assert s["unchanged_false_hits"]["frac"] == 0.1
    assert fc.pass_rule(s)["passed"]
    lo, hi = fc.wilson(8, 10)
    assert lo < 0.8 < hi


def test_add_star_flux_and_mag_scale():
    import arc_diff as ad

    img = np.zeros((41, 41))
    ad.add_star(img, 20.0, 20.0, 1.5, 3.0)
    assert np.isclose(img.sum(), 3.0, rtol=1e-6)
    f1 = ad.mag_to_mjysr(26.0, 2.3e-14)
    f2 = ad.mag_to_mjysr(26.75, 2.3e-14)
    assert np.isclose(f2 / f1, 10 ** (-0.4 * 0.75))
