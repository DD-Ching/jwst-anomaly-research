"""Frozen outcome rule and arc model of scripts/m2_embeddings.py (offline)."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import m2_embeddings as m  # noqa: E402


def test_decide_beats_needs_margin_and_every_seed():
    assert m.decide([0.5, 0.5, 0.5], [0.3, 0.3, 0.3]) == "beats"
    # Mean margin met but one seed below the baseline: tie.
    assert m.decide([0.9, 0.6, 0.2], [0.3, 0.3, 0.3]) == "tie"
    assert m.decide([0.35, 0.35, 0.35], [0.3, 0.3, 0.3]) == "tie"


def test_decide_worse_and_errors():
    assert m.decide([0.1, 0.1, 0.1], [0.3, 0.3, 0.3]) == "worse"
    assert m.decide([0.2, 0.2, 0.2], [0.3, 0.3, 0.3]) == "worse"  # exactly at the margin
    with pytest.raises(ValueError):
        m.decide([], [])
    with pytest.raises(ValueError):
        m.decide([0.1], [0.1, 0.2])


def test_arc_model_flux_and_radius():
    rng = np.random.default_rng(1)
    img, par = m.arc_model(rng, 64, 0.0312, 10.0)
    assert img.sum() == pytest.approx(10.0)
    assert 0.35 / 0.0312 <= par["r_px"] <= 0.70 / 0.0312
    assert 90 <= par["extent_deg"] <= 180
    c = 31.5
    yy, xx = np.mgrid[:64, :64]
    rr = np.hypot(yy - c, xx - c)
    w = img / img.sum()
    assert float((rr * w).sum()) == pytest.approx(par["r_px"], abs=1.0)
