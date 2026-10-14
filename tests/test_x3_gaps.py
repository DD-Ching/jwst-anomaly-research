"""Offline tests for scripts/x3_gaps.py (X3 G0: W3 in unselected Gaia epoch photometry)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

_P = Path(__file__).resolve().parents[1] / "scripts" / "x3_gaps.py"
_spec = importlib.util.spec_from_file_location("x3_gaps", _P)
x3 = importlib.util.module_from_spec(_spec)
sys.modules.setdefault("x3_gaps", x3)
_spec.loader.exec_module(x3)


def test_gaps_pixels_cover_the_cone():
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astropy_healpix import HEALPix

    pix = set(x3.gaps_pixels())
    hp = HEALPix(nside=2**x3.HPX_LEVEL, order="nested")
    c = SkyCoord(x3.M31[0] * u.deg, x3.M31[1] * u.deg)
    rng = np.random.default_rng(1)
    pa = rng.uniform(0, 360, 400) * u.deg
    sep = np.sqrt(rng.uniform(0, 1, 400)) * 5.5 * u.deg
    pts = c.directional_offset_by(pa, sep)
    assert set(int(p) for p in hp.lonlat_to_healpix(pts.ra, pts.dec)) <= pix


def test_source_id_ranges_follow_gaia_healpix_encoding():
    # Gaia: HEALPix level-12 nested index = source_id // 2**35
    shift = 35 + 2 * (12 - x3.HPX_LEVEL)
    sid = 375970751214564608  # a GAPS source
    lvl12 = sid >> 35
    assert sid >> shift == lvl12 >> (2 * (12 - x3.HPX_LEVEL))
    assert "MOD(source_id, 10) = 0" in x3.DEV_QUERY  # frozen dev split


def _tab(kind, flagged, te=np.nan, rho=np.nan):
    n = len(kind)
    return Table(
        {
            "event_id": [str(i) for i in range(n)],
            "kind": kind,
            "flagged": flagged,
            "tE": np.full(n, te) if np.isscalar(te) else te,
            "rho": np.full(n, rho) if np.isscalar(rho) else rho,
            "error": [""] * n,
        }
    )


def test_gate_verdicts():
    from astropy.table import vstack

    real_ok = _tab(["real"] * 200, [False] * 199 + [True])  # 0.5 %
    real_bad = _tab(["real"] * 100, [False] * 90 + [True] * 10)  # 10 %
    good = [
        _tab(["W3"] * 10, [True] * 5 + [False] * 5, te, 0.1) for te in (10.0, 30.0, 100.0, 300.0)
    ]
    dead = [_tab(["W3"] * 10, [False] * 10, te, 0.1) for te in (10.0, 30.0, 100.0, 300.0)]
    assert x3.summarise(vstack([real_ok, *good]))["gate"] == "continue to G1"
    assert x3.summarise(vstack([real_ok, *dead]))["gate"].startswith("stop")
    assert x3.summarise(vstack([real_bad, *good]))["gate"].startswith("one recorded")


def test_v1_veto_counts_vetoed_injections_as_missed():
    from astropy.table import vstack

    real = _tab(["real"] * 100, [False] * 98 + [True] * 2)
    real["event_id"] = [str(i) for i in range(100)]
    w3 = [_tab(["W3"] * 10, [True] * 10, te, 0.1) for te in (10.0, 30.0, 100.0, 300.0)]
    for t in w3:
        t["event_id"] = [str(i) for i in range(10)]
    tab = vstack([real, *w3])
    veto = {str(i): i < 5 or i in (98, 99) for i in range(100)}  # vetoes both real flags
    s = x3.summarise(tab, veto)
    assert s["flag_fraction_real"] == 0.0
    assert abs(s["mean_eff_te30_300"] - 0.5) < 1e-12  # half the injection bases vetoed
    assert abs(s["v1_vetoed_fraction"] - 0.07) < 1e-12
