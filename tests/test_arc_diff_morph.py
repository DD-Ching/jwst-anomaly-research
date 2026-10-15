"""Unit tests for scripts/arc_diff_morph.py on synthetic significance maps (offline)."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.wcs import WCS

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("arc_diff_morph", _DIR / "arc_diff_morph.py")
adm = importlib.util.module_from_spec(_spec)
sys.modules["arc_diff_morph"] = adm
_spec.loader.exec_module(adm)

N = 121
SCALE = 0.063  # arcsec per pixel
BANDS = ("F277W", "F356W")


def _wcs_text() -> str:
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [64.0, -24.0]
    w.wcs.crpix = [N // 2 + 1, N // 2 + 1]
    w.wcs.cdelt = [-SCALE / 3600, SCALE / 3600]
    return w.to_header().tostring()


def _blob(y: float, x: float, sigma: float = 1.5, amp: float = 1.0) -> np.ndarray:
    yy, xx = np.mgrid[:N, :N]
    return amp * np.exp(-((yy - y) ** 2 + (xx - x) ** 2) / (2 * sigma**2))


def _run(tmp_path: Path, s: np.ndarray, name: str) -> dict:
    rng = np.random.default_rng(0)
    maps = {}
    for b in BANDS:
        sb = s + rng.normal(0, 1, (N, N))
        sb[np.isnan(s)] = np.nan
        maps |= {f"{b}_s": sb, f"{b}_wcs": np.array(_wcs_text())}
    npz = tmp_path / f"{name}.npz"
    np.savez(npz, **maps)
    hits = tmp_path / f"{name}.ecsv"
    Table({"hit_id": ["h000"], "ra": [64.0], "dec": [-24.0], "near_bright_core": [False]}).write(
        hits
    )
    out = tmp_path / f"{name}.json"
    adm.main(["--maps", str(npz), "--hits", str(hits), "--out", str(out)])
    return json.loads(out.read_text())["targets"][0]


def test_point_transient_passes(tmp_path):
    t = _run(tmp_path, _blob(N // 2, N // 2, amp=12.0), "pt")
    assert t["verdict"] == "pass"
    assert t["bands"]["F277W"]["extent_arcsec"] < 1.0


def test_spike_streak_is_extended(tmp_path):
    yy, xx = np.mgrid[:N, :N]
    s = -8.0 * np.exp(-((yy - xx) ** 2) / (2 * 1.5**2))  # diagonal line through the centre
    t = _run(tmp_path, s, "streak")
    assert t["verdict"] == "extended"
    assert t["bands"]["F277W"]["extent_arcsec"] > 3.0


def test_alternating_lobes_are_a_rosette(tmp_path):
    c = N // 2
    s = _blob(c, c, amp=8.0)
    for k in range(6):
        phi = k * np.pi / 3
        s += _blob(c + 7 * np.sin(phi), c + 7 * np.cos(phi), amp=8.0 * (-1) ** (k + 1))
    t = _run(tmp_path, s, "rosette")
    assert t["verdict"] == "rosette"
    assert t["bands"]["F356W"]["n_lobes"] >= 4 and t["bands"]["F356W"]["both_signs"]


def test_core_dipole_is_flagged(tmp_path):
    c = N // 2
    s = _blob(c, c, amp=12.0) + _blob(c, c + 3, amp=-10.0)
    t = _run(tmp_path, s, "dipole")
    assert t["verdict"] == "dipole"
    assert t["bands"]["F277W"]["dipole_ratio"] > 0.5


def test_brightening_beside_a_weak_negative_lobe_passes(tmp_path):
    c = N // 2
    s = _blob(c, c, amp=12.0) + _blob(c, c + 3, amp=-4.0)
    t = _run(tmp_path, s, "weakdip")
    assert t["verdict"] == "pass"
    assert t["bands"]["F277W"]["dipole_ratio"] < 0.5


def test_hit_next_to_missing_data_is_edge(tmp_path):
    s = _blob(N // 2, N // 2, amp=12.0)
    s[:, N // 2 + 4 :] = np.nan
    t = _run(tmp_path, s, "edge")
    assert t["verdict"] == "edge"
    assert t["bands"]["F277W"]["edge_arcsec"] < 0.5


def test_component_extent_of_a_line_is_its_length():
    m = np.zeros((50, 50), bool)
    m[25, 5:45] = True
    assert abs(adm.component_extent(m, 25, 20) - 40) < 1.0
    assert adm.component_extent(m, 10, 10) == 0.0
