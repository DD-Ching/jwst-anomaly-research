"""Unit tests for scripts/arc_diff_vet.py on synthetic arc scenes (offline)."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from astropy.wcs import WCS

_PATH = Path(__file__).resolve().parents[1] / "scripts" / "arc_diff_vet.py"
_spec = importlib.util.spec_from_file_location("arc_diff_vet", _PATH)
adv = importlib.util.module_from_spec(_spec)
sys.modules["arc_diff_vet"] = adv
_spec.loader.exec_module(adv)

N = 121
SCALE = 0.063  # arcsec per pixel


def _wcs_text() -> str:
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [64.0, -24.0]
    w.wcs.crpix = [N // 2 + 1, N // 2 + 1]
    w.wcs.cdelt = [-SCALE / 3600, SCALE / 3600]
    return w.to_header().tostring()


def _blob(y: float, x: float, sigma: float = 1.0, amp: float = 1.0) -> np.ndarray:
    yy, xx = np.mgrid[:N, :N]
    return amp * np.exp(-((yy - y) ** 2 + (xx - x) ** 2) / (2 * sigma**2))


def _arc() -> np.ndarray:
    """A vertical arc of knots through the centre, every 6 px (0.38″)."""
    return sum(_blob(N // 2 + k, N // 2) for k in range(-24, 25, 6)) + 0.05


def _maps(diff: np.ndarray, rng: np.random.Generator, band: str = "F277W") -> dict:
    d = diff + rng.normal(0, 0.01, (N, N))
    return {
        f"{band}_d": d,
        f"{band}_s": d / 0.01,
        f"{band}_e1img": _arc(),
        f"{band}_wcs": np.array(_wcs_text()),
    }


def _write(tmp_path: Path, name: str, maps: dict) -> Path:
    p = tmp_path / f"{name}.npz"
    np.savez(p, **maps)
    return p


def test_point_transient_is_not_arc_coherent(tmp_path):
    rng = np.random.default_rng(0)
    m = _write(tmp_path, "pt", _maps(_blob(N // 2, N // 2, amp=0.08), rng))
    out = tmp_path / "pt.json"
    adv.main(["--maps", "p", str(m), "--ra", "64.0", "--dec", "-24.0", "--out", str(out)])
    r = json.loads(out.read_text())["pairs"]["p"]
    band = r["bands"]["F277W"]
    assert band["s_hit"] > 5
    assert band["knots"], "neighbouring knots expected on the arc"
    assert band["coherent_fraction"] < 0.5 and not r["arc_coherent"]
    assert band["flux_njy"] > 3 * band["sky_sigma_njy"]


def test_arc_wide_scale_residual_is_arc_coherent(tmp_path):
    rng = np.random.default_rng(1)
    m = _write(tmp_path, "arc", _maps(0.06 * (_arc() - 0.05), rng))
    out = tmp_path / "arc.json"
    adv.main(["--maps", "p", str(m), "--ra", "64.0", "--dec", "-24.0", "--out", str(out)])
    r = json.loads(out.read_text())["pairs"]["p"]
    assert r["bands"]["F277W"]["coherent_fraction"] >= 0.5 and r["arc_coherent"]


def test_coherence_counts_same_sign_and_size():
    assert adv.coherence(6.0, [5.0, 3.5, -4.0, 1.0], 0.5) == 0.5
    assert np.isnan(adv.coherence(6.0, [], 0.5))
    assert adv.coherence(-4.0, [-3.0, -1.0], 0.5) == 0.5


def test_position_off_the_map_or_all_nan_is_skipped(tmp_path):
    rng = np.random.default_rng(2)
    maps = _maps(_blob(N // 2, N // 2, amp=0.08), rng)
    maps["F277W_s"][N // 2 - 3 : N // 2 + 4, N // 2 - 3 : N // 2 + 4] = np.nan
    m = _write(tmp_path, "nan", maps)
    out = tmp_path / "nan.json"
    adv.main(["--maps", "p", str(m), "--ra", "64.0", "--dec", "-24.0", "--out", str(out)])
    assert json.loads(out.read_text())["pairs"] == {}
    adv.main(["--maps", "p", str(m), "--ra", "64.1", "--dec", "-24.0", "--out", str(out)])
    assert json.loads(out.read_text())["pairs"] == {}
