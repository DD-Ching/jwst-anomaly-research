"""Unit tests for scripts/arc_diff_knots.py on synthetic knot fields (offline)."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from astropy.wcs import WCS

_DIR = Path(__file__).resolve().parents[1] / "scripts"
for _name in ("arc_diff_vet", "arc_diff_knots"):
    _spec = importlib.util.spec_from_file_location(_name, _DIR / f"{_name}.py")
    _mod = importlib.util.module_from_spec(_spec)
    sys.modules[_name] = _mod
    _spec.loader.exec_module(_mod)
adk = sys.modules["arc_diff_knots"]

N = 121
SCALE = 0.063  # arcsec per pixel


def _wcs() -> WCS:
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [64.0, -24.0]
    w.wcs.crpix = [N // 2 + 1, N // 2 + 1]
    w.wcs.cdelt = [-SCALE / 3600, SCALE / 3600]
    return w


def _blob(y: float, x: float, amp: float, sigma: float = 1.0) -> np.ndarray:
    yy, xx = np.mgrid[:N, :N]
    return amp * np.exp(-((yy - y) ** 2 + (xx - x) ** 2) / (2 * sigma**2))


def _scene(tmp_path: Path, deficit: float, extra: float) -> tuple[Path, float, float]:
    """A 10 × 10 grid of knots; ``d`` = deficit × e1 plus ``extra`` × the centre knot."""
    rng = np.random.default_rng(1)
    centres = [(y, x) for y in range(15, 110, 10) for x in range(15, 110, 10)]
    amps = rng.uniform(1000, 3000, len(centres))
    e1 = sum(_blob(y, x, a) for (y, x), a in zip(centres, amps, strict=True))
    d = deficit * e1 + extra * _blob(*centres[55], amps[55])
    d = d + rng.normal(0, 0.05, (N, N))
    p = tmp_path / "maps.npz"
    np.savez(p, F150W_d=d, F150W_e1img=e1, F150W_wcs=np.array(_wcs().to_header().tostring()))
    ra, dec = (float(v) for v in _wcs().all_pix2world(centres[55][1], centres[55][0], 0))
    return p, ra, dec


def _run(tmp_path: Path, deficit: float, extra: float) -> dict:
    maps, ra, dec = _scene(tmp_path, deficit, extra)
    out = tmp_path / "k.json"
    args = [str(maps), "--pair", "t", "--bands", "F150W", "--target", "x", str(ra), str(dec)]
    adk.main(args + ["--min-njy", "50", "--max-njy", "1e9", "--out", str(out)])
    return json.loads(out.read_text())["bands"]["F150W"]


def test_uniform_deficit_is_the_knot_median_and_target_is_ordinary(tmp_path):
    r = _run(tmp_path, 0.02, 0.0)
    assert r["n_knots"] == 100
    assert abs(r["median_frac_comparable"] - 0.02) < 0.002
    assert abs(r["targets"]["x"]["z_vs_knots"]) < 3


def test_extra_brightening_stands_out_of_the_knots(tmp_path):
    r = _run(tmp_path, 0.02, 0.10)
    t = r["targets"]["x"]
    assert abs(t["frac"] - 0.12) < 0.01
    assert t["rank_from_top"] == 1
    assert t["z_vs_knots"] > 10


def test_robust_z_needs_five_references():
    assert np.isnan(adk.robust_z(1.0, np.array([0.0, 1.0]))[2])
    med, sig, z = adk.robust_z(3.0, np.array([0.0, 1.0, -1.0, 2.0, -2.0]))
    assert med == 0.0 and sig > 0 and z > 0


def test_holed_knot_excluded_and_edge_target_recorded(tmp_path):
    maps, ra, dec = _scene(tmp_path, 0.02, 0.0)
    with np.load(maps) as m:
        arrs = dict(m)
    arrs["F150W_d"][15, 15] = np.nan  # inside the aperture of the corner knot at (15, 15)
    np.savez(maps, **arrs)
    out = tmp_path / "k.json"
    ra_edge, dec_edge = (float(v) for v in _wcs().all_pix2world(0, 0, 0))
    adk.main(
        [str(maps), "--pair", "t", "--bands", "F150W", "--target", "x", str(ra), str(dec)]
        + ["--target", "edge", str(ra_edge), str(dec_edge)]
        + ["--min-njy", "50", "--max-njy", "1e9", "--out", str(out)]
    )
    r = json.loads(out.read_text())["bands"]["F150W"]
    assert r["n_knots"] == 99
    assert r["targets"]["edge"] == {"skipped": "outside the map interior"}
    assert r["targets"]["x"]["aperture_has_nan"] is False


def test_finite_or_none_gives_strict_json():
    obj = {"a": float("nan"), "b": [1.0, np.float64("inf")], "c": "s"}
    assert json.loads(json.dumps(adk.finite_or_none(obj), allow_nan=False)) == {
        "a": None,
        "b": [1.0, None],
        "c": "s",
    }
