"""Unit tests for scripts/arc_diff_morph.py on synthetic significance maps (offline)."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table
from astropy.wcs import WCS
from scipy import ndimage

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


def _noise(rng: np.random.Generator, corr_px: float) -> np.ndarray:
    """Unit-variance noise; corr_px > 0 mimics a matched-filter S map (correlated noise)."""
    n = rng.normal(0, 1, (N, N))
    if corr_px > 0:
        n = ndimage.gaussian_filter(n, corr_px)
        n /= n.std()
    return n


def _run(tmp_path: Path, s: np.ndarray, name: str, corr_px: float = 0.0) -> dict:
    rng = np.random.default_rng(0)
    maps = {}
    for b in BANDS:
        sb = s + _noise(rng, corr_px)
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


def test_point_transient_passes_in_correlated_noise(tmp_path):
    # Matched-filter S has correlated noise; a white-noise rescale let it percolate into "extended".
    t = _run(tmp_path, _blob(N // 2, N // 2, amp=12.0), "ptc", corr_px=1.5)
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


def test_dipole_ratio_ignores_a_hit_below_the_level():
    s = np.zeros((21, 21))
    s[10, 10] = 2.0
    s[10, 12:15] = -5.0
    s[11, 12:15] = -5.0
    assert adm.dipole_ratio(s, 10, 10, 4.0, 3.0, 4) == 0.0
    s[10, 10] = 6.0
    assert adm.dipole_ratio(s, 10, 10, 4.0, 3.0, 4) > 0.5


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


def _band(s_hit, extent=0.0, n_lobes=0, both=False, dip=0.0, edge=5.0) -> dict:
    return {
        "s_hit": s_hit,
        "extent_arcsec": extent,
        "n_lobes": n_lobes,
        "both_signs": both,
        "dipole_ratio": dip,
        "edge_arcsec": edge,
    }


def test_classify_extent_counts_every_detected_band():
    # A spike broken up in the top band but continuous in another detected band is still extended.
    p = adm.Params()
    res = {"F150W": _band(20.0, extent=0.4), "F200W": _band(8.0, extent=9.0)}
    assert adm.classify(res, p) == ("extended", "F150W")
    res["F200W"]["s_hit"] = 2.0  # not detected there: its extent no longer counts
    assert adm.classify(res, p) == ("pass", "F150W")


def test_classify_ignores_lobes_and_edges_of_undetected_bands():
    p = adm.Params()
    noisy = _band(1.5, n_lobes=6, both=True, edge=0.1)
    res = {"F150W": _band(12.0), "F277W": noisy, "F356W": dict(noisy)}
    assert adm.classify(res, p)[0] == "pass"
    res["F277W"]["s_hit"] = res["F356W"]["s_hit"] = 4.0
    assert adm.classify(res, p)[0] == "edge"


def test_classify_without_a_detected_band_uses_the_top_band():
    p = adm.Params()
    res = {"F150W": _band(2.5, edge=0.2), "F200W": _band(1.0)}
    assert adm.classify(res, p) == ("edge", "F150W")


def test_dipole_ratio_uses_the_local_region_peak_not_a_far_streak_peak():
    s = np.zeros((41, 41))
    s[20, 20] = 10.0
    s[20:22, 23:40] = -4.0  # touches the r = 3.5 disk only with its -4 end
    s[20, 25] = -9.0  # a lobe peak just outside the disk: counts
    assert abs(adm.dipole_ratio(s, 20, 20, 3.5, 3.0, 4) - 0.9) < 1e-9
    s[20, 25] = -4.0
    s[20, 38] = -12.0  # the far end of a long streak: does not count
    assert abs(adm.dipole_ratio(s, 20, 20, 3.5, 3.0, 4) - 0.4) < 1e-9


def test_rejudge_refuses_other_measurement_thresholds(tmp_path):
    path = tmp_path / "m.json"
    path.write_text(json.dumps({"thresholds_assumption": {"level": 4.0}, "targets": []}))
    with pytest.raises(ValueError, match="level"):
        adm.rejudge(path, adm.Params())


def test_hit_at_the_map_border_is_edge_not_dropped(tmp_path):
    w = WCS(_wcs_text())
    s = _blob(N // 2, 1, amp=12.0)
    rng = np.random.default_rng(1)
    maps = {}
    for b in BANDS:
        maps |= {f"{b}_s": s + _noise(rng, 0.0), f"{b}_wcs": np.array(_wcs_text())}
    npz = tmp_path / "border.npz"
    np.savez(npz, **maps)
    ra, dec = (float(v) for v in w.all_pix2world(0.0, N // 2, 0))
    hits = tmp_path / "border.ecsv"
    Table({"hit_id": ["h000"], "ra": [ra], "dec": [dec]}).write(hits)
    out = tmp_path / "border.json"
    adm.main(["--maps", str(npz), "--hits", str(hits), "--out", str(out)])
    assert json.loads(out.read_text())["targets"][0]["verdict"] == "edge"


def test_rejudge_reclassifies_stored_values(tmp_path):
    rows = [
        {
            "id": "a",
            "verdict": "rosette",
            "top_band": "F150W",
            "bands": {
                "F150W": _band(12.0),
                "F277W": _band(1.0, n_lobes=5, both=True),
                "F356W": _band(2.0, n_lobes=5, both=True),
            },
        },
    ]
    path = tmp_path / "m.json"
    path.write_text(json.dumps({"counts": {"rosette": 1}, "targets": rows}))
    rep = adm.rejudge(path, adm.Params())
    assert rep["changed"] == [{"id": "a", "old": "rosette", "new": "pass"}]
    d = json.loads(path.read_text())
    assert d["counts"] == {"pass": 1} and d["rejudged"] and d["dipole_measured"]


def test_rosette_lobe_ratio_counts_only_strong_lobes():
    # A bright star beside weak residual lobes of a neighbour is not a rosette under the ratio rule.
    weak = _band(20.0, n_lobes=5, both=True)
    weak["lobe_peaks"] = [20.0, 3.5, -3.2, 4.0, -3.8]
    res = {"F277W": weak, "F356W": dict(weak)}
    assert adm.classify(res, adm.Params(rosette_lobe_ratio=0.0))[0] == "rosette"
    assert adm.classify(res, adm.Params(rosette_lobe_ratio=0.3))[0] == "pass"
    strong = _band(8.0, n_lobes=4, both=True)
    strong["lobe_peaks"] = [8.0, -7.0, 6.5, -5.0]
    res = {"F277W": strong, "F356W": dict(strong)}
    assert adm.classify(res, adm.Params(rosette_lobe_ratio=0.3))[0] == "rosette"
    old = _band(8.0, n_lobes=4, both=True)  # outputs without lobe_peaks keep the original rule
    assert adm.classify({"A": old, "B": dict(old)}, adm.Params(rosette_lobe_ratio=0.3))[0] == (
        "rosette"
    )


def test_lobe_peaks_are_signed_region_peaks():
    s = np.zeros((40, 40))
    s[18:23, 18:23] = 9.0
    s[10:13, 18:23] = -4.0
    pk = sorted(adm.lobe_peaks(s, 20, 20, 12, 3.0, 4))
    assert pk == [-4.0, 9.0]


def test_strong_opposite_signs_between_bands_are_a_sign_flip():
    # PSF-mismatch rings flip sign between SW and LW; a new point source changes every band one way.
    res = {"F150W": _band(-6.0), "F200W": _band(-4.0), "F356W": _band(7.0)}
    assert adm.classify(res, adm.Params())[0] == "sign_flip"
    assert adm.classify(res, adm.Params(sign_flip_level=0.0))[0] == "pass"
    # Yan+2023 D21-S1 in o211/o212: a weak opposite-sign band (-3.2) does not veto it.
    res = {"F200W": _band(-3.2), "F277W": _band(6.9), "F356W": _band(6.7)}
    assert adm.classify(res, adm.Params())[0] == "pass"


def _ring_maps(with_ring: bool):
    """Epoch-1 galaxy at the centre; S has a point source 0.6″ east and optionally a full ring."""
    yy, xx = np.mgrid[:101, :101]
    e1 = 100.0 * np.exp(-(np.hypot(yy - 50, xx - 50) ** 2) / 8.0)
    s = np.zeros((101, 101))
    r = np.hypot(yy - 50, xx - 50)
    if with_ring:
        s += 6.0 * np.exp(-((r - 20.0) ** 2) / 2.0)
    s += 8.0 * np.exp(-((yy - 50) ** 2 + (xx - 70) ** 2) / 4.0)
    return s, e1


def test_ring_fill_separates_a_galaxy_ring_from_a_point_source():
    p = adm.Params()
    s, e1 = _ring_maps(with_ring=True)
    assert adm.ring_fill(s, e1, 50, 70, 0.03, p) == 1.0
    s, e1 = _ring_maps(with_ring=False)
    assert adm.ring_fill(s, e1, 50, 70, 0.03, p) == 0.0
    # A hit on the epoch-1 core is not judged as a ring.
    assert adm.ring_fill(s, e1, 50, 51, 0.03, p) == 0.0


def test_ring_needs_min_ring_bands_detected_bands():
    p = adm.Params()
    ring = _band(8.0)
    ring["ring_fill"] = 0.8
    assert adm.classify({"F277W": ring, "F356W": dict(ring)}, p)[0] == "ring"
    assert adm.classify({"F277W": ring, "F356W": _band(8.0)}, p)[0] == "pass"
    undetected = dict(ring, s_hit=2.0)  # a band without the hit cannot add a ring vote
    assert adm.classify({"F277W": ring, "F356W": undetected}, p)[0] == "pass"
