"""Unit tests for scripts/arc_chain.py (offline, synthetic)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np


def _load():
    path = Path(__file__).resolve().parents[1] / "scripts" / "arc_chain.py"
    spec = importlib.util.spec_from_file_location("arc_chain", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["arc_chain"] = mod
    spec.loader.exec_module(mod)
    return mod


ac = _load()


def _image(n=300):
    rng = np.random.default_rng(0)
    img = rng.normal(0, 1.0, (n, n)) + 5.0  # sky pedestal
    yy, xx = np.indices(img.shape)
    # a thin curved arc (seeded) and a compact galaxy far away (not seeded)
    r = np.hypot(xx - 150, yy + 200)
    img += 20 * np.exp(-((r - 300) ** 2) / (2 * 2.0**2)) * (np.abs(xx - 150) < 80)
    img += 50 * np.exp(-((xx - 40) ** 2 + (yy - 40) ** 2) / (2 * 3.0**2))
    return img


def test_arc_mask_keeps_seeded_arc_and_drops_other_sources():
    img = _image()
    scale = 0.03
    seeds = [(150.0, 100.0)]  # on the arc (r = 300 at x = 150 -> y = 100)
    m = ac.arc_mask(img, scale, seeds)
    assert m[100, 150] and m[100, 200]  # along the arc
    assert not m[40, 40]  # the unseeded compact galaxy
    assert m.mean() < 0.25


def test_arc_mask_empty_without_seed_on_a_source():
    img = _image()
    m = ac.arc_mask(img, 0.03, [(20.0, 280.0)])  # blank sky
    assert not m.any()


def _wcs(n=100, scale=0.03):
    from astropy.wcs import WCS

    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [10.0, -20.0]
    w.wcs.crpix = [n / 2 + 1, n / 2 + 1]
    w.wcs.cdelt = [-scale / 3600, scale / 3600]
    return w


def _sky(w, x, y):
    ra, dec = w.all_pix2world(x, y, 0)
    return float(ra), float(dec)


def test_in_mask_inside_outside_and_off_image():
    w = _wcs()
    mask = np.zeros((100, 100), bool)
    mask[40:60, 40:60] = True
    assert ac.in_mask(mask, w, *_sky(w, 50, 50))
    assert not ac.in_mask(mask, w, *_sky(w, 10, 10))
    assert not ac.in_mask(mask, w, *_sky(w, -5, 50))  # negative pixel: no wrap-around
    assert not ac.in_mask(mask, w, *_sky(w, 50, 150))


def test_filtered_snr_rejects_blank_box():
    import pytest

    with pytest.raises(ValueError):
        ac.filtered_snr(np.zeros((50, 50)), 0.03)


def _fake_mask(monkeypatch, w):
    mask = np.zeros((100, 100), bool)
    mask[40:60, 40:60] = True
    snr = np.where(mask, 5.0, 0.0)
    monkeypatch.setattr(
        ac,
        "load_mask",
        lambda maps, field, with_snr=False: (
            (mask, w, "F444W", 0.03, snr) if with_snr else (mask, w, "F444W", 0.03)
        ),
    )
    monkeypatch.setattr(ac, "image_positions", lambda field: (np.array([10.0]), np.array([-20.0])))


def test_known_chain_stages(tmp_path, monkeypatch):
    import json

    from astropy.table import Table

    w = _wcs()
    _fake_mask(monkeypatch, w)
    pos = {"in_pass": (50, 50), "in_veto": (45, 45), "out_pass": (10, 10), "missed": (55, 55)}
    sky = {k: _sky(w, *v) for k, v in pos.items()}
    hits = Table(
        {
            "hit_id": ["h1", "h2", "h3"],
            "ra": [sky["in_pass"][0], sky["in_veto"][0], sky["out_pass"][0]],
            "dec": [sky["in_pass"][1], sky["in_veto"][1], sky["out_pass"][1]],
            "near_bright_core": [False, False, False],
        }
    )
    hits.write(tmp_path / "hits.ecsv")
    morph = {"targets": [{"id": "h1", "verdict": "pass"}, {"id": "h2", "verdict": "rosette"},
                         {"id": "h3", "verdict": "pass"}]}  # fmt: skip
    (tmp_path / "morph.json").write_text(json.dumps(morph))
    (tmp_path / "known.csv").write_text(
        "name,ra,dec\n" + "".join(f"{k},{v[0]!r},{v[1]!r}\n" for k, v in sky.items())
    )
    res = ac.known_chain(
        tmp_path / "m.npz", tmp_path / "hits.ecsv", tmp_path / "morph.json",
        tmp_path / "known.csv", "smacs",
    )  # fmt: skip
    by = {r["name"]: r for r in res["known"]}
    assert by["in_pass"]["full_chain"] and by["in_pass"]["mask_distance_arcsec"] == 0
    assert by["in_veto"]["morph"] == "rosette" and not by["in_veto"]["full_chain"]
    assert by["out_pass"]["morph"] == "pass" and not by["out_pass"]["in_mask"]
    assert by["out_pass"]["mask_distance_arcsec"] > 0 and not by["out_pass"]["full_chain"]
    assert by["missed"]["hit"] is None and by["missed"]["in_mask"]
    assert res["mask_area_ok"] and {h["id"] for h in res["hits_in_mask"]} == {"h1", "h2"}


def test_inject_chain_counts_and_pools(tmp_path, monkeypatch):
    import json

    w = _wcs()
    _fake_mask(monkeypatch, w)
    (tmp_path / "morph.json").write_text(
        json.dumps(
            {
                "targets": [
                    {"id": "h1", "verdict": "pass"},
                    {"id": "h2", "verdict": "edge"},
                    {"id": "h3", "verdict": "pass"},
                ]
            }
        )  # fmt: skip
    )
    rows = [  # (x, y, mag, hit, near_bright_core)
        (50, 50, "27.5", "h1", "False"),
        (45, 45, "27.5", "h2", "False"),
        (55, 55, "27.5", "", "False"),
        (52, 52, "27.5", "h3", "True"),
        (10, 10, "27.5", "h1", "False"),  # outside the mask: not counted
    ]
    lines = ["name,ra,dec,mag,hit,near_bright_core"]
    for i, (x, y, mag, hit, core) in enumerate(rows):
        ra, dec = _sky(w, x, y)
        lines.append(f"inj{i},{ra!r},{dec!r},{mag},{hit},{core}")
    (tmp_path / "inj.csv").write_text("\n".join(lines) + "\n")
    run = (tmp_path / "m.npz", tmp_path / "inj.csv", tmp_path / "morph.json")
    c = ac.inject_chain([run, run], "smacs")["per_mag"]["27.5"]
    assert (c["n_in_mask"], c["hit"], c["chain"], c["chain_unflagged"]) == (8, 6, 4, 2)
    assert c["eff_chain"] == 0.5 and c["eff_chain_unflagged"] == 0.25
