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
    morph = {
        "targets": [
            {"id": "h1", "verdict": "pass"},
            {"id": "h2", "verdict": "rosette"},
            {"id": "h3", "verdict": "pass"},
        ]
    }
    (tmp_path / "morph.json").write_text(json.dumps(morph))
    (tmp_path / "known.csv").write_text(
        "name,ra,dec\n" + "".join(f"{k},{v[0]!r},{v[1]!r}\n" for k, v in sky.items())
    )
    res = ac.known_chain(
        tmp_path / "m.npz",
        tmp_path / "hits.ecsv",
        tmp_path / "morph.json",
        tmp_path / "known.csv",
        "smacs",
    )
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
        )
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


def test_review_list_keeps_off_arc_passes_and_labels_them(tmp_path, monkeypatch):
    import json

    from astropy.table import Table

    w = _wcs()
    _fake_mask(monkeypatch, w)  # one catalogued image at the box centre (pixel 50, 50)
    pos = [
        (50, 50, False, "pass"),
        (10, 10, False, "pass"),
        (52, 52, True, "pass"),
        (45, 45, False, "edge"),
    ]
    sky = [_sky(w, x, y) for x, y, _, _ in pos]
    Table(
        {
            "hit_id": ["h1", "h2", "h3", "h4"],
            "ra": [s[0] for s in sky],
            "dec": [s[1] for s in sky],
            "near_bright_core": [p[2] for p in pos],
            "max_abs_sigma": [6.0, 9.0, 7.0, 8.0],
        }
    ).write(tmp_path / "hits.ecsv")
    (tmp_path / "morph.json").write_text(
        json.dumps({"targets": [{"id": f"h{i + 1}", "verdict": p[3]} for i, p in enumerate(pos)]})
    )
    (tmp_path / "known.csv").write_text(f"name,ra,dec\nK,{sky[0][0]!r},{sky[0][1]!r}\n")
    res = ac.review_list(
        tmp_path / "m.npz",
        tmp_path / "hits.ecsv",
        tmp_path / "morph.json",
        "smacs",
        tmp_path / "known.csv",
    )
    by = {r["id"]: r for r in res["review"]}
    assert set(by) == {"h1", "h2"}  # h3 flagged near a bright core, h4 not a morph pass
    assert by["h1"]["priority"] == "A" and by["h1"]["old_rule"] and by["h1"]["known"] == "K"
    # 1.70" from the image: outside the old 1" rule, priority B, still listed
    assert by["h2"]["priority"] == "B" and not by["h2"]["old_rule"] and not by["h2"]["in_mask"]
    assert (
        res["n_added"] == 1
        and res["n_unflagged"] == 3
        and res["by_priority"] == {"A": 1, "B": 1, "C": 0}
    )


def test_inject_review_counts_every_position(tmp_path, monkeypatch):
    import json

    w = _wcs()
    _fake_mask(monkeypatch, w)
    (tmp_path / "morph.json").write_text(
        json.dumps({"targets": [{"id": "h1", "verdict": "pass"}, {"id": "h2", "verdict": "pass"}]})
    )
    rows = [
        (50, 50, "h1", "False"),
        (10, 10, "h2", "False"),
        (12, 12, "", "False"),
        (55, 55, "h2", "True"),
    ]
    lines = ["name,ra,dec,mag,hit,near_bright_core"]
    for i, (x, y, hit, core) in enumerate(rows):
        ra, dec = _sky(w, x, y)
        lines.append(f"inj{i},{ra!r},{dec!r},27.5,{hit},{core}")
    (tmp_path / "inj.csv").write_text("\n".join(lines) + "\n")
    res = ac.inject_review(
        [(tmp_path / "m.npz", tmp_path / "inj.csv", tmp_path / "morph.json")], "smacs"
    )
    c = res["per_mag"]["27.5"]
    assert (c["all"]["n"], c["all"]["hit"], c["all"]["review"]) == (4, 3, 2)
    assert (c["in_mask"]["n"], c["in_mask"]["review"]) == (
        2,
        1,
    )  # the flagged in-mask star is not reviewed
    assert (c["out_mask"]["n"], c["out_mask"]["review"]) == (2, 1)
    assert c["level_ge3"]["n"] == 2 and c["level_lt3"]["n"] == 2


def _review_inputs(tmp_path, w, pos):
    import json

    from astropy.table import Table

    sky = [_sky(w, x, y) for x, y, _, _, _ in pos]
    Table(
        {
            "hit_id": [f"h{i + 1}" for i in range(len(pos))],
            "ra": [s[0] for s in sky],
            "dec": [s[1] for s in sky],
            "near_bright_core": [p[2] for p in pos],
            "max_abs_sigma": [p[4] for p in pos],
        }
    ).write(tmp_path / "hits.ecsv", overwrite=True)
    targets = [{"id": f"h{i + 1}", "verdict": p[3]} for i, p in enumerate(pos) if p[3]]
    (tmp_path / "morph.json").write_text(json.dumps({"targets": targets}))


def test_review_list_priority_c_sort_and_labels(tmp_path, monkeypatch):
    w = _wcs(n=300)
    mask = np.zeros((300, 300), bool)
    mask[140:160, 140:160] = True
    snr = np.where(mask, 5.0, 0.0)
    monkeypatch.setattr(
        ac, "load_mask", lambda maps, field, with_snr=False: (mask, w, "F444W", 0.03, snr)
    )
    monkeypatch.setattr(ac, "image_positions", lambda field: (np.array([10.0]), np.array([-20.0])))
    # centre pixel 150: A in the mask; C at 120 px = 3.6" away; B at 50 px = 1.5" away
    _review_inputs(
        tmp_path,
        w,
        [
            (150, 150, False, "pass", 6.0),
            (270, 150, False, "pass", 9.0),
            (150, 200, False, "pass", 7.0),
        ],
    )
    res = ac.review_list(
        tmp_path / "m.npz", tmp_path / "hits.ecsv", tmp_path / "morph.json", "smacs"
    )
    assert [r["priority"] for r in res["review"]] == ["A", "B", "C"]
    a, b, c = res["review"]
    assert a["mask_distance_arcsec"] == 0 and a["arc_level_sigma"] == 5.0
    assert b["mask_distance_arcsec"] > 0 and b["arc_level_sigma"] == 0.0 and not b["old_rule"]
    assert c["nearest_image_arcsec"] > 2.0 and c["id"] == "h2"


def test_review_list_and_inject_review_refuse_missing_verdicts(tmp_path, monkeypatch):
    import pytest

    w = _wcs()
    _fake_mask(monkeypatch, w)
    _review_inputs(tmp_path, w, [(50, 50, False, "pass", 6.0), (45, 45, False, None, 6.0)])
    with pytest.raises(ValueError):
        ac.review_list(tmp_path / "m.npz", tmp_path / "hits.ecsv", tmp_path / "morph.json", "smacs")
    ra, dec = _sky(w, 45, 45)
    (tmp_path / "inj.csv").write_text(
        f"name,ra,dec,mag,hit,near_bright_core\ni0,{ra!r},{dec!r},27.5,h2,False\n"
    )
    with pytest.raises(ValueError):
        ac.inject_review(
            [(tmp_path / "m.npz", tmp_path / "inj.csv", tmp_path / "morph.json")], "smacs"
        )


def test_inject_review_outcomes_and_magnitudes(tmp_path, monkeypatch):
    import json

    w = _wcs()
    _fake_mask(monkeypatch, w)
    verdicts = {"h1": "pass", "h2": "rosette", "h3": "pass"}
    (tmp_path / "morph.json").write_text(
        json.dumps({"targets": [{"id": k, "verdict": v} for k, v in verdicts.items()]})
    )
    rows = [
        (50, 50, "27.5", "h1", "False"),
        (10, 10, "27.5", "h2", "False"),
        (12, 12, "27.5", "h3", "True"),
        (14, 14, "28.0", "", "False"),
    ]
    lines = ["name,ra,dec,mag,hit,near_bright_core"]
    for i, (x, y, mag, hit, core) in enumerate(rows):
        ra, dec = _sky(w, x, y)
        lines.append(f"inj{i},{ra!r},{dec!r},{mag},{hit},{core}")
    (tmp_path / "inj.csv").write_text("\n".join(lines) + "\n")
    out = tmp_path / "out.json"
    ac.main(
        [
            "inject-review",
            "--field",
            "smacs",
            "--maps",
            "m.npz",
            "--injected",
            str(tmp_path / "inj.csv"),
            "--morph",
            str(tmp_path / "morph.json"),
            "--out",
            str(out),
        ]
    )
    per = json.loads(out.read_text())["per_mag"]
    a = per["27.5"]["all"]
    assert (a["n"], a["hit"], a["review"], a["eff_hit"], a["eff_review"]) == (3, 3, 1, 1.0, 0.333)
    assert a["hit_outcomes"] == {"pass": 1, "rosette": 1, "near_bright_core": 1}
    assert per["27.5"]["level_ge3"]["review"] == 1 and per["27.5"]["level_lt3"]["review"] == 0
    assert per["28.0"]["all"] == {
        "n": 1,
        "hit": 0,
        "review": 0,
        "hit_outcomes": {},
        "eff_hit": 0.0,
        "eff_review": 0.0,
    }


def test_review_cli_without_known(tmp_path, monkeypatch):
    import json

    w = _wcs()
    _fake_mask(monkeypatch, w)
    _review_inputs(tmp_path, w, [(50, 50, False, "pass", 6.0)])
    out = tmp_path / "r.json"
    ac.main(
        [
            "review",
            "--field",
            "smacs",
            "--maps",
            "m.npz",
            "--hits",
            str(tmp_path / "hits.ecsv"),
            "--morph",
            str(tmp_path / "morph.json"),
            "--out",
            str(out),
        ]
    )
    res = json.loads(out.read_text())
    assert res["n_review"] == 1 and res["review"][0]["known"] is None
