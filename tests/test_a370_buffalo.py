"""Tests of the Abell 370 BUFFALO model A gate-K script: list, selection, gate L (PLAN.md)."""

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import a370_buffalo as m  # noqa: E402
import lens_consistency as lc  # noqa: E402

# Lines in the CDS byte layout (copied from the format, not all from the real table)
ROWS = [
    " 1.1 39.967047 -1.5769172 0.8041             gold   -",
    " 3.3 39.978925 -1.5674624 1.9553             bronze -",
    "10.1 39.970000 -1.5800000 2.7512             gold   1",
    "16.2 39.971000 -1.5810000 3.7743             gold   1,2",
    "21.1 39.972000 -1.5820000 1.2567             gold   2",
    "30.1 39.973000 -1.5830000        1.500 0.100 gold   -",
]


def test_constants_name_the_planned_field_and_model():
    assert (m.FIELD, m.MODEL, m.GRID_STEP_K, m.L_MAX_MEDIAN_REL) == (
        "abell370",
        "abell370-buffalo-a-maps",
        0.25,
        0.05,
    )
    assert m.OUT == Path("results/lens_consistency/abell370_buffalo")
    spec = lc.MODELS[m.MODEL]
    assert spec["kind"] == "maps" and spec["z_lens"] == 0.375
    assert {"alpha_x", "alpha_y", "kappa_map"} <= set(spec["files"])


def test_parse_and_select_follow_the_plan(tmp_path):
    p = tmp_path / "table1.dat"
    p.write_text("\n".join(ROWS) + "\n")
    t = m.parse_table1(p)
    assert list(t["image_id"]) == ["1.1", "3.3", "10.1", "16.2", "21.1", "30.1"]
    assert t["ra"][0] == 39.967047 and t["dec"][0] == -1.5769172
    assert np.isnan(t["zspec"][5]) and t["zspec"][0] == 0.8041
    assert list(t["class"]) == ["gold", "bronze", "gold", "gold", "gold", "gold"]
    # gold with zspec, Note without 1: 1.1 and 21.1 (Note 2 only) are kept
    assert list(t["image_id"][m.select(t)]) == ["1.1", "21.1"]


def test_gate_l_threshold():
    assert m.gate_l({"median_rel_diff": 0.05})["passes"]
    assert not m.gate_l({"median_rel_diff": 0.051})["passes"]
    assert not m.gate_l({"median_rel_diff": float("nan")})["passes"]


@pytest.mark.network
def test_written_images_are_the_102_selected_gold_spectroscopic_images(tmp_path):
    tab = Table.read(m.write_images(tmp_path), format="ascii.ecsv")
    assert len(tab) == 102 and bool(np.all(tab["z_spec"]))
    assert np.all(np.isfinite(np.asarray(tab["z"], float)))


def test_k1_overlap_splits_by_distance_to_gold(tmp_path):
    lag = Table({"image_id": ["1.1", "2.1"], "ra": [10.0, 10.01], "dec": [0.0, 0.0]})
    gold = Table({"image_id": ["1.1"], "ra": [10.0 + 0.3 / 3600.0], "dec": [0.0]})
    k1 = Table({"image_id": ["1.1", "2.1"], "dtheta_arcsec": [0.4, np.nan]})
    lag.write(tmp_path / m.IMAGES, format="ascii.ecsv")
    gold.write(tmp_path / m.IMAGES_K2, format="ascii.ecsv")
    k1.write(tmp_path / f"imageplane_{m.FIELD}_{m.MODEL}.ecsv", format="ascii.ecsv")
    out = m.k1_overlap(tmp_path)
    assert out["shared"]["n"] == 1 and out["shared"]["within_3arcsec"] == 1
    assert np.isclose(out["shared"]["rms_arcsec"], 0.4)
    assert out["not_shared"]["image_ids"] == ["2.1"] and out["not_shared"]["n_finite"] == 0
    assert out["nearest_gold_arcsec"]["1.1"] == pytest.approx(0.3, abs=1e-3)


def test_k2_decision():
    assert m.k2_decision({"passing_models": [m.MODEL]}).startswith("passes")
    assert m.k2_decision({"passing_models": []}).startswith("fails")
