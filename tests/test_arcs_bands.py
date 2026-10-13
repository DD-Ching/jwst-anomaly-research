"""Tests for scripts/arcs_bands.py (synthetic tables; no network)."""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("arcs_bands", _DIR / "arcs_bands.py")
ab = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ab)


def test_s3_url_and_manifests():
    fn = "jw01208-o002_t001_nircam_clear-f444w_cat.ecsv"
    assert ab.s3_url(fn) == (
        "https://stpubdata.s3.amazonaws.com/jwst/public/jw01208/L3/t/o002/" + fn
    )
    for manifest, filename in ab.FIELDS.values():
        assert len(ab.manifest_sha(manifest, filename)) == 64


def test_match_second_band_radius_and_axis_offset():
    stored = Table({"ra": [10.0, 10.001], "dec": [0.0, 0.0], "pa_obs": [170.0, 20.0]})
    shapes = Table(
        {
            "ra": [10.0 + 0.1 / 3600, 10.001 + 2.0 / 3600],
            "dec": [0.0, 0.0],
            "pa_obs": [5.0, 20.0],
            "ellipticity": [0.4, 0.6],
            "snr": [30.0, 40.0],
        }
    )
    t = ab.match_second_band(stored, shapes, match_arcsec=0.3)
    assert np.isclose(t["dpa_bands"][0], 15.0)  # axes wrap at 180 deg
    assert np.isclose(t["e_444"][0], 0.4)
    assert np.isnan(t["pa_444"][1])  # 2" away: no match
    assert np.isclose(t["sep_nearest"][1], 2.0, atol=0.01)


def test_match_second_band_flags_shared_match():
    stored = Table({"ra": [10.0, 10.0 + 0.2 / 3600], "dec": [0.0, 0.0], "pa_obs": [0.0, 0.0]})
    shapes = Table(
        {
            "ra": [10.0 + 0.1 / 3600],
            "dec": [0.0],
            "pa_obs": [0.0],
            "ellipticity": [0.5],
            "snr": [50.0],
        }
    )
    t = ab.match_second_band(stored, shapes, match_arcsec=0.3)
    assert list(t["shared_444"]) == [True, True]
    assert np.all(np.isnan(t["pa_444"]))


def test_summarize_counts_control_and_anti():
    t = Table(
        {
            "label": [1, 2, 3, 4],
            "strong_shear": [True, True, True, True],
            "population": ["background"] * 4,
            "orientation_class": ["aligned", "aligned", "anti", "anti"],
            "class_444": ["aligned", "mixed", "anti", "aligned"],
            "pa_444": [10.0, 20.0, 30.0, np.nan],
            "dpa_bands": [1.0, 40.0, 2.0, np.nan],
            "e_444": [0.5, 0.5, 0.5, np.nan],
            "snr_444": [50.0, 50.0, 50.0, np.nan],
            "g_min": [0.5, 0.5, 0.5, 0.5],
            "shared_444": [False, False, False, False],
            "sep_nearest": [0.1, 0.1, 0.1, 1.5],
        }
    )
    args = argparse.Namespace(
        min_e444=0.3, min_snr444=10.0, null_e=0.5, aligned_deg=30.0, anti_deg=60.0, sigma_eps=0.26
    )
    s = ab.summarize(t, args)
    assert s["usable_444"] == 3
    assert (s["control_aligned_200"], s["control_stay_aligned_444"]) == (2, 1)
    assert s["anti_444"] == [3]
    assert [a["class_444"] for a in s["anti_200"]] == ["anti", "no_match"]
