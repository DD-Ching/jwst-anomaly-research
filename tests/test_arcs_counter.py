"""Tests for scripts/arcs_counter.py (an SIS stand-in model; no network)."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

_DIR = Path(__file__).resolve().parents[1] / "scripts"
_spec = importlib.util.spec_from_file_location("arcs_counter", _DIR / "arcs_counter.py")
acn = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(acn)


class _SIS:
    """Singular isothermal sphere at the frame origin; sky = frame (x West, y North, arcsec)."""

    theta_e = 10.0

    def to_frame(self, ra, dec):
        return np.asarray(ra, float), np.asarray(dec, float)

    def source_points(self, x, y, z_s):
        r = np.hypot(x, y)
        return x - self.theta_e * x / r, y - self.theta_e * y / r

    def lens_map(self, x, y, z_s):
        x, y = np.asarray(x, float), np.asarray(y, float)
        r = np.hypot(x, y)
        bx, by = self.source_points(x, y, z_s)
        k = self.theta_e / r**3  # A = I - theta_e / r (I - n n^T)
        return {
            "beta_x": bx,
            "beta_y": by,
            "a11": 1.0 - k * y**2,
            "a12": k * x * y,
            "a21": k * x * y,
            "a22": 1.0 - k * x**2,
        }


def test_nearest_image_finds_the_sis_counter_image():
    # source at beta = (1, 0): images at x = 11 (mu = 11) and x = -9 (mu = -9)
    off, mu = acn.nearest_image(_SIS(), 11.0, 0.0, 2.0, -9.3, 0.2, 3.0, 0.05)
    assert off == pytest.approx(np.hypot(0.3, 0.2), abs=1e-4)
    assert mu == pytest.approx(-9.0, rel=1e-3)


def test_nearest_image_none_in_window():
    off, mu = acn.nearest_image(_SIS(), 11.0, 0.0, 2.0, 0.0, 20.0, 3.0, 0.05)
    assert np.isnan(off) and np.isnan(mu)


def test_catalogued_mask_and_catalog_match():
    arcs = Table({"ra": [10.0, 10.01], "dec": [0.0, 0.0]})
    images = Table({"ra": [10.0 + 1.5 / 3600], "dec": [0.0]})
    assert list(acn.catalogued_mask(arcs, images, 2.0)) == [True, False]
    assert list(acn.catalogued_mask(arcs, images, 1.0)) == [False, False]
    shapes = Table({"label": [1, 2], "ra": [10.0, 10.0], "dec": [0.0, 0.5 / 3600], "mag": [27, 24]})
    lab, mag, sep, ok = acn.catalog_match(shapes, 10.0, 0.0, 1.0, 25.0)
    assert (lab, mag, ok) == (2, 24.0, True)
    assert sep == pytest.approx(0.5, abs=1e-3)
    assert acn.catalog_match(shapes, 10.0, 0.0, 1.0, 23.0)[3] is False
    assert acn.catalog_match(shapes, 10.1, 0.0, 1.0, 30.0)[0] == -1


def test_mag_to_image_units():
    # AB 8.90 is ~1 Jy: 1e6 MJy/sr summed over 1e-12 sr pixels
    assert acn.mag_to_image_units(8.9, 1e-12) == pytest.approx(1e6, rel=1e-3)


def test_known_images_adds_the_cats_list_even_when_gated_off(monkeypatch):
    # abell370-canucs has a gated-off own list; the CATS list of the same cluster is added
    lists = {"abell370-canucs": [10.0], "abell370-cats": [10.01]}
    monkeypatch.setattr(
        acn.lc,
        "MODELS",
        {
            "abell370-canucs": {"files": {"arcs.dat": None}, "image_list_ok": False},
            "abell370-cats": {"files": {"arcs.dat": None}},
            "elgordo-caminha23": {"files": {}},
        },
    )
    monkeypatch.setattr(
        acn.lc, "image_list", lambda n, f, p: (Table({"ra": lists[n], "dec": [0.0]}), {})
    )
    monkeypatch.setattr(acn.lc, "shift_images", lambda n, t: t)
    t = acn.known_images("abell370-canucs", {}, None)
    assert sorted(t["ra"]) == [10.0, 10.01]
    assert list(acn.known_images("abell370-cats", {}, None)["ra"]) == [10.01]
    assert len(acn.known_images("elgordo-caminha23", {}, None)) == 0
