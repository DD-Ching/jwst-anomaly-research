"""Tests for the pure helpers of scripts/vet_evidence.py."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS

_PATH = Path(__file__).resolve().parents[1] / "scripts" / "vet_evidence.py"
_spec = importlib.util.spec_from_file_location("vet_evidence", _PATH)
vet = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(vet)


@pytest.mark.parametrize(
    "a, b, expected",
    [
        (10.0, 10.0, 0.0),
        (0.0, 90.0, 90.0),
        (170.0, 10.0, 20.0),
        (45.0, 225.0, 0.0),
        (5.0, 95.0, 90.0),
    ],
)
def test_axis_offset_deg(a, b, expected):
    assert vet.axis_offset_deg(a, b) == pytest.approx(expected)


def _elongated(path: Path, pa_deg: float) -> None:
    """Elongated Gaussian whose major axis points pa_deg east of north (north up, east left)."""
    n = 81
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [150.0, 2.0]
    w.wcs.crpix = [(n + 1) / 2, (n + 1) / 2]
    w.wcs.cdelt = [-0.05 / 3600, 0.05 / 3600]
    yy, xx = np.mgrid[0:n, 0:n] - (n - 1) / 2
    theta = np.deg2rad(pa_deg)
    # unit vector of the major axis in pixel coords: east = -x, north = +y
    ux, uy = -np.sin(theta), np.cos(theta)
    along = xx * ux + yy * uy
    across = -xx * uy + yy * ux
    img = np.exp(-0.5 * ((along / 12.0) ** 2 + (across / 4.0) ** 2))
    img += np.random.default_rng(0).normal(0, 1e-3, img.shape)
    fits.PrimaryHDU(img.astype(np.float32), header=w.to_header()).writeto(path)


@pytest.mark.parametrize("pa", [0.0, 30.0, 75.0, 120.0])
def test_moment_orientation_recovers_position_angle(tmp_path, pa):
    path = tmp_path / f"e{int(pa)}.fits"
    _elongated(path, pa)
    m = vet.moment_orientation(path)
    assert m is not None
    assert vet.axis_offset_deg(m["pa_deg"], pa) < 2.0
    assert m["axis_ratio"] == pytest.approx(4.0 / 12.0, rel=0.25)


def test_load_lens_images_parses_lenstool_arcs(tmp_path):
    f = tmp_path / "arcs.dat"
    f.write_text(
        "# comment\n1.1 110.84 -73.45 0.1 0.1 0 1.449 0.0\n\n2.3 110.83 -73.46 0.5 0.5 0 0 0\n"
    )
    t = vet.load_lens_images(str(f))
    assert list(t["image_id"]) == ["1.1", "2.3"]
    assert t["ra"][1] == pytest.approx(110.83)


def test_moments_growth_reports_the_size_measured(monkeypatch, tmp_path):
    from astropy.table import Table

    sizes = []

    def fake_cutouts(uri, targets, size_arcsec, out_dir):
        sizes.append(size_arcsec)
        return Table({"path": ["x.fits"]})

    monkeypatch.setattr(vet.cutouts, "make_cutouts", fake_cutouts)
    monkeypatch.setattr(
        vet, "moment_orientation", lambda path: {"pa_deg": 0.0, "touches_border": True}
    )
    moments, err = vet._moments_with_growth("uri", Table(), 4.0, tmp_path)
    assert err is None
    assert sizes == [4.0, 8.0, 16.0, 32.0]  # MAX_GROWTH = 3 doublings
    assert moments["cutout_arcsec"] == 32.0 and moments["touches_border"]


def _zout(path, **extra):
    from astropy.table import Table

    t = Table(
        {
            "id": [1, 2],
            "ra": [10.0, 10.01],
            "dec": [-5.0, -5.0],
            "z_phot": [2.5, 4.8],
            "z025": [2.3, 4.6],
            "z160": [2.4, 4.7],
            "z500": [2.5, 4.8],
            "z840": [2.6, 4.9],
            "z975": [2.8, 5.4],
            "z_phot_chi2": [10.0, 3.5],
            "z_spec": [-1.0, 4.81],
            **extra,
        }
    )
    t.meta.update(VERSION="0.8.3", TEMPLATES_FILE="templates/x.param")
    t.write(path, overwrite=True)
    return path


def test_photoz_evidence_nearest_within_radius(tmp_path):
    from astropy.coordinates import SkyCoord

    table = vet.load_photoz(_zout(tmp_path / "x.zout.fits", nusefilt=[6, 16]))
    assert table.meta["provenance"] == "model_prediction" and len(table.meta["sha256"]) == 64
    near = SkyCoord(10.0, -5.0 + 0.05 / 3600, unit="deg")
    ev = vet.photoz_evidence(table, near, 0.2)
    assert ev["match"]["id"] == 1 and ev["match"]["z_spec"] is None  # -1 = no z_spec
    assert ev["match"]["sep_arcsec"] == pytest.approx(0.05, abs=1e-3)
    assert ev["match"]["nusefilt"] == 6 and ev["eazy_version"] == "0.8.3"
    second = vet.photoz_evidence(table, SkyCoord(10.01, -5.0, unit="deg"), 0.2)
    assert second["match"]["z_spec"] == pytest.approx(4.81)
    far = SkyCoord(10.0, -5.0 + 1.0 / 3600, unit="deg")
    assert vet.photoz_evidence(table, far, 0.2)["match"] is None


def test_load_photoz_rejects_other_tables(tmp_path):
    from astropy.table import Table

    Table({"ra": [1.0], "dec": [2.0]}).write(tmp_path / "x.fits")
    with pytest.raises(ValueError, match="zout"):
        vet.load_photoz(tmp_path / "x.fits")
