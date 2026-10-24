import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import spt0546_inputs as sp  # noqa: E402

MUL = (
    "#REFERENCE 0\n#fk5\n"
    "1.1 86.6548684 -53.7582990 0.2 0.2 0. 0. 25.\n"
    "1.2 86.6541568 -53.7583831 0.2 0.2 0. 0. 25.\n"
    "c2.3 86.6594888 -53.7638632 0.2 0.2 0. 0. 25.\n"
    "10.2 86.6662699 -53.7606063 0.2 0.2 0. 0. 25.\t\t\n"
)


def _header(n=101):
    h = fits.Header()
    h.update(
        CTYPE1="RA---TAN",
        CTYPE2="DEC--TAN",
        CRVAL1=86.65609133,
        CRVAL2=-53.76055403,
        CRPIX1=(n + 1) / 2,
        CRPIX2=(n + 1) / 2,
        CDELT1=-sp.PIXEL_ARCSEC / 3600,
        CDELT2=sp.PIXEL_ARCSEC / 3600,
    )
    return h


def test_parse_mul_ids_systems_candidates():
    t = sp.parse_mul(MUL)
    assert list(t["image_id"]) == ["1.1", "1.2", "c2.3", "10.2"]
    assert list(t["system"]) == ["1", "1", "2", "10"]
    assert list(t["candidate"]) == [False, False, True, False]
    with pytest.raises(ValueError):
        sp.parse_mul("#only a comment\n")


def test_specz_systems_rule():
    images = sp.parse_mul(MUL)
    cat = Table(
        rows=[
            # 0.3" from 1.1, behind the lens, conf 4: counts
            (1, 86.6548684 + 0.3 / 3600 / np.cos(np.radians(53.758)), -53.7582990, 3.1, 4),
            # on 10.2 but a cluster member: does not count
            (2, 86.6662699, -53.7606063, 1.12, 4),
            # on c2.3 (a candidate): does not count
            (3, 86.6594888, -53.7638632, 2.5, 3),
        ],
        names=("id", "ra", "dec", "z", "conf"),
    )
    r = sp.specz_systems(images, cat)
    assert r["n_systems"] == 1 and list(r["systems"]) == ["1"]
    assert r["systems"]["1"][0]["sep"] == pytest.approx(0.3, abs=0.01)
    cat["conf"][0] = 0  # unknown confidence never counts
    assert sp.specz_systems(images, cat)["n_systems"] == 0


def test_unit_check_needs_exactly_one_reading():
    one = sp.unit_check({"A_arcsec": {"1": 0.1, "2": 0.6, "3": 0.3}, "P_pixels": {"1": 9.0}})
    assert one["passes"] and one["reading"] == "A_arcsec"
    both = sp.unit_check({"A_arcsec": {"1": 0.1}, "P_pixels": {"1": 0.5}})
    assert not both["passes"] and both["reading"] is None
    none = sp.unit_check({"A_arcsec": {"1": 0.51}, "P_pixels": {"1": 0.1, "2": float("inf")}})
    assert not none["passes"]  # an unfitted system (inf) never accepts a reading


def test_build_model_normalises_to_dlsds1_and_covers():
    n = 101
    ax = np.full((n, n), 2.0)
    ay = np.zeros((n, n))
    m = sp.build_model(ax, ay, _header(n), sp.READINGS["P_pixels"], "x", "t")
    x, y = m.deflection_xy(np.array([0.0]), np.array([0.0]))
    assert float(x[0]) == pytest.approx(2.0 * 0.04 / sp.dls_ds(sp.Z_NORM), rel=1e-5)
    images = sp.parse_mul(MUL)
    c = sp.coverage(m, images)  # a 4" map: only the images near the centre are inside
    assert not c["passes"] and c["inside"] + len(c["outside"]) == 4
    assert sp.dls_ds(3.5) == pytest.approx(0.48931, abs=1e-4)


def test_unit_minima_point_lens_recovers_z():
    # a point lens of known Einstein radius: two images on a line through the centre map to one
    # source only at the true z, so the minimum rms sits at (near) that grid z
    n = 401
    g = (np.arange(n) - (n - 1) / 2) * sp.PIXEL_ARCSEC
    xx, yy = np.meshgrid(g, g)
    r2 = xx**2 + yy**2 + 1e-6
    te2 = 4.0  # theta_E^2 at D_LS/D_S = 1 (arcsec^2)
    k = sp.dls_ds(sp.Z_NORM)
    # map values in arcsec at z_s 3.5 (alpha = theta_E^2 theta / |theta|^2 along +i West, +j North)
    m = sp.build_model(te2 * xx / r2 * k, te2 * yy / r2 * k, _header(n), 1.0, "x", "pl")
    zs = sp.z_grid()
    z_true = float(zs[30])
    te = np.sqrt(te2 * sp.dls_ds(z_true))
    b = 0.8  # source offset (arcsec, North)
    t1 = 0.5 * (b + np.sqrt(b * b + 4 * te * te))
    t2 = 0.5 * (b - np.sqrt(b * b + 4 * te * te))
    dec = m.dec0 + np.array([t1, t2]) / 3600.0
    images = Table(
        {"image_id": ["1.1", "1.2"], "system": ["1", "1"], "ra": [m.ra0] * 2, "dec": dec}
    )
    images["a"] = np.ones(2)
    out = sp.unit_minima(m, images, zs)
    assert out["1"]["min_rms_arcsec"] < 0.05
    assert out["1"]["z_at_min"] == pytest.approx(z_true, rel=0.05)
    assert not out["1"]["z_at_min_on_edge"]


def test_dja_check():
    html = "<a>abell2744clu-grizli-v7.4</a> <a>spt0546-grizli-v7.0</a>"
    assert sp.dja_check(html)["matches"] == ["spt0546"]
    assert not sp.dja_check("<a>macs0416-grizli-v7.0</a>")["passes"]
    with pytest.raises(ValueError):
        sp.dja_check("nothing")


def test_unit_check_none_minimum_rejects_reading():
    r = sp.unit_check({"A_arcsec": {"1": 0.1, "2": None}, "P_pixels": {"1": 0.2, "2": 0.3}})
    assert r["passes"] and r["reading"] == "P_pixels"


def test_parse_specz_comments_and_short_rows():
    good = "# Id RA Dec\n  # indented comment\n1 86.6 -53.7 0.5 0.5 0 1.2 22.0 4\n"
    assert len(sp.parse_specz(good)) == 1
    with pytest.raises(ValueError):
        sp.parse_specz("1 86.6 -53.7 0.5\n")
