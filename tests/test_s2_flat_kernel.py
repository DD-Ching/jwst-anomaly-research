"""Tests for scripts/s2_flat_kernel.py (synthetic tables, offline)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("s2_flat_kernel", _DIR / "s2_flat_kernel.py")
s2 = importlib.util.module_from_spec(_spec)
sys.modules["s2_flat_kernel"] = s2  # dataclasses look their module up
_spec.loader.exec_module(s2)

P = s2.Params()


def test_lens_weight_vanishes_at_both_ends_and_behind():
    c = FlatLambdaCDM(H0=70, Om0=0.3)
    w = s2.lens_weight(np.array([1e-4, 0.3, 0.8, 1.2]), 0.8, c)
    assert w[0] < 1e-3 and w[1] > 0 and w[2] == 0 and w[3] == 0


def test_box_terms_wraps_ra():
    assert len(s2.box_terms(0.01, 0.0, 120)) == 2
    assert len(s2.box_terms(180.0, 0.0, 120)) == 1


def test_wls_recovers_coefficients():
    rng = np.random.default_rng(1)
    n = 2000
    z = rng.uniform(0.1, 1, n)
    xl, xf = rng.normal(0, 0.3, n), rng.normal(0, 0.3, n)
    err = np.full(n, 0.1)
    y = 0.2 - 0.1 * z - 0.05 * xl + 0.03 * xf + rng.normal(0, 0.1, n)
    beta, sig, rchi2 = s2.wls(y, err, z, xl, xf)
    assert abs(beta[2] + 0.05) < 4 * sig[2] and abs(beta[3] - 0.03) < 4 * sig[3]
    assert 0.8 < rchi2 < 1.2


def test_scramble_keeps_redshift_bins():
    rng = np.random.default_rng(2)
    z = np.repeat([0.12, 0.33, 0.71], 50)
    xl = np.arange(150.0)
    sl, sf = s2.scramble(z, xl, -xl, rng)
    for zb in (0.12, 0.33, 0.71):
        m = z == zb
        assert set(sl[m]) == set(xl[m])
    assert np.array_equal(sf, -sl)  # pairs move together


def test_columns_flag_overdense_line_of_sight():
    rng = np.random.default_rng(3)
    sne = Table({"ra": np.arange(20) * 1.0 + 10, "dec": np.zeros(20), "z": np.full(20, 0.8)})
    ras, decs, zs = [], [], []
    for i, s in enumerate(sne):
        n = 60 if i == 0 else 30  # SN 0 has twice the foreground counts
        r = 60 * np.sqrt(rng.uniform(0, 1, n)) / 3600
        t = rng.uniform(0, 2 * np.pi, n)
        ras.append(s["ra"] + r * np.cos(t))
        decs.append(r * np.sin(t))
        zs.append(rng.uniform(0.05, 0.7, n))
    n = sum(len(a) for a in ras)
    gal = Table(
        {
            "ra": np.concatenate(ras),
            "dec": np.concatenate(decs),
            "release": np.full(n, 9011),
            "z_spec": np.full(n, -99.0),
            "z_phot_median": np.concatenate(zs),
        }
    )
    c = s2.columns(sne, gal, P)
    assert c["xf"][0] > 0.5 and c["xl"][0] > 0.4
    assert np.all(c["xf"][1:] < 0.5)


def test_alpha2_column_is_unbiased_for_poisson():
    rng = np.random.default_rng(4)
    sne = Table({"ra": np.arange(200) * 0.5 + 10, "dec": np.zeros(200), "z": np.full(200, 0.8)})
    ras, decs, zs = [], [], []
    for s in sne:
        n = rng.poisson(40)
        r = 60 * np.sqrt(rng.uniform(0, 1, n)) / 3600
        t = rng.uniform(0, 2 * np.pi, n)
        ras.append(s["ra"] + r * np.cos(t))
        decs.append(r * np.sin(t))
        zs.append(rng.uniform(0.05, 0.7, n))
    n = sum(len(a) for a in ras)
    gal = Table(
        {
            "ra": np.concatenate(ras),
            "dec": np.concatenate(decs),
            "release": np.full(n, 9011),
            "z_spec": np.full(n, -99.0),
            "z_phot_median": np.concatenate(zs),
        }
    )
    c = s2.columns(sne, gal, s2.Params(alpha=2))
    assert abs(np.nanmean(c["xf"])) < 0.05  # unclustered field: <(1 + delta)^2> = 1
    # a pure Poisson field has no signal in the column: dilution factor near 0
    sel = np.isfinite(c["xf"])
    lam_l, lam_f = s2.attenuation(c, sne, s2.Params(alpha=2), rng, sel, n_sim=10)
    assert abs(lam_l) < 0.3 and abs(lam_f) < 0.3


def test_unknown_release_is_ignored():
    sne = Table({"ra": [10.0], "dec": [0.0], "z": [0.8]})
    gal = Table(
        {
            "ra": [10.0, 10.001],
            "dec": [0.0, 0.0],
            "release": [9999, 9999],
            "z_spec": [-99.0, -99.0],
            "z_phot_median": [0.3, 0.4],
        }
    )
    c = s2.columns(sne, gal, P)
    assert c["region"][0] == -1 and np.isnan(c["xf"][0])


def test_z_perm_stays_in_bin():
    z = np.repeat([0.12, 0.33], 30)
    perm = s2.z_perm(z, np.random.default_rng(5))
    assert np.array_equal(z[perm], z) and sorted(perm) == list(range(60))


def test_load_des_joins_hd_meta_and_positions(tmp_path, monkeypatch):
    hd = (
        "# comment\nVARNAMES: CID IDSURVEY zHD zHEL MU MUERR PROBIA_BEAMS\n"
        "SN: 101 10 0.50 0.50 42.2 0.15 0.99\n"
        "SN: 102 10 0.60 0.60 42.6 0.20 0.10\n"  # likely core collapse: cut
        "SN: 103 10 0.05 0.05 36.7 0.12 1.00\n"  # below z_min: cut
        "SN: 101 150 0.03 0.03 35.4 0.12 1.00\n"  # same CID, other survey: ignored
        "SN: 104 10 0.70 0.70 43.0 0.18 0.90\n"  # no position: ignored
    )
    meta = (
        "VARNAMES: CID IDSURVEY c\n"
        "SN: 101 10 -0.05\nSN: 102 10 0.1\nSN: 103 10 0.0\nSN: 101 150 0.2\nSN: 104 10 0.0\n"
    )
    (tmp_path / "DES-Dovekie_HD.csv").write_text(hd)
    (tmp_path / "DES-Dovekie_Metadata.csv").write_text(meta)
    head = Table(
        {"SNID": ["101", "102", "103"], "RA": [35.1, 35.2, 9.0], "DEC": [-5.0, -5.1, -43.0]}
    )
    head.write(tmp_path / "DES-SN5YR_DES_HEAD.FITS.gz")
    monkeypatch.setattr(s2, "cache", lambda: tmp_path)
    t = s2.load_des(s2.Params(sample="des"))
    assert list(t["cid"]) == ["101"]
    assert t["ra"][0] == 35.1 and t["m"][0] == 42.2 and t["c"][0] == -0.05 and t["err"][0] == 0.15
