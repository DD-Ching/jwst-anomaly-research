"""Offline tests for scripts/sunrise_counter.py (synthetic SIS maps and tables; no network)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.wcs import WCS

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import sunrise_counter as sc  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402

THETA_E = 5.0  # SIS Einstein radius at D_LS/D_S = 1, arcsec
PIX = 0.05


def _sis_model(half: float = 20.0) -> lensmodel.MapLensModel:
    n = int(round(2 * half / PIX)) + 1
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = list(sc.CENTRE)
    w.wcs.crpix = [(n + 1) / 2, (n + 1) / 2]
    w.wcs.cdelt = [-PIX / 3600.0, PIX / 3600.0]
    g = (np.arange(n) - (n - 1) / 2) * PIX
    x, y = np.meshgrid(g, g)  # +x West along +i, +y North along +j
    r = np.hypot(x, y)
    r[r == 0] = 1e-6
    return lensmodel.MapLensModel(
        THETA_E * x / r, THETA_E * y / r, w, sc.Z_LENS, sc.COSMO, centre=sc.CENTRE
    )


def test_predict_recovers_sis_pair(monkeypatch):
    model = _sis_model()
    grid = lensmodel.DeflectionGrid.compute(model, 15.0, 0.1)
    z = 3.0
    te = THETA_E * model.dls_ds(z)
    beta = np.array([1.0, 0.5])
    u = beta / np.hypot(*beta)
    a, b = beta + te * u, beta - te * u
    ra, dec = model.to_sky(np.array([a[0]]), np.array([a[1]]))
    # only image "2a" is catalogued; "2b" is the SIS counter-image
    monkeypatch.setattr(sc, "IMAGES", {"2a": (float(ra[0]), float(dec[0]))})
    pred = sc.predict(model, grid, "2", z)
    assert pred["src_scatter"] < 1e-6
    sys_xy = sc.system_xy(model, "2")
    known, counter, recovered, rms = sc.split_images(pred, sys_xy, 2.0)
    # the map's finite-difference core adds a demagnified central image (as real maps do)
    counter = counter[np.abs(counter["magnification"]) >= 0.5]
    assert len(known) == 1 and len(counter) == 1 and recovered.all()
    assert rms < 0.05
    assert np.hypot(counter["x"][0] - b[0], counter["y"][0] - b[1]) < 0.05
    mu_b = float(counter["magnification"][0])
    assert mu_b < 0  # odd parity inside the Einstein ring


def test_split_images_marks_unrecovered():
    imgs = Table({"x": [0.0, 10.0], "y": [0.0, 0.0], "magnification": [5.0, -2.0]})
    pred = {"images": imgs, "ref_xy": np.array([[0.1, 0.0], [0.0, 4.0]])}
    known, counter, recovered, rms = sc.split_images(pred, pred["ref_xy"], 2.0)
    assert list(known["x"]) == [0.0] and list(counter["x"]) == [10.0]
    assert recovered.tolist() == [True, False]
    assert np.isclose(rms, 0.1)  # the missed image does not widen its own copy radius


def _res(model: str, pts):
    return {
        "model": model,
        "counter_images": [
            {"knot": "1.1", "ra": ra, "dec": dec, "class": cls, "inject_class": "recovered"}
            for ra, dec, cls in pts
        ],
    }


def test_cross_model_flags_only_images_every_model_predicts():
    d = 1.0 / 3600.0
    results = [
        _res("a", [(10.0, 0.0, "absent"), (10.1, 0.0, "absent")]),
        _res("b", [(10.0 + 0.5 * d, 0.0, "undetectable")]),
        _res("c", [(10.0, 0.5 * d, "confused"), (10.1 + 10 * d, 0.0, "absent")]),
    ]
    sc.cross_model(results, 3.0)
    a = results[0]["counter_images"]
    assert [c["n_models"] for c in a] == [3, 1]
    assert [c["flag"] for c in a] == [True, False]
    assert not any(c["flag"] for r in results[1:] for c in r["counter_images"])


def test_cross_model_ignores_demagnified_and_single_model_runs():
    d = 1.0 / 3600.0
    results = [
        _res("a", [(10.0, 0.0, "absent")]),
        _res("b", [(10.0 + 0.5 * d, 0.0, "demagnified")]),
    ]
    sc.cross_model(results, 3.0)
    assert results[0]["counter_images"][0]["n_models"] == 1
    assert not results[0]["counter_images"][0]["flag"]
    alone = [_res("a", [(10.0, 0.0, "absent")])]
    sc.cross_model(alone, 3.0)
    assert not alone[0]["counter_images"][0]["flag"]


def _noise_stamp(sigma=1.0, nan=False):
    g = (np.arange(-60, 61)) * 0.03  # 0.03" pixels, +-1.8"
    xx, yy = np.meshgrid(g, g)
    rng = np.random.default_rng(1)
    sci = rng.normal(0.0, sigma, xx.shape)
    err = np.full(xx.shape, sigma / sc.lc.ERR_SCALE)
    if nan:
        sci[:] = np.nan
    return lambda ra, dec: (sci, err, xx, yy)


def test_classify_absent_and_injection_recovers():
    v = sc.classify(_noise_stamp(), 0.0, 0.0, 200.0, 1.0)
    assert v["class"] == "absent"
    assert v["inject_class"] == "recovered"
    assert sc.classify(_noise_stamp(nan=True), 0.0, 0.0, 200.0, 1.0)["class"] == "off_image"
    assert sc.classify(_noise_stamp(), 0.0, 0.0, float("nan"), 1.0)["class"] == "no_reference"


def test_read_block_returns_the_requested_rows(monkeypatch):
    arr = np.arange(7 * 5, dtype=">f4").reshape(7, 5)
    offset = 2880
    blob = b"\0" * offset + arr.tobytes()
    monkeypatch.setattr(sc, "_get", lambda session, url, a, b: blob[a:b])
    monkeypatch.setattr(sc, "CHUNK", 12)  # several range requests, not aligned to rows
    out = sc.read_block(None, "u", {"NAXIS1": 5, "BITPIX": -32}, offset, (2, 6))
    assert np.array_equal(out, arr[2:6])
