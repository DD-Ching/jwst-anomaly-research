import io
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import arcs_muse_2157 as m  # noqa: E402


def _fake_cube(nl=50, ny=12, nx=10):
    data = np.arange(nl * ny * nx, dtype=np.float32).reshape(nl, ny, nx)
    hl = fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(data, name="DATA")])
    buf = io.BytesIO()
    hl.writeto(buf)
    raw = buf.getvalue()
    return data, [raw[i : i + 1000] for i in range(0, len(raw), 1000)]


def test_stream_cube_boxes_and_white():
    data, chunks = _fake_cube()
    out = m.stream_cube(chunks, {"a": (2, 5, 3, 7)})
    np.testing.assert_array_equal(out["boxes"]["a"], data[:, 3:7, 2:5])
    np.testing.assert_allclose(out["white"], data.mean(axis=0), rtol=1e-6)


def test_scan_recovers_injected_lines():
    lam = 4750 + 1.25 * np.arange(3681)
    rng = np.random.default_rng(1)
    noise = np.ones_like(lam)
    flux = rng.normal(0, 1, lam.size)
    z = 0.756
    for rest in (3727.4, 4861.3, 5006.8):
        flux += 8 * np.exp(-0.5 * ((lam - rest * (1 + z)) / 2.0) ** 2)
    resid = m.continuum_subtract(lam, flux)
    scan = m.scan_redshift(lam, resid, noise, np.arange(0.0, 1.6, 0.0005))
    assert abs(scan["z_best"] - z) < 0.002
    feats = m.features_at(lam, resid, noise, scan["z_best"])
    assert m.n_features(feats) >= 2
    assert m.classify(scan["z_best"], m.n_features(feats)) == "background"


def test_classify_rules():
    assert m.classify(0.36, 2) == "member"
    assert m.classify(0.80, 3) == "background"
    assert m.classify(0.20, 2) == "foreground"
    assert m.classify(0.80, 1) == "single_line"
    assert m.classify(0.80, 0) == "undetermined"
