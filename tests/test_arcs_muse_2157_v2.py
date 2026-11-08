import io
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import arcs_muse_2157_v2 as m  # noqa: E402


def test_stream_data_stat_boxes():
    nl, ny, nx = 7, 9, 8
    data = np.arange(nl * ny * nx, dtype=np.float32).reshape(nl, ny, nx)
    stat = data * 2
    hl = fits.HDUList(
        [fits.PrimaryHDU(), fits.ImageHDU(data, name="DATA"), fits.ImageHDU(stat, name="STAT")]
    )
    buf = io.BytesIO()
    hl.writeto(buf)
    raw = buf.getvalue()
    chunks = (raw[i : i + 1000] for i in range(0, len(raw), 1000))
    out = m.stream_data_stat(chunks, {"a": (1, 4, 2, 6)})
    np.testing.assert_array_equal(out["DATA"]["boxes"]["a"], data[:, 2:6, 1:4])
    np.testing.assert_array_equal(out["STAT"]["boxes"]["a"], stat[:, 2:6, 1:4])


def test_recentre_finds_peak_within_radius():
    yy, xx = np.mgrid[0:40, 0:40]
    img = np.exp(-((xx - 22.0) ** 2 + (yy - 18.0) ** 2) / 4.0)
    x, y = m.recentre(img, 20.0, 20.0, 5.0)
    assert abs(x - 22.0) < 0.3 and abs(y - 18.0) < 0.3
    x, y = m.recentre(img, 5.0, 5.0, 3.0)  # peak outside the search radius: stays within it
    assert np.hypot(x - 5.0, y - 5.0) <= 3.0 + 1e-9


def test_peak_contrast_and_rules():
    z = np.arange(0, 1, 0.001)
    score = np.ones_like(z)
    score[356] = 3.0
    assert m.peak_contrast(z, score) == 3.0
    assert m.gate_m2(0.356, 3.0, 0.356) and not m.gate_m2(0.356, 1.4, 0.356)
    assert not m.gate_m2(0.37, 3.0, 0.356)
    assert m.classify_v2(0.36, 1.6) == "member"
    assert m.classify_v2(0.8, 1.6) == "background"
    assert m.classify_v2(0.2, 1.6) == "foreground"
    assert m.classify_v2(0.8, 1.49) == "undetermined"
    assert m.classify_v2(0.8, float("nan")) == "undetermined"


def test_oii_model_integrates_and_injection_recovers():
    lam = np.arange(4750.0, 9350.0, 1.25)
    model = m.oii_model(lam, 0.7, 500.0)
    assert abs(model.sum() * 1.25 - 500.0) < 1.0
    rng = np.random.default_rng(1)
    noise = np.full(lam.size, 5.0)
    spec = 100 + rng.normal(0, 5.0, lam.size)
    mask = np.zeros(lam.size, bool)
    zg = np.arange(0.0, 1.6, 0.0002)
    det, zid = m.recovered(lam, spec, noise, mask, zg, 0.7, 2000.0)
    assert det and not zid  # a lone line is z-degenerate under the scan (amendment 1)
    assert not m.recovered(lam, spec, noise, mask, zg, 0.7, 1.0)[0]
    assert m.f90([100, 50, 200], [0.95, 0.5, 1.0]) == 100
    assert m.f90([1, 2], [0.1, 0.2]) is None
