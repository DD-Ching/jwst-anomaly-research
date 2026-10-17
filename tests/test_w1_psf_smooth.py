import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import w1_psf_smooth as w  # noqa: E402


def test_smooth_adds_variances_in_quadrature():
    y, x = np.mgrid[-60:61, -60:61].astype(float)
    im = np.exp(-0.5 * (x**2 / 9.0 + y**2 / 4.0))
    out = w.smooth(im, (1.0, np.sqrt(3.0)))  # total σ² = 4
    c = w.wsf.isophotal_chi(out, 1e-8)
    assert np.isclose(c.real, (13 - 8) / (13 + 8), atol=2e-3)
    assert w.smooth(im, ()) is im


def test_s3_star_per_component():
    assert w.s3_star([0.05, -0.01], [0.10, -0.01])
    assert not w.s3_star([0.05, -0.01], [0.101, -0.01])
    assert not w.s3_star([0.15, 0.0], [0.15, 0.051])
