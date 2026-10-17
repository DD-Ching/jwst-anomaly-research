import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import w1_stpsf_floor as w  # noqa: E402


def test_isophotal_chi_elongated_gaussian():
    y, x = np.mgrid[-40:41, -40:41].astype(float)
    im = np.exp(-0.5 * (x**2 / 9.0 + y**2 / 4.0))
    c = w.isophotal_chi(im, 1e-6)
    assert abs(c.imag) < 1e-9
    assert np.isclose(c.real, (9 - 4) / (9 + 4), atol=1e-3)


def test_pixel_chi_to_sky_identity_and_rotation():
    # x = West, y = North (PC = [[-1, 0], [0, 1]]): an x-elongated χ has PA 90°
    sky = w.pixel_chi_to_sky(0.2 + 0j, pc=((-1.0, 0.0), (0.0, 1.0)))
    assert np.isclose(abs(sky), 0.2)
    assert np.isclose(np.mod(0.5 * np.degrees(np.angle(sky)), 180), 90.0)
    # the Abell 2744 mosaic: the pixel x axis points to PA ≈ 130.8°
    sky = w.pixel_chi_to_sky(0.1 + 0j)
    assert np.isclose(np.mod(0.5 * np.degrees(np.angle(sky)), 180), 130.8, atol=0.1)


def test_s1_consistent_per_component():
    assert w.s1_consistent(complex(-0.015, -0.156))
    assert w.s1_consistent(complex(0.03, -0.11))
    assert not w.s1_consistent(complex(-0.015, -0.10))  # imag differs by 0.056
    assert not w.s1_consistent(complex(0.04, -0.156))  # real differs by 0.055
