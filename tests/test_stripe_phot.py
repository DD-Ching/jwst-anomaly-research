import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from stripe_phot import classify_excess, perp_along, stripe_angle, stripe_model  # noqa: E402


def test_classify_excess():
    assert classify_excess(41.0, 3.0, [40.0, 45.0], [3.0, 5.0], 1.0) == "explained"
    assert classify_excess(76.0, 2.2, [40.0, 50.0], [2.6, 4.6], 1.0) == "open"
    assert classify_excess(50.0, 2.0, [42.0], [2.0], 1.0) == "undetermined"
    assert classify_excess(41.0, 3.0, [40.0], [3.0], 0.5) == "invalid"
    assert classify_excess(41.0, 3.0, [], [], 1.0) == "invalid"
    assert classify_excess(float("nan"), 3.0, [40.0], [3.0], 1.0) == "invalid"
    assert classify_excess(41.0, 3.0, [40.0], [float("nan")], 1.0) == "invalid"


def test_stripe_found_and_removed():
    scale = 0.06
    yy, xx = np.mgrid[-50:51, -50:51] * scale
    e, n = -xx, yy  # east to the left, as on the sky
    pa = 37.0
    d, _ = perp_along(e, n, pa)
    rng = np.random.default_rng(0)
    img = 5.0 * np.exp(-(d**2) / (2 * 0.1**2)) + rng.normal(0, 0.1, e.shape)
    got, _ = stripe_angle(img, e, n, scale)
    assert abs(got - pa) <= 1.0
    resid = img - np.nan_to_num(stripe_model(img, e, n, got, scale))
    assert abs(np.median(resid[np.hypot(e, n) < 0.3])) < 0.1
