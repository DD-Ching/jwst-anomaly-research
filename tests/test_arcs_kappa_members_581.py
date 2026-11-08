import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import arcs_kappa_members_581 as akm  # noqa: E402


def _map(peak_xy=None, amp=5.0):
    k = np.ones((100, 100))
    if peak_xy is not None:
        k[peak_xy[1], peak_xy[0]] = amp
    return k


def test_peak_found_on_position():
    st = akm.peak_at(_map((50, 50)), 50.0, 51.0, 0.1)
    assert st["peak"] and st["dist_arcsec"] == pytest.approx(0.1)


def test_peak_too_far_or_weak():
    assert not akm.peak_at(_map((50, 56)), 50.0, 50.0, 0.1)["peak"]  # 0.6 arcsec away
    assert not akm.peak_at(_map((50, 50), amp=1.2), 50.0, 50.0, 0.1)["peak"]


def test_edge_skipped():
    assert akm.peak_at(_map(), 5.0, 50.0, 0.1) is None


def test_gate_and_label():
    g = akm.gate([True] * 8 + [False] * 2, [False] * 9 + [True])
    assert g["passes"]
    assert not akm.gate([True] * 7, [False] * 9)["passes"]
    assert not akm.gate([True] * 6 + [False] * 4, [False] * 10)["passes"]
    assert not akm.gate([True] * 10, [True] * 2 + [False] * 8)["passes"]
    assert akm.label(True, True) == "GLAFIC member (model input)"
    assert akm.label(True, False) == "GLAFIC background (model input)"
    assert akm.label(False, True) == "undetermined"


def test_classes():
    assert akm.z_class(0.43) == "member"
    assert akm.z_class(0.48) == "background"
    assert akm.z_class(0.30) == "foreground"
    assert akm.spectroscopic("SLS", "G") and not akm.spectroscopic("SLS?", "G")
    assert not akm.spectroscopic("PUN", "G_Lens")
