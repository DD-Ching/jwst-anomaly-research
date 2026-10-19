import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from dimming_exposure_check import chromatic_veto  # noqa: E402


def test_chromatic_veto_flags_sw_only_point():
    assert chromatic_veto(4.27, 0.0, 0.01)  # 4.27 / 20 = 0.21 > 0.03


def test_chromatic_veto_keeps_flat_spectrum_source():
    assert not chromatic_veto(0.16, 0.057, 0.005)  # ratio 2.8, ordinary colour


def test_chromatic_veto_respects_lw_noise():
    assert not chromatic_veto(1.0, 0.0, 0.02)  # 0 + 3 * 0.02 = 0.06 > 0.05


def test_chromatic_veto_not_applicable():
    assert not chromatic_veto(np.nan, 0.0, 0.01)
    assert not chromatic_veto(-1.0, 0.0, 0.01)
    assert not chromatic_veto(1.0, np.nan, 0.01)
