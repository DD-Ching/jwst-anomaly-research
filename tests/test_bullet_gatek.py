import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import bullet_gatek as bg  # noqa: E402

CHIRES = """chi multiples
 N    ID    z   Narcs    chip    chix    chiy    chia   rmss     rmsi    dx      dy    nwarn
 0  N11.2 0.897   1      1.54    0.00    0.00    0.00   0.126    0.62    0.30    0.40  0
 0  N11.1 0.897   1      1.31    0.00    0.00    0.00   0.236    0.57   -0.60   -0.80  0
 0  N11.2 0.897   2      2.85    0.00    0.00    0.00   0.189    0.60    N/A     N/A   0
chimul                  92.10    0.00    0.00    0.00   0.166    0.41    N/A     N/A   0
log(Likelihood)       -107.02
"""


def test_parse_chires():
    r = bg.parse_chires(CHIRES)
    assert r["total"]["chi2"] == pytest.approx(92.10)
    assert r["total"]["rms_image"] == pytest.approx(0.41)
    assert r["total"]["n_images"] == 2  # the Narcs = 2 system row is not an image
    assert r["images"]["N11.1"] == (-0.60, -0.80)
    assert r["total"]["rms_from_rows"] == pytest.approx((0.5**2 / 2 + 1.0**2 / 2) ** 0.5)


def test_parse_chires_needs_total():
    with pytest.raises(ValueError):
        bg.parse_chires(CHIRES.replace("chimul", "other"))


def test_parse_chires_na_offsets():
    text = CHIRES.replace("-0.60   -0.80", "N/A     N/A")
    r = bg.parse_chires(text)
    assert r["total"]["n_images"] == 2
    assert r["total"]["rms_from_rows"] == pytest.approx(0.5)
