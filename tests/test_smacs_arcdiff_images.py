import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from smacs_arcdiff_images import classify  # noqa: E402


def test_known_wins_over_verdict():
    assert classify("dipole", 5.0, 0.05) == "known_not_pass"
    assert classify("pass", 5.0, 0.05) == "known"
    assert classify("pass", 0.2, 0.15) == "known"


def test_pass_on_image_is_vetted():
    assert classify("pass", 1.0, 3.0) == "vet"
    assert classify("pass", 1.01, 3.0) == "off_image"


def test_non_pass_is_vetoed():
    assert classify("edge", 0.1, 3.0) == "vetoed"
