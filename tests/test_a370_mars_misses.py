import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import a370_mars_misses as am  # noqa: E402


def test_classify():
    assert am.classify(0.5, 2.9) == "both"
    assert am.classify(0.5, 3.1) == "MARS only"
    assert am.classify(20.9, 0.4) == "BUFFALO only"
    assert am.classify(np.inf, np.nan) == "neither"
    assert am.classify(3.0, 3.0) == "both"
