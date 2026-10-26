import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import macs0647_g85_host as h  # noqa: E402


def test_sep_and_host_selection():
    ra0, dec0 = h.G85
    dra = 0.4 / 3600 / np.cos(np.radians(dec0))
    sep = h.sep_arcsec([ra0 + dra, ra0, ra0 + 3 * dra], [dec0, dec0 + 0.2 / 3600, dec0])
    assert sep == pytest.approx([0.4, 0.2, 1.2], abs=1e-6)
    assert h.select_host(sep) == 1
    assert h.select_host(np.array([0.6, 1.0])) is None
    assert h.select_host(np.array([])) is None


def test_classify_photoz():
    assert h.classify_photoz(0.55, 0.60, 0.65) == "member"
    assert h.classify_photoz(0.70, 0.74, 0.80) == "member"  # |0.74 - 0.591| = 0.149 <= 0.159
    assert h.classify_photoz(0.70, 0.76, 0.80) == "background"  # 0.169 > 0.159
    assert h.classify_photoz(1.8, 2.0, 2.2) == "background"
    assert h.classify_photoz(0.2, 0.3, 0.4) == "foreground"
    assert h.classify_photoz(0.5, 1.0, 2.1) == "unconstrained"
    assert h.classify_photoz(np.nan, 1.0, 1.2) == "unconstrained"


def test_spec_class():
    rows = [{"grade": "3", "z_best": "2.10"}, {"grade": "2", "z_best": "0.59"}]
    assert h.spec_class(rows) == ("background", 2.10)
    assert h.spec_class([{"grade": "2", "z_best": "0.59"}]) == ("none", None)
    assert h.spec_class([{"grade": "3.0", "z_best": "0.60"}])[0] == "member"


def test_abs_mag():
    # DM(z=0.591) for H0 70, Om 0.3 is about 42.70; 2.5 log10(1.591) = 0.504
    m = h.abs_mag(26.28, 0.591)
    assert m == pytest.approx(26.28 - 42.70 + 0.504, abs=0.03)
