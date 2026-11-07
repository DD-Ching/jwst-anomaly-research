import sys
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import miri_confirm as m  # noqa: E402


def test_verdict_follows_plan_gate():
    assert m.verdict(0.4, 0.05) == "adopt"
    assert m.verdict(0.4, 0.09) == "not adopted: chance rate too high"
    assert m.verdict(0.0, 0.0) == "not adopted: no real matches"
    assert m.verdict(float("nan"), 0.0) == "not adopted: no real matches"


def test_nircam_detected_needs_snr_in_one_band():
    t = Table()
    for b in m.NIRCAM_BANDS:
        t[f"{b}_flux_aper_1"] = [1.0, 1.0, np.nan]
        t[f"{b}_fluxerr_aper_1"] = [1.0, 1.0, 1.0]
    t["f444w_flux_aper_1"][1] = 6.0
    assert list(m.nircam_detected(t)) == [False, True, False]


def test_match_fraction_real_and_shifted():
    ref = SkyCoord([10.0] * u.deg, [-70.0] * u.deg)
    hit, sep = m.match_fraction(np.array([10.0]), np.array([-70.0]), ref, 0.3)
    assert hit[0] and sep[0] < 1e-6
    hit, sep = m.match_fraction(np.array([10.0]), np.array([-70.0]), ref, 0.3, dra_arcsec=5.0)
    assert not hit[0] and abs(sep[0] - 5.0) < 1e-3
    hit, sep = m.match_fraction(np.array([]), np.array([]), ref, 0.3)
    assert hit.size == 0
