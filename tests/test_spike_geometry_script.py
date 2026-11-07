import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import spike_geometry as sg  # noqa: E402


def _table(pa_streak: float, ell: float = 0.95, n_bands: int = 1, dmag: float = 5.0) -> Table:
    # Bright source at the origin; streak 4.6" due north of it (bearing from the streak = 180 deg).
    dec0 = 0.0
    t = Table()
    t["source_uid"] = ["x_f200w_1", "x_f200w_2"]
    t["ra"] = [10.0, 10.0]
    t["dec"] = [dec0, dec0 + 4.6 / 3600]
    t["n_bands"] = [3, n_bands]
    t["ref_band"] = ["f200w", "f200w"]
    t["anchor_band"] = ["f200w", "f200w"]
    t["f200w_detected"] = [True, True]
    t["f200w_aper_total_abmag"] = [20.0, 20.0 + dmag]
    t["f200w_isophotal_abmag"] = [20.0, 20.0 + dmag]
    t["f200w_ellipticity"] = [0.1, ell]
    t["f200w_sky_orientation"] = [0.0, pa_streak]
    return t


def _flag(t: Table, rotate: float = 0.0) -> bool:
    return sg.spike_flags(t, np.array([1]), np.array(["f200w"], dtype=object), rotate=rotate)[0]["flag"]


def test_axis_diff_is_mod_180():
    assert sg.axis_diff(5, 185) == 0
    assert sg.axis_diff(10, 170) == 20
    assert sg.axis_diff(0, 90) == 90


def test_radial_streak_is_flagged_and_rotation_null_is_not():
    t = _table(pa_streak=0.0)  # axis N-S, pointing at the bright source
    assert _flag(t)
    assert _flag(_table(pa_streak=189.0))  # within 10 deg as an axis
    assert not _flag(t, rotate=90.0)


def test_tangential_or_round_or_multiband_or_faint_neighbour_not_flagged():
    assert not _flag(_table(pa_streak=90.0))  # tangential, like an arc
    assert not _flag(_table(pa_streak=0.0, ell=0.7))
    assert not _flag(_table(pa_streak=0.0, n_bands=2))
    assert not _flag(_table(pa_streak=0.0, dmag=2.0))
