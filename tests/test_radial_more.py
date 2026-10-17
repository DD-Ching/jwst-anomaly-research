import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import radial_more as rm  # noqa: E402


def test_flagged_uses_strict_threshold():
    s = {"peaks": [{"p_random": 0.049}, {"p_random": 0.05}, {"p_random": None}, {"p_random": 0.2}]}
    assert rm.flagged(s) == [{"p_random": 0.049}]


class _Frame:
    def to_frame(self, ra, dec):
        return np.asarray(ra, float), np.asarray(dec, float)


def test_contributors_match_line_counts():
    # arcs at (0, 5) pointing N-S pass the origin; one at (5, 0) pointing N-S does not
    cand = Table({"ra": [0.0, 5.0, 0.0], "dec": [5.0, 0.0, -20.0], "pa_obs": [0.0, 0.0, 0.0]})
    args = SimpleNamespace(line_tol_arcsec=1.0, max_len_arcsec=15.0)
    use = rm.contributors(cand, _Frame(), {"x": 0.0, "y": 0.0}, args)
    assert use.tolist() == [True, False, False]  # the third is 20" away, beyond max_len
