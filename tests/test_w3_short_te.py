"""Offline tests for scripts/w3_short_te.py (frozen dev / validation split)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
_spec = importlib.util.spec_from_file_location("w3_short_te", _DIR / "w3_short_te.py")
st = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(st)


def test_split_depends_on_the_injection_not_on_row_order():
    rows = [
        {
            "event_id": f"gb12-R-1-1-{i}",
            "tE": 3.0,
            "rho": 0.01,
            "u0": i / 997,
            "t0": 5000.0 + i,
            "Is": 19.0,
        }  # fmt: skip
        for i in range(2000)
    ]
    a = [st.split_of(r, 3003) for r in rows]
    b = [st.split_of(r, 3003) for r in reversed(rows)][::-1]
    assert a == b  # order-free and deterministic
    assert abs(a.count("dev") - 1000) < 100
    assert a != [st.split_of(r, 3004) for r in rows]  # the seed enters


def test_tracked_sample_uses_the_job_key():
    from astropy.table import Table

    p = st.OUT / "injections_gb12_te3_seed3003.ecsv.gz"
    tab = Table.read(p, format="ascii.ecsv")
    assert list(tab["split"]) == [st.split_of(r, 3003) for r in tab]


def _nightly(n_nights=300, per_night=6, s_night=0.0, seed=1):
    import numpy as np

    rng = np.random.default_rng(seed)
    nt = np.repeat(np.arange(n_nights), per_night)
    w = np.full(nt.size, 1.0 / 0.5**2)  # σ = 0.5 per exposure
    r = rng.normal(0.0, 0.5, nt.size) + rng.normal(0.0, s_night, n_nights)[nt]
    return r, w, nt


def test_night_chi2_reduces_to_white_noise_without_a_night_term():
    import numpy as np

    r, w, nt = _nightly()
    assert np.isclose(st.night_chi2(r, w, nt, 0.0), float((w * r * r).sum()))


def test_night_chi2_absorbs_a_single_offset_night():
    import numpy as np

    r, w, nt = _nightly()
    bad = r.copy()
    bad[nt == 7] += 5.0  # one night 10σ high in every exposure
    white = float((w * bad * bad).sum()) - float((w * r * r).sum())
    robust = st.night_chi2(bad, w, nt, 1.0) - st.night_chi2(r, w, nt, 1.0)
    assert white > 500.0
    assert robust < 0.1 * white  # the night term caps what one coherent night can contribute
    assert np.isfinite(robust)


def test_night_variance_recovers_known_cases():
    assert st.night_variance(*_nightly(s_night=0.0)) < 0.01  # white noise: ~0
    s2 = st.night_variance(*_nightly(s_night=0.4))
    assert 0.25**2 < s2 < 0.55**2  # injected σ_night = 0.4
