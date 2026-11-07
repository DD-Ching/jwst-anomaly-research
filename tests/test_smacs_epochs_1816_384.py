import importlib.util
import sys
from pathlib import Path

_path = Path(__file__).resolve().parents[1] / "scripts" / "smacs_epochs_1816_384.py"
sys.path.insert(0, str(_path.parent))
_spec = importlib.util.spec_from_file_location("smacs_epochs_1816_384", _path)
mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mod)


def _r(sig, total, vec, dm=0.0, dms=0.05):
    return {
        "target_residual_significance": sig,
        "target_residual_total_mas": total,
        "target_residual_mas": vec,
        "target_dmag": dm,
        "reference_dmag_sigma": dms,
    }


def test_missing_band_is_undetermined():
    d = mod.decide({"f150w": _r(10, 20, [20, 0]), "f444w": None})
    assert d == {"motion": "undetermined", "variability": "undetermined"}


def test_moves_needs_same_direction():
    assert (
        mod.decide({"f150w": _r(5, 20, [20, 0]), "f444w": _r(4, 18, [18, 2])})["motion"]
        == "moves: star"
    )
    assert (
        mod.decide({"f150w": _r(5, 20, [20, 0]), "f444w": _r(4, 18, [-18, 0])})["motion"]
        == "undetermined"
    )


def test_no_motion_needs_tight_bound():
    assert (
        mod.decide({"f150w": _r(1, 3, [3, 0]), "f444w": _r(1, 4, [0, 4])})["motion"]
        == "no motion detected"
    )
    # 1 sigma at 8 mas -> 3 x error 24 mas > 15 mas bound
    assert (
        mod.decide({"f150w": _r(1, 8, [8, 0]), "f444w": _r(1, 4, [0, 4])})["motion"]
        == "undetermined"
    )


def test_variability_both_bands():
    both = mod.decide({"f150w": _r(1, 3, [3, 0], dm=0.5), "f444w": _r(1, 3, [3, 0], dm=-0.4)})
    assert both["variability"] == "variable"
    one = mod.decide({"f150w": _r(1, 3, [3, 0], dm=0.5), "f444w": _r(1, 3, [3, 0], dm=0.0)})
    assert one["variability"] == "undetermined"
    none = mod.decide({"f150w": _r(1, 3, [3, 0]), "f444w": _r(1, 3, [3, 0])})
    assert none["variability"] == "not variable"


def _cats(n=40, target_shift_mas=50.0, seed=1):
    import astropy.units as u
    import numpy as np
    from astropy.coordinates import SkyCoord
    from astropy.table import Table

    rng = np.random.default_rng(seed)
    ra = 110.65 + rng.uniform(-0.01, 0.01, n)
    dec = -73.49 + rng.uniform(-0.003, 0.003, n)
    c1 = SkyCoord(ra, dec, unit="deg")
    # every source moves by a common 10 mas east frame offset plus 2 mas noise; row 0 also moves
    dra = 10.0 + rng.normal(0, 2.0, n)
    dra[0] += target_shift_mas
    c2 = c1.spherical_offsets_by(dra * u.mas, rng.normal(0, 2.0, n) * u.mas)
    cat1 = Table({"aper50_abmag": np.full(n, 25.7), "semimajor_sigma": np.full(n, 1.0)})
    return c1, c2, cat1


def test_empirical_null_uses_given_tie_and_drops_target():
    c1, c2, cat1 = _cats()
    out = mod.empirical_null(c1, c2, cat1, 0, 0, [10.0, 0.0])
    assert out["n_peers"] == 39  # the target's rows are excluded
    assert 45 < out["target_residual_mas"] < 55
    assert out["n_peers_with_larger_residual"] == 0
    assert out["peer_residual_median_mas"] < 6
    # a wrong tie shifts every residual, the target's included
    off = mod.empirical_null(c1, c2, cat1, 0, 0, [-40.0, 0.0])
    assert off["peer_residual_median_mas"] > 40


def test_zero_significance_gives_no_bound():
    d = mod.decide({"f150w": _r(0, 0, [0, 0]), "f444w": _r(1, 3, [3, 0])})
    assert d["motion"] == "undetermined"
