import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import w2_desi_pairs as wd  # noqa: E402


def _spec(z, rng, *, ew=1.0, shift_kms=0.0, width_kms=2500.0, noise=0.02, line="MgII"):
    wave = np.arange(3600.0, 9824.0, 0.8)
    rest = wd.LINES[line]
    v = wd.C_KMS * (wave / (rest * (1 + z)) - 1)
    cont = 10.0 * (wave / 6000.0) ** -0.5
    prof = ew * 0.8 * np.exp(-0.5 * ((v - shift_kms) / width_kms) ** 2)
    flux = cont * (1 + prof)
    sig = noise * cont
    return wave, flux + rng.normal(0, sig), np.full(len(wave), 1 / sig**2)


def test_same_source_with_scaled_line_is_consistent():
    rng = np.random.default_rng(1)
    w, fa, iva = _spec(1.0, rng)
    _, fb, ivb = _spec(1.0, rng, ew=0.7)  # microlensing-like equivalent-width change, same profile
    st = wd.pair_statistic(w, 0.6 * fa, iva / 0.36, fb, ivb, 1.0)  # flux ratio 0.6
    assert st["lines"]["MgII"]["used"]
    assert st["D"] < 1.5
    assert abs(st["dv_kms"]) < 100


def test_velocity_offset_or_profile_change_is_two_spectra():
    rng = np.random.default_rng(2)
    w, fa, iva = _spec(1.0, rng)
    _, fb, ivb = _spec(1.0, rng, shift_kms=800.0)
    st = wd.pair_statistic(w, fa, iva, fb, ivb, 1.0)
    assert st["D"] > wd.T_D
    assert abs(abs(st["dv_kms"]) - 800) < 150
    _, fc, ivc = _spec(1.0, rng, width_kms=4500.0)
    assert wd.pair_statistic(w, fa, iva, fc, ivc, 1.0)["D"] > wd.T_D


def test_lines_outside_coverage_or_faint_are_not_counted():
    rng = np.random.default_rng(3)
    w, fa, iva = _spec(1.0, rng, ew=0.0)  # no line at all: line S/N below 10
    _, fb, ivb = _spec(1.0, rng, ew=0.0)
    st = wd.pair_statistic(w, fa, iva, fb, ivb, 1.0)
    assert st["D"] is None
    assert not st["lines"]["MgII"]["used"]
    assert "Halpha" not in st["lines"]  # 13129 A is outside DESI


def test_best_pair_matches_catalogued_separation():
    rows = [
        {"mean_fiber_ra": "10.0", "mean_fiber_dec": "0.0"},
        {"mean_fiber_ra": str(10.0 + 2.0 / 3600), "mean_fiber_dec": "0.0"},
        {"mean_fiber_ra": str(10.0 + 3.4 / 3600), "mean_fiber_dec": "0.0"},
    ]
    best = wd.best_pair(rows, 2.05)
    assert best is not None and abs(best[3] - 2.0) < 0.01
    assert wd.best_pair(rows, 5.0) is None
