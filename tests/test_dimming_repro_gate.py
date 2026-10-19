import sys
from pathlib import Path

from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from dimming_repro_gate import repro_gate  # noqa: E402

REF = {0.1: 10.0, 1.0: 50.0}


def _lim(values, calibrated=True, rho=0.1):
    t = Table(
        {
            "t_e_yr": list(values),
            "rho": [rho] * len(values),
            "exposure_source_yr": list(values.values()),
        }
    )
    t.meta["calibrated"] = calibrated
    return t


def test_passes_within_tolerance():
    assert repro_gate(_lim({0.1: 10.9, 1.0: 46.0}), REF)["passed"]


def test_fails_when_one_t_e_outside():
    res = repro_gate(_lim({0.1: 10.0, 1.0: 44.0}), REF)
    assert not res["passed"] and [r["within"] for r in res["rows"]] == [True, False]


def test_fails_when_uncalibrated():
    assert not repro_gate(_lim({0.1: 10.0, 1.0: 50.0}, calibrated=False), REF)["passed"]


def test_ignores_other_rho_rows():
    assert not repro_gate(_lim({0.1: 10.0, 1.0: 50.0}, rho=0.01), REF)[
        "passed"
    ]  # no ρ = 0.1 rows → NaN
