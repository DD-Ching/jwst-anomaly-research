"""Offline test for scripts/s1_blind_spots.py: TTE light curve from a synthetic event file."""

import importlib.util
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))
_spec = importlib.util.spec_from_file_location(
    "s1_blind_spots", ROOT / "scripts" / "s1_blind_spots.py"
)
bs = importlib.util.module_from_spec(_spec)
sys.modules["s1_blind_spots"] = bs
_spec.loader.exec_module(bs)
from jwst_anomaly import burst_twins as bt  # noqa: E402


def _tte(path, rng):
    trig = 5.0e8
    # flat background 1,000 counts/s over [-30, 40] s plus two 20-ms pulses at 0.10 and 0.30 s
    bkg = rng.uniform(-30, 40, size=70_000)
    p1 = rng.normal(0.10, 0.01, size=1500)
    p2 = rng.normal(0.30, 0.01, size=1500)
    t = np.sort(np.r_[bkg, p1, p2]) + trig
    pha = rng.integers(0, 128, size=t.size)
    emin = np.geomspace(8.0, 1000.0, 129)[:-1]
    emax = np.geomspace(8.0, 1000.0, 129)[1:]
    h0 = fits.PrimaryHDU()
    h0.header["TRIGTIME"] = trig
    eb = fits.BinTableHDU.from_columns(
        [
            fits.Column("CHANNEL", "I", array=np.arange(128)),
            fits.Column("E_MIN", "E", array=emin),
            fits.Column("E_MAX", "E", array=emax),
        ],
        name="EBOUNDS",
    )
    ev = fits.BinTableHDU.from_columns(
        [fits.Column("TIME", "D", array=t), fits.Column("PHA", "I", array=pha)], name="EVENTS"
    )
    fits.HDUList([h0, eb, ev]).writeto(path)


def test_tte_lightcurve_resolves_two_pulses(tmp_path):
    f = tmp_path / "tte.fit"
    _tte(f, np.random.default_rng(0))
    flux, err, dt = bs.tte_lightcurve([f], 0.05, 0.3)
    assert dt <= 0.004
    assert len(bt.find_pulses(flux, err, bt.Params())) == 2
