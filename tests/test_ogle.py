"""Offline tests for the OGLE-IV Mróz-sample adapter on tiny fixture files."""

from __future__ import annotations

import hashlib
import io
import tarfile

import numpy as np
import pytest

from jwst_anomaly import ogle
from jwst_anomaly.signatures import LightCurveSurvey

TABLE3 = """\
# Table 3. Best-fitting parameters (fixture: two rows of the published layout)
#    1- 16 S16    name     Star name
#   18- 26 S9     field    OGLE-IV subfield name
#   28- 33 I6     star_id  Star ID in the OGLE-IV database
#   35- 36 I2     ra_h     Right ascension (J2000) (hours)
#   38- 39 I2     ra_m     Right ascension (J2000) (minutes)
#   41- 45 F5.2   ra_s     Right ascension (J2000) (seconds)
#   47- 49 I3     dec_d    Declination (J2000) (degrees)
#   51- 52 I2     dec_m    Declination (J2000) (minutes)
#   54- 57 F4.1   dec_s    Declination (J2000) (seconds)
#   59- 67 F9.5   ra       Right ascension (J2000) (deg)
#   69- 77 F9.5   dec      Declination (J2000) (deg)
#   79- 87 F9.5   glon     Galactic longitude (deg)
#   89- 97 F9.5   glat     Galactic latitude (deg)
#   99-109 F11.3  t0_best  Best-fit t_0 (HJD)
#  111-117 F7.3   tE_best  Best-fit t_E (d)
#  119-123 F5.3   u0_best  Best-fit u_0
#  125-130 F6.3   Is_best  Best-fit I_s (mag)
#  132-136 F5.3   fs_best  Best-fit f_s
#  138-148 F11.3  t0_med   Median t_0 (HJD)
#  150-156 F7.3   t0_err1  Negative error bar t_0
#  158-163 F6.3   t0_err2  Positive error bar t_0
#  165-171 F7.3   tE_med   Median t_E (d)
#  173-179 F7.3   tE_err1  Negative error bar t_E (d)
#  181-186 F6.3   tE_err2  Positive error bar t_E (d)
#  188-192 F5.3   u0_med   Median u_0
#  194-199 F6.3   u0_err1  Negative error bar u_0
#  201-205 F5.3   u0_err2  Positive error bar u_0
#  207-212 F6.3   Is_med   Median I_s (mag)
#  214-219 F6.3   Is_err1  Negative error bar I_s (mag)
#  221-225 F5.3   Is_err2  Positive error bar I_s (mag)
#  227-231 F5.3   fs_med   Median f_s
#  233-238 F6.3   fs_err1  Negative error bar f_s
#  240-244 F5.3   fs_err2  Positive error bar f_s
#  246-251 F6.2   weight   Inverse of the detection efficiency
#  253-270 S18    ews_id   OGLE EWS ID
BLG617.16.73378  BLG617.16  73378 17:13:08.00 -29:48:13.0 258.28333 -29.80361  -4.63170   5.42560 2455434.908  24.294 0.206 20.104 0.911 2455434.908  -0.218  0.229  24.414  -2.846  4.555 0.206 -0.045 0.041 20.105 -0.234 0.308 0.910 -0.224 0.219   4.10 X
BLG617.24.41328  BLG617.24  41328 17:13:54.30 -29:36:35.5 258.47625 -29.60986  -4.37550   5.40191 2457462.780 1201.029 0.138 20.298 0.245 2457462.791  -0.657  0.572 199.662 -17.422 18.486 0.139 -0.016 0.019 20.286 -0.156 0.144 0.248 -0.031 0.038   2.80 OGLE-2016-BLG-0231
"""  # noqa: E501 (published row layout)

TABLE7 = """\
#    1-  6 S6     field        OGLE-IV subfield name
#    8- 15 F8.4   glon         Galactic longitude (deg)
#   17- 24 F8.4   glat         Galactic latitude (deg)
#   26- 29 F4.2   tau          Optical depth (10^{-6})
#   31- 34 F4.2   tau_err      Uncertainty of tau (10^{-6})
#   36- 39 F4.1   gam          Event rate per star (10^{-6} per year)
#   41- 44 F4.1   gam_err      Uncertainty of gam (10^{-6} per year)
#   46- 50 F5.1   gam_deg2     Event rate per unit area (per year per deg^2)
#   52- 55 F4.1   gam_deg2_err Uncertainty of gam_deg2 (per year per deg^2)
#   57- 60 F4.1   t_E_mean     Mean Einstein timescale (d)
#   62- 65 F4.1   t_E_mean_err Uncertainty of t_E_mean (d)
#   67- 69 I3     N_events     Number of events
#   71- 75 F5.2   N_stars      Number of sources (10^6)
BLG500   0.9999  -1.0293 1.93 0.21 23.9  2.0 168.8 13.7 18.8  1.6 164  6.78
BLG617  -4.2100   4.9609 0.70 0.14  5.2  1.1  20.9  4.4 30.8  7.0  43  6.32
"""

EFF = """\
# Detection efficiencies for field BLG617
# log_tE_min log_tE_max efficiency
0.000 0.125 0.010000
0.125 0.250 0.020000
0.250 0.375 0.040000
"""

PHOT = "2455642.79215 20.151 0.215\n2455640.81799 20.082 0.211\n2455649.81597 19.737 0.121\n"


def _tar(members: dict[str, str]) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        for name, text in members.items():
            data = text.encode()
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return buf.getvalue()


@pytest.fixture
def sample(tmp_path, monkeypatch):
    files = {
        "table3.dat": TABLE3.encode(),
        "table7.dat": TABLE7.encode(),
        "eff.tar.gz": _tar({"eff/BLG617.eff": EFF, "eff/BLG500.eff": EFF}),
        "phot.tar.gz": _tar({"phot/BLG617.16.73378.dat": PHOT}),
    }
    pins = dict(ogle.FILES)
    base = ogle.SAMPLES["bulge2019"].base
    for name, data in files.items():
        sha = hashlib.sha256(data).hexdigest()
        (tmp_path / f"{sha[:12]}_{name}").write_bytes(data)  # pre-cached: no download
        pins[base + name] = (sha, len(data))
    monkeypatch.setattr(ogle, "FILES", pins)
    return ogle.OgleMrozSample("bulge2019", cache_dir=tmp_path)


def test_parse_byte_table_handles_sexagesimal_and_overflowing_columns():
    t = ogle.parse_byte_table(TABLE3)
    assert len(t) == 2 and t["ra_h"][0] == 17 and t["dec_d"][1] == -29
    assert t["tE_best"][1] == pytest.approx(1201.029)  # wider than its F7.3 byte range
    assert t["ews_id"][1] == "OGLE-2016-BLG-0231"
    with pytest.raises(ValueError, match="no byte-by-byte"):
        ogle.parse_byte_table("1 2 3\n")


def test_efficiency_bins_and_out_of_grid():
    lo, hi, eff = ogle.parse_efficiency(EFF)
    got = ogle.eff_at(lo, hi, eff, [1.1, 1.5, 2.0, 0.5, 10.0])
    assert list(got) == [0.01, 0.02, 0.04, 0.0, 0.0]


def test_adapter_satisfies_the_protocol(sample):
    assert isinstance(sample, LightCurveSurvey)
    ev = sample.events()
    assert list(ev["event_id"]) == ["BLG617.16.73378", "BLG617.24.41328"]
    assert ev["field"][0] == "BLG617" and ev["tE_pub"][0] == pytest.approx(24.294)
    assert ev.meta["provenance"] == "observed"
    lc = sample.light_curve("BLG617.16.73378")
    assert np.all(np.diff(lc["time"]) > 0) and lc.meta["time_system"] == "HJD"
    with pytest.raises(KeyError, match="no photometry"):
        sample.light_curve("BLG617.24.41328")
    f = sample.fields()
    assert list(f["field"]) == ["BLG617"]  # high-cadence BLG500 excluded
    assert f["n_sources"][0] == pytest.approx(6.32e6)
    assert sample.efficiency(1.5) == pytest.approx(0.02)
    expected = 6.32e6 * 0.02 * 2741.0 / 365.25
    assert sample.exposure_star_years(1.5) == pytest.approx(expected)
