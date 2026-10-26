import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import d115b_lenstool_refit as m  # noqa: E402

BEST = """#Chi2pos: 344.298308
runmode
\treference     3 64.038217 -24.067501
\timage     1 gold.dat
\tend
potentiel O1
\tprofil       81
\tv_disp     574.4
\tend
fini
"""
INPUT = """runmode
    reference 3 64.0 -24.0
    inverse   3 0.015 10000
    end

image
    multfile    1 gold.dat
    forme       -1
    sigposArcsec\t0.49
    end

finish
"""


def test_build_v0_replaces_runmode_and_adds_image():
    out = m.build_v0(BEST, INPUT)
    lines = [ln.strip() for ln in out.splitlines()]
    i = lines.index("runmode")
    assert lines[i + 1].split() == ["reference", "3", "64.038217", "-24.067501"]
    assert lines[i + 2 : i + 5] == ["inverse 0", "chi2 1", "end"]
    assert lines[i + 5] == "image" and "forme       -1" in lines
    assert "multfile    1 arcs.dat" in lines and "multfile    1 gold.dat" not in lines
    assert "sigposArcsec\t0.49" in lines and lines.count("end") == 3
    assert "image     1 gold.dat" not in lines  # the old runmode body is dropped
    assert "v_disp     574.4" in lines and lines[-1] == "fini"
    assert out.splitlines()[0] == "#Chi2pos: 344.298308"


def test_build_v0_needs_one_reference():
    with pytest.raises(ValueError):
        m.build_v0(BEST.replace("\treference     3 64.038217 -24.067501\n", ""), INPUT)


def test_build_v0_needs_one_multfile():
    with pytest.raises(ValueError):
        m.build_v0(BEST, INPUT.replace("    multfile    1 gold.dat\n", ""))


def test_read_chi2pos_takes_last():
    log = "INFO: chi2 1.500000  Lhood -1\nfoo\nINFO: chi2 344.300000  Lhood -408.6\n"
    assert m.read_chi2pos(log) == pytest.approx(344.30)
    with pytest.raises(ValueError):
        m.read_chi2pos("nothing here")


def test_gate_v():
    assert m.gate_v(344.298308)
    assert m.gate_v(344.298308 * 1.0099)
    assert not m.gate_v(344.298308 * 1.0101)
    assert not m.gate_v(330.81)


CHIRES = """chi multiples
 N    ID    z   Narcs    chip    chix    chiy    chia   rmss     rmsi    dx      dy    nwarn
 0  A.2 0.940   1      0.01    0.00    0.00    0.00   0.008    0.04   -0.02    0.04  0
 0  A.1 0.940   1      0.01    0.00    0.00    0.00   0.008    0.05   -0.04    0.03  1
 0  A.2 0.940   2      0.02    0.00    0.00    0.00   0.008    0.05    N/A     N/A   0
 1  D115a.1 2.300   1      0.99    0.00    0.00    0.00   0.044    0.49    0.38   -0.30  0
 1  D115a.3 2.300   1      4.60    0.00    0.00    0.00   0.060    0.61    N/A     N/A   0
chimul                 4.62    0.00    0.00    0.00   0.211    0.53    N/A     N/A   0
chitot           9.99
"""


def test_read_chires():
    chimul, fam, nwarn = m.read_chires(CHIRES)
    assert chimul == pytest.approx(4.62)
    assert fam == {"A": pytest.approx(0.02), "D115a": pytest.approx(4.60)}
    assert nwarn == 1
    with pytest.raises(ValueError):
        m.read_chires("chitot 1.0\n")


BEST2 = BEST.replace(
    "fini\n", "potentiel 8515\n\tprofil       81\n\tv_disp     120.600977\n\tend\nfini\n"
)
INPUT2 = INPUT.replace("finish\n", "limit 1\n    v_disp      1 350. 1000. 0.1\n    end\nfinish\n")


def test_build_vprime_swaps_only_the_reference():
    v0 = m.build_v0(BEST, INPUT)
    vp = m.build_vprime(v0, BEST, INPUT)
    assert "reference 3 64.0 -24.0" in vp and "64.038217" not in vp
    assert vp.replace("    reference 3 64.0 -24.0", "\treference     3 64.038217 -24.067501") == v0


def test_set_vdisp_and_image_mode():
    v0 = m.build_v0(BEST2, INPUT)
    par = m.add_image_mode(m.set_vdisp(v0, "8515", 0.75 * m.SIGMA_8515))
    lines = [ln.strip() for ln in par.splitlines()]
    assert "v_disp     90.450733" in lines and "v_disp     574.4" in lines
    assert lines[lines.index("inverse 0") + 1] == "image 1 arcs.dat"
    with pytest.raises(ValueError):
        m.set_vdisp(v0, "9999", 1.0)


def test_rule_f():
    fs = m.F_GRID
    flag = {f: f >= 0.80 for f in fs}
    d = {f: 0.0 for f in fs}
    assert m.rule_f({**d, 0.75: 6.0}, flag) == ("cheap", 0.75)
    assert m.rule_f({**d, 0.75: 6.01}, flag) == ("costly", 0.75)
    assert m.rule_f(d, {f: True for f in fs}) == ("undetermined", None)
    assert m.rule_f(d, {**flag, 1.0: False}) == ("undetermined", None)
    # f_L* is the largest flag-free factor, even if the flag returns below it
    assert m.rule_f({**d, 0.85: 9.0}, {**flag, 0.85: False}) == ("costly", 0.85)


def test_build_f2_frees_halos_and_8515():
    best = BEST2.replace("fini\n", "").replace(
        "potentiel 8515",
        "potentiel O2\n\tend\npotentiel O3\n\tend\npotentiel O4\n\tend\npotentiel 8515",
    )
    best += "fini\n"
    inp = INPUT2.replace(
        "finish\n",
        "".join(f"limit {i}\n    v_disp      1 1. 2. 0.1\n    end\n" for i in (2, 3, 4))
        + "finish\n",
    )
    vp = m.build_vprime(m.build_v0(best, inp), best, inp)
    out = m.build_f2(vp, best, inp)
    lines = [ln.strip() for ln in out.splitlines()]
    assert "inverse   3 0.015 10000" in lines and "inverse 0" not in lines
    assert [ln for ln in lines if ln.startswith("limit")] == [f"limit {i}" for i in (1, 2, 3, 4, 5)]
    assert lines[lines.index("limit 5") + 1] == "v_disp      1 60. 200. 0.1" and lines[-1] == "fini"


def test_flag_images(tmp_path):
    ra0, dec0 = 64.03821668520617, -24.06750122335985
    fx = -(64.048551 - ra0) * np.cos(np.radians(dec0)) * 3600
    fy = (-24.064825 - dec0) * 3600
    (tmp_path / "image.all").write_text(
        f"#REFERENCE 3 {ra0} {dec0}\nD115b.1 {fx + 0.6} {fy} 0\nD115b.1 0 0 0\n"
        f"D115b.2 {fx + 0.8} {fy} 0\n"
    )
    assert m.flag_images(tmp_path / "image.all") == {"D115b.1": True, "D115b.2": False}
