import sys
from pathlib import Path

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
