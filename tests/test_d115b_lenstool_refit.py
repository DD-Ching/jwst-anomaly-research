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


def test_family_rows():
    assert m._family_rows(CHIRES, "A") == (0.02, 0.02)
    with pytest.raises(ValueError):
        m._family_rows(CHIRES, "Z")


def test_read_penalty_and_images():
    assert m.read_penalty(CHIRES) == 0.0
    assert m.read_penalty(CHIRES + "chipenalty 3748.44 355.70\n") == pytest.approx(3748.44)
    log = "x\nWARNING: chi2 penalty for image D14a.2 = 3748.443149\nWARNING: other\n"
    assert m.penalty_images(log) == ["D14a.2"]
    assert m.PENALTY == pytest.approx(3748.44, abs=0.005)


def test_rule_p():
    P = m.PENALTY
    ok = {
        1.0: {"penalty": 0.0, "penalty_images": [], "D14a_row": 2.01, "D14a_images": 2.0},
        0.9: {
            "penalty": round(P, 2),
            "penalty_images": ["D14a.1"],
            "D14a_row": 3759.12,
            "D14a_images": 10.68,
        },
    }
    assert m.rule_p(ok) == ("confirmed", [])
    # a penalty that is not a multiple of the unit
    bad = {**ok, 0.8: {**ok[0.9], "penalty": 3000.0}}
    assert m.rule_p(bad)[0] == "rejected"
    # the D14a excess without a D14a penalty warning
    bad = {**ok, 0.8: {**ok[0.9], "penalty_images": ["D52b.1"]}}
    out, why = m.rule_p(bad)
    assert out == "rejected" and why[0].startswith("f 0.80")
    # two D14a penalties need twice the excess
    two = {
        0.9: {
            "penalty": round(2 * P, 2),
            "penalty_images": ["D14a.1", "D14a.3"],
            "D14a_row": 10.68 + 2 * P,
            "D14a_images": 10.68,
        }
    }
    assert m.rule_p(two)[0] == "confirmed"


def test_posthoc_d14a(tmp_path, monkeypatch):
    monkeypatch.setattr(m, "F_GRID", (1.0,))
    d = tmp_path / "p1.00"
    d.mkdir()
    dec = -24.0
    dra = 1.0 / 3600.0 / np.cos(np.radians(dec))  # D14a.2 1" east of D14a.1 -> x = -1
    (d / "arcs.dat").write_text(
        f"#REFERENCE 0\nD14a.1 64.0 {dec} 0.3 0.3 0 2.09 25\n"
        f"D14a.2 {64.0 + dra} {dec} 0.3 0.3 0 2.09 25\n"
    )
    # chires prints -dx and dy; D14a.1 predicted at (-1, 0.5), D14a.2 at (-1, 0.5)
    (d / "chires.dat").write_text(
        "29 D14a.1 2.095   1   1.0 0 0 0 0.1 1.1  1.00  0.50  0\n"
        "29 D14a.2 2.095   1   1.0 0 0 0 0.1 0.5  0.00  0.50  0\n"
    )
    r = m.posthoc_d14a(tmp_path)["per_f"]["1.00"]
    assert r["pred_D14a.1_xy_arcsec"] == pytest.approx([-1.0, 0.5])
    assert r["pred_separation_arcsec"] == pytest.approx(0.0, abs=1e-3)
    # an image row without offsets gives no separation, not an error
    (tmp_path / "p1.00" / "chires.dat").write_text(
        "29 D14a.1 2.095   1   1.0 0 0 0 0.1 1.1  N/A  N/A  0\n"
    )
    r = m.posthoc_d14a(tmp_path)["per_f"]["1.00"]
    assert r["pred_D14a.1_xy_arcsec"] is None and r["pred_separation_arcsec"] is None


def test_build_f2b_frees_o1_to_o6_only():
    best = BEST2.replace("fini\n", "").replace(
        "potentiel 8515",
        "".join(f"potentiel O{i}\n\tend\n" for i in range(2, 7)) + "potentiel 8515",
    )
    best = best.replace("potentiel O1", "grille\n\tnombre      128\n\tend\npotentiel O1") + "fini\n"
    inp = INPUT2.replace(
        "finish\n",
        "".join(f"limit {i}\n    v_disp      1 1. 2. 0.1\n    end\n" for i in range(2, 7))
        + "finish\n",
    )
    vp = m.build_vprime(m.build_v0(best, inp), best, inp)
    out = m.build_f2b(vp, best, inp)
    lines = [ln.strip() for ln in out.splitlines()]
    assert "inverse 3 0.5 20" in lines and "inverse 0" not in lines and "chi2 1" in lines
    assert lines[lines.index("grille") + 1] == "nlens_opt   6"
    assert [ln for ln in lines if ln.startswith("limit")] == [f"limit {i}" for i in range(1, 7)]
    assert "60. 200." not in out and lines[-1] == "fini"
    with pytest.raises(ValueError):
        m.build_f2b(vp, BEST2, inp)  # fewer than six O potentials


def test_read_bayes_and_gate_v2():
    text = "#Nsample\n#ln(Lhood)\n#Chi2\n#O1 : v_disp\n1 -400 800.5 600\n2 -200 345.1 610\n3 -1\n"
    assert m.read_bayes(text) == (2, 345.1)
    assert m.read_bayes("#Nsample\n#Chi2\n") == (0, None)
    with pytest.raises(ValueError):
        m.read_bayes("#Nsample\n#Lhood\n1 2\n")
    assert m.gate_v2("finished", 344.42 * 1.01) == "pass"
    assert m.gate_v2("finished", 344.42 * 1.0101) == "fail"
    assert m.gate_v2("wall cap reached", 300.0) == "undetermined"
    assert m.gate_v2("finished", None) == "undetermined"


def test_posthoc_f2b(tmp_path, monkeypatch):
    import gzip

    monkeypatch.setattr(m, "OUT", tmp_path)
    (tmp_path / "f2b").mkdir()
    (tmp_path / "f2b" / "burnin.dat").write_text("#Nsample\n1 -900.5 3\n1 -100.25 4\n")
    with gzip.open(tmp_path / "lenstool_v0.log.gz", "wt") as g:
        g.write("INFO: chi2 356.674111  Lhood -302.923774\n")
    (tmp_path / "result_f2b.json").write_text('{"runtime_s": 600.0}')
    out = m.posthoc_f2b(tmp_path)
    assert out["n_burnin_lines"] == 2 and out["lnLhood_max_burnin"] == -100.25
    assert (
        out["lnLhood_stage_v_info_line"] == -302.923774 and out["seconds_per_burnin_line"] == 300.0
    )
    assert (tmp_path / "burnin_f2b.dat.gz").exists()


LIMITS = """limit 1
    v_disp      1 350. 1000. 0.1
    ellipticity 1 0.2 0.9 0.05
    end

limit 2
    cut_radius  1 5. 100. 0.1
    end
"""
POT = """potentiel O1
\tellipticite     0.800814
\tcore_radius         7.389947
\tcore_radius_kpc     39.462830
\tv_disp     574.453350
\tend
fini"""


def test_free_params_reads_limit_blocks():
    assert m.free_params(LIMITS, 2) == [
        ("O1", "v_disp", 350.0, 1000.0),
        ("O1", "ellipticity", 0.2, 0.9),
        ("O2", "cut_radius", 5.0, 100.0),
    ]
    with pytest.raises(ValueError):
        m.free_params(LIMITS.replace("cut_radius  1", "cut_radius  3"), 2)


def test_get_set_param_maps_keyword_and_drops_kpc():
    assert m.get_param(POT, "O1", "ellipticity") == 0.800814
    out = m.set_param(POT, "O1", "core_radius", 5.0)
    assert m.get_param(out, "O1", "core_radius") == 5.0
    assert "core_radius_kpc" not in out
    out = m.set_param(out, "O1", "ellipticity", 0.5)
    assert m.get_param(out, "O1", "ellipticity") == 0.5
    assert m.get_param(out, "O1", "v_disp") == 574.45335


def test_objective_feasibility():
    head = "chimul 344.42\n N ID z Narcs chip chix chiy chia rmss rmsi dx dy nwarn\n"
    ok = head + " 0 D1f.2 0.94 2 0.02 0 0 0 0.008 0.05 N/A N/A 0\n"
    assert m.objective(ok) == (344.42, True)
    lost = head + " 2 D1d.2 0.94 2 0.00 0 0 0 0.004 0.02 N/A N/A 1\n"
    assert m.objective(lost) == (m.F2L_INFEASIBLE, False)
    assert m.objective(ok + "chipenalty 3748.44\n") == (m.F2L_INFEASIBLE, False)


def test_perturbed_start_alternates_and_clips():
    u = m.perturbed_start(np.array([0.5, 0.5, 0.99, 0.0]))
    assert np.allclose(u, [0.53, 0.47, 0.999, 0.001])


def test_gate_v3():
    lim = m.V2_REF * 1.01
    assert m.gate_v3((lim - 0.1, True), lim - 0.1, [lim], False) == "uninformative"
    assert m.gate_v3((500.0, False), None, [m.F2L_INFEASIBLE], False) == "uninformative"
    assert m.gate_v3((500.0, True), lim, [500.0, lim], False) == "pass"
    falling = list(np.linspace(500, 400, 400))
    assert m.gate_v3((500.0, True), 400.0, falling, True) == "undetermined"
    assert m.gate_v3((500.0, True), 400.0, falling, False) == "fail"
    flat = [500.0] * 50 + [400.0] * 350
    assert m.gate_v3((500.0, True), 400.0, flat, True) == "fail"


def test_warn_families():
    t = (
        " 0 D1f.2 0.94 2 0.02 0 0 0 0.008 0.05 N/A N/A 0\n"
        " 2 D1d.2 0.94 2 0.00 0 0 0 0.004 0.02 N/A N/A 1\n"
        " 2 D1d.1 0.94 1 0.00 0 0 0 0.004 0.02 0.1 0.2 1\n"
    )
    assert m._warn_families(t) == {"D1d": 1}
