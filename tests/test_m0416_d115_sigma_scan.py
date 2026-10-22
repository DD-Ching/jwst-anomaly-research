import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import m0416_d115_sigma_scan as m  # noqa: E402


def _fits(chi2, shared=None, flag=None):
    """chi2[f][fam] -> fit records; shared / flag override per (f, fam)."""
    shared, flag = shared or {}, flag or {}
    out = []
    for f, fams in chi2.items():
        for fam, c in fams.items():
            r = {"f": f, "family": fam, "chi2_fit": c, "chi2_validate_source": c,
                 "shared": shared.get((f, fam), False)}  # fmt: skip
            if fam == "D115b":
                r["flag_image"] = flag.get(f, f == 1.0)
            out.append(r)
    return out


def _chi2(delta):
    return {f: {"D115a": 1.0, "D115b": 0.5, "X": 3.0 + delta.get(f, 0.0)} for f in m.FACTORS}


def test_verdict_picks_largest_factor_without_flag():
    fits = _fits(_chi2({0.85: 2.0, 0.8: 9.0}), flag={0.85: True})
    scan = m.summarize(fits)
    assert scan["0.85"]["flag_image"] and not scan["0.80"]["flag_image"]
    assert m.verdict(scan) == ("inconclusive", "0.80")  # f* 0.80 costs 9 > 6
    scan = m.summarize(_fits(_chi2({0.85: 5.9})))
    assert m.verdict(scan) == ("ordinary", "0.85")


def test_lost_match_and_no_fstar():
    scan = m.summarize(_fits(_chi2({}), shared={(0.85, "X"): True}))
    assert scan["0.85"]["lost_match"] == ["X"] and scan["0.85"]["n_families_compared"] == 2
    assert m.verdict(scan) == ("inconclusive", "0.85")
    flags = {f: True for f in m.FACTORS}
    assert m.verdict(m.summarize(_fits(_chi2({}), flag=flags))) == ("inconclusive", None)
    shared = {(f, "D115a"): True for f in m.FACTORS if f < 1}
    assert m.verdict(m.summarize(_fits(_chi2({}), shared=shared))) == ("inconclusive", None)


def test_gates():
    c = _chi2({})
    c[0.70].update({"D115a": 2.78, "D115b": 0.53})
    g = m.gates(_fits(c))
    assert g["G1"] and g["G2"]
    c[0.70]["D115a"] = 2.9
    fits = _fits(c)
    fits[0]["chi2_validate_source"] = 0.5  # f = 1 D115a fit worse than its start
    g = m.gates(fits)
    assert not g["G1"] and g["G2_failures"] == ["D115a"]
