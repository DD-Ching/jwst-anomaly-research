"""No-signalling check on the NIST 2015 loophole-free Bell test (goal round 10).

Does one party's spacelike-separated setting choice change the other party's local click
statistics? Pre-registered in docs/goal_round10/preregistration.md. Input: NIST belltestdata
`Peter's Data Calculations/Diagnostic Statistics/*.xlsx` (not committed; sha256 prefixes are
written to the results). Each pulse-grouping tab's first row holds 16 whole-run counts
N(oA, oB | x, y) in the order ++, +0, 0+, 00 for the setting pairs ab, ab', a'b, a'b'.

The tests are written independently of the spreadsheet formulas (binomial G-test and Fisher exact
test). NIST's own pooled two-proportion z (cells in the same tab) is read back only as a
reproduction check.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy import stats

RUNS = {
    "19_45": "diag_19-45.xlsx",
    "02_54": "diag_02-54.xlsx",
    "03_43": "diag_03-43.xlsx",
    "XOR_1": "diag_Classical XOR 1.xlsx",
    "XOR_2": "diag_Classical XOR 2.xlsx",
    "XOR_3": "diag_Diagnostics - Classical XOR 3.xlsx",
}
GROUPINGS = ("5", "456", "34567", "2345678", "123456789")
PRIMARY = "2345678"  # largest grouping reported as spacelike in the paper (pre-registered)
PAIRS = ("ab", "ab'", "a'b", "a'b'")
TRIALS_PER_S = 99_100  # Shalm et al. 2015, trial rate


def read_tab(ws) -> dict:
    rows = [list(r) for r in ws.iter_rows(min_row=1, max_row=30, values_only=True)]
    counts = [int(v) for v in rows[0][:16]]
    n = {}
    for k, pair in enumerate(PAIRS):
        pp, p0, zp, zz = counts[4 * k : 4 * k + 4]
        n[pair] = {"++": pp, "+0": p0, "0+": zp, "00": zz, "N": pp + p0 + zp + zz}
    # the sheet's own trial totals per setting pair (a row labelled ab, ab', a'b, a'b'), if present
    totals = {}
    for r in rows:
        cells = [c for c in r if c is not None]
        if cells and cells[0] == "ab" and len(cells) >= 8:
            totals = {str(cells[i]): int(cells[i + 1]) for i in range(0, 8, 2)}
            break
    # NIST's own pooled z values (reproduction check only): the row after the first 'p1-hat' row
    # holds Alice's two z values (columns E, M); Bob's are in the first row whose first cell is 'z'.
    nist = {}
    for i, r in enumerate(rows):
        if r[1] == "p1-hat" and "A|a" not in nist:
            nist["A|a"], nist["A|a'"] = float(rows[i + 1][4]), float(rows[i + 1][12])
        cells = [(j, c) for j, c in enumerate(r) if c is not None]
        if cells and cells[0][1] == "z" and "B|b" not in nist:
            vals = [c for _, c in cells if isinstance(c, (int, float))]
            nist["B|b"], nist["B|b'"] = float(vals[0]), float(vals[1])
    return {"n": n, "totals": totals, "nist_z": nist}


def comparisons(n: dict) -> dict:
    """Successes and trials for each no-signalling comparison: (k1, n1, k2, n2)."""

    def a_plus(p):
        return n[p]["++"] + n[p]["+0"]

    def b_plus(p):
        return n[p]["++"] + n[p]["0+"]

    return {
        "A|a": (a_plus("ab"), n["ab"]["N"], a_plus("ab'"), n["ab'"]["N"]),  # Bob's setting varies
        "A|a'": (a_plus("a'b"), n["a'b"]["N"], a_plus("a'b'"), n["a'b'"]["N"]),
        "B|b": (b_plus("ab"), n["ab"]["N"], b_plus("a'b"), n["a'b"]["N"]),  # Alice's setting varies
        "B|b'": (b_plus("ab'"), n["ab'"]["N"], b_plus("a'b'"), n["a'b'"]["N"]),
    }


def g_test(k1, n1, k2, n2) -> tuple[float, float]:
    """Binomial likelihood-ratio test of p1 = p2 (G statistic, chi^2_1)."""

    def ll(k, n, p):
        return k * math.log(p) + (n - k) * math.log1p(-p)

    p1, p2, p = k1 / n1, k2 / n2, (k1 + k2) / (n1 + n2)
    g = 2 * (ll(k1, n1, p1) + ll(k2, n2, p2) - ll(k1, n1, p) - ll(k2, n2, p))
    return g, float(stats.chi2.sf(g, 1))


def pooled_z(k1, n1, k2, n2) -> float:
    p1, p2, p = k1 / n1, k2 / n2, (k1 + k2) / (n1 + n2)
    return (p1 - p2) / math.sqrt(p * (1 - p) * (1 / n1 + 1 / n2))


def bern_mi(p0: float, p1: float) -> float:
    """Mutual information (bits): uniform binary input switching a Bernoulli output p0 <-> p1."""

    def h(p):
        return 0.0 if p in (0.0, 1.0) else -(p * math.log2(p) + (1 - p) * math.log2(1 - p))

    return h((p0 + p1) / 2) - (h(p0) + h(p1)) / 2


def analyse(data_dir: Path) -> dict:
    import openpyxl

    out = {"runs": {}, "primary": PRIMARY}
    prim_p = []
    for run, fname in RUNS.items():
        wb = openpyxl.load_workbook(data_dir / fname, data_only=True)
        out["runs"][run] = {}
        for g in GROUPINGS:
            tab = read_tab(wb[g])
            n = tab["n"]
            tm = all(tab["totals"][p] == n[p]["N"] for p in PAIRS) if tab["totals"] else None
            row = {"totals_match": tm, "tests": {}}
            for name, (k1, n1, k2, n2) in comparisons(n).items():
                gs, pg = g_test(k1, n1, k2, n2)
                _, pf = stats.fisher_exact([[k1, n1 - k1], [k2, n2 - k2]])
                p1, p2 = k1 / n1, k2 / n2
                se = math.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
                d = p1 - p2
                row["tests"][name] = {
                    "k1": k1,
                    "n1": n1,
                    "k2": k2,
                    "n2": n2,
                    "delta": d,
                    "se": se,
                    "delta_rel": d / ((k1 + k2) / (n1 + n2)),
                    "G": gs,
                    "p_G": pg,
                    "p_fisher": float(pf),
                    "z_pooled": pooled_z(k1, n1, k2, n2),
                    "z_nist": tab["nist_z"][name],
                }
                if g == PRIMARY:
                    prim_p.append(pg)
            out["runs"][run][g] = row
    m = len(prim_p)
    out["primary_family"] = {
        "m": m,
        "bonferroni_alpha": 0.05 / m,
        "min_p": min(prim_p),
        "n_below_bonferroni": sum(p < 0.05 / m for p in prim_p),
        "n_below_0.05": sum(p < 0.05 for p in prim_p),
        "fisher_combined_p": float(stats.combine_pvalues(prim_p, method="fisher").pvalue),
        "ks_uniform_p": float(stats.kstest(prim_p, "uniform").pvalue),
    }
    # inverse-variance combination of delta over the six independent runs (primary grouping)
    comb = {}
    for name in ("A|a", "A|a'", "B|b", "B|b'"):
        ds = np.array([out["runs"][r][PRIMARY]["tests"][name]["delta"] for r in RUNS])
        ses = np.array([out["runs"][r][PRIMARY]["tests"][name]["se"] for r in RUNS])
        ks = sum(out["runs"][r][PRIMARY]["tests"][name]["k1"] for r in RUNS)
        ns = sum(out["runs"][r][PRIMARY]["tests"][name]["n1"] for r in RUNS)
        w = 1 / ses**2
        d = float(np.sum(w * ds) / np.sum(w))
        se = float(1 / math.sqrt(np.sum(w)))
        pbar = ks / ns
        bound = abs(d) + 1.96 * se
        mi = bern_mi(pbar, pbar + bound)
        comb[name] = {
            "delta": d,
            "se": se,
            "z": d / se,
            "p_click": pbar,
            "abs_delta_95_upper": bound,
            "rel_delta_95_upper": bound / pbar,
            "mi_bits_per_trial_upper": mi,
            "mi_bits_per_s_upper": mi * TRIALS_PER_S,
            "trials_per_bit_lower": 1 / mi,
        }
    out["combined_primary"] = comb
    # post hoc (not pre-registered): heterogeneity of delta across runs (Cochran Q), and the
    # combination without the run holding the single family-wise hit
    het = {}
    for name in ("A|a", "A|a'", "B|b", "B|b'"):
        ds = np.array([out["runs"][r][PRIMARY]["tests"][name]["delta"] for r in RUNS])
        ses = np.array([out["runs"][r][PRIMARY]["tests"][name]["se"] for r in RUNS])
        w = 1 / ses**2
        d = np.sum(w * ds) / np.sum(w)
        q = float(np.sum(w * (ds - d) ** 2))
        het[name] = {"cochran_Q": q, "dof": len(ds) - 1, "p": float(stats.chi2.sf(q, len(ds) - 1))}
    keep = [r for r in RUNS if r != "XOR_2"]
    ds = np.array([out["runs"][r][PRIMARY]["tests"]["A|a"]["delta"] for r in keep])
    ses = np.array([out["runs"][r][PRIMARY]["tests"]["A|a"]["se"] for r in keep])
    w = 1 / ses**2
    d_wo, se_wo = float(np.sum(w * ds) / np.sum(w)), float(1 / math.sqrt(np.sum(w)))
    x2 = out["runs"]["XOR_2"][PRIMARY]["tests"]["A|a"]
    out["post_hoc"] = {
        "heterogeneity": het,
        "A|a_without_XOR_2": {"delta": d_wo, "se": se_wo, "z": d_wo / se_wo},
        "A|a_XOR_2_minus_rest_z": (x2["delta"] - d_wo) / math.sqrt(x2["se"] ** 2 + se_wo**2),
    }
    z_dev = [
        abs(t["z_pooled"] - t["z_nist"])
        for r in out["runs"].values()
        for row in r.values()
        for t in row["tests"].values()
    ]
    out["reproduction_max_abs_z_difference"] = max(z_dev)
    tms = [row["totals_match"] for r in out["runs"].values() for row in r.values()]
    out["totals_checked"] = sum(t is not None for t in tms)
    out["totals_all_match"] = all(t for t in tms if t is not None)
    return out


def file_hashes(data_dir: Path) -> dict:
    import hashlib

    return {f: hashlib.sha256((data_dir / f).read_bytes()).hexdigest()[:16] for f in RUNS.values()}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument(
        "data_dir", type=Path, help="directory with the six NIST diagnostic .xlsx files"
    )
    args = ap.parse_args(argv)
    out = analyse(args.data_dir)
    out["provenance"] = (
        "derived from NIST belltestdata (s3://nist-belltestdata/belldata/"
        "Peter's Data Calculations/Diagnostic Statistics/, accessed 2026-10-11; "
        "NIST public-service licence, attribution required)"
    )
    out["sha256_prefix"] = file_hashes(args.data_dir)
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
