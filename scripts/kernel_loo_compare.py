"""Full-chain comparison of two ``arc_diff.py`` runs (refit vs shared stamp kernels).

PLAN ``results/transients/kernel_loo/PLAN.md`` §3 tests 2-4, tolerances fixed there:
- maps: max|ΔS| ≤ 0.01 per band over pixels finite in both (and max|Δd| reported);
- hits: same positions (0.01″), flags, signs, band sets, morphology verdicts and review list;
  every hit with 4.9 ≤ |S| ≤ 5.1 in any band, and every changed item, listed;
- injections (optional): each injected row's hit / sign / flag / morph / review outcome.

    python scripts/kernel_loo_compare.py --a A --b B [--a-inj DIR --b-inj DIR] --out cmp.json

Each run directory holds ``hits.ecsv``, ``maps.npz``, ``morph_all.json`` (and ``injected_pre.csv``).
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

TOL_S = 0.01
BORDER = (4.9, 5.1)


def _verdicts(d: Path) -> dict:
    return {
        t["id"]: t["verdict"] for t in json.loads((d / "morph_all.json").read_text())["targets"]
    }


def _review(h, v) -> set:
    return {
        str(r["hit_id"]) for r in h if not r["near_bright_core"] and v.get(r["hit_id"]) == "pass"
    }


def maps_diff(a: Path, b: Path) -> dict:
    za, zb = np.load(a / "maps.npz", allow_pickle=True), np.load(b / "maps.npz", allow_pickle=True)
    out = {}
    for k in sorted(za.keys()):
        if not (k.endswith("_s") or k.endswith("_d")):
            continue
        x, y = za[k], zb[k]
        ok = np.isfinite(x) & np.isfinite(y)
        out[k] = {
            "max_abs_diff": float(np.max(np.abs(x[ok] - y[ok]))) if ok.any() else None,
            "finite_mismatch_px": int((np.isfinite(x) != np.isfinite(y)).sum()),
        }
        if k.endswith("_s"):
            m = ok & (np.abs(x) >= BORDER[0]) & (np.abs(x) <= BORDER[1])
            out[k]["border_px"] = int(m.sum())
            out[k]["threshold_crossings_px"] = int(
                ((np.abs(x) >= 5.0) != (np.abs(y) >= 5.0))[ok].sum()
            )
    return out


def hits_diff(a: Path, b: Path) -> dict:
    ha, hb = Table.read(a / "hits.ecsv"), Table.read(b / "hits.ecsv")
    va, vb = _verdicts(a), _verdicts(b)

    def key(r):
        return (round(float(r["ra"]), 6), round(float(r["dec"]), 6))

    ma, mb = {key(r): r for r in ha}, {key(r): r for r in hb}
    only_a = [dict(ra=k[0], dec=k[1], id=str(ma[k]["hit_id"])) for k in ma if k not in mb]
    only_b = [dict(ra=k[0], dec=k[1], id=str(mb[k]["hit_id"])) for k in mb if k not in ma]
    changed, border = [], []
    for k in set(ma) & set(mb):
        ra_, rb_ = ma[k], mb[k]
        fields = {
            "near_bright_core": (bool(ra_["near_bright_core"]), bool(rb_["near_bright_core"])),
            "sign": (str(ra_["sign"]), str(rb_["sign"])),
            "bands": (str(ra_["bands"]), str(rb_["bands"])),
            "morph": (va.get(ra_["hit_id"]), vb.get(rb_["hit_id"])),
        }
        diff = {f: v for f, v in fields.items() if v[0] != v[1]}
        if diff:
            changed.append({"id": [str(ra_["hit_id"]), str(rb_["hit_id"])], **diff})
        s = float(ra_["max_abs_sigma"])
        if BORDER[0] <= s <= BORDER[1]:
            border.append(
                {"id": str(ra_["hit_id"]), "max_abs_sigma": [s, float(rb_["max_abs_sigma"])]}
            )
    ds = [
        abs(float(ma[k]["max_abs_sigma"]) - float(mb[k]["max_abs_sigma"]))
        for k in set(ma) & set(mb)
    ]
    return {
        "hits": [len(ha), len(hb)],
        "only_a": only_a,
        "only_b": only_b,
        "changed": changed,
        "border_hits_a": border,
        "max_abs_d_maxsigma": max(ds) if ds else None,
        "review": [sorted(_review(ha, va)), sorted(_review(hb, vb))],
    }


def inj_diff(a: Path, b: Path) -> dict:
    va, vb = _verdicts(a), _verdicts(b)
    ra = list(csv.DictReader(open(a / "injected_pre.csv")))
    rb = list(csv.DictReader(open(b / "injected_pre.csv")))

    def stage(r, v):
        hit = r["hit"] or None
        flag = r["near_bright_core"] == "True"
        morph = v.get(hit) if hit else None
        return {
            "hit": bool(hit),
            "sign": r["hit_sign"],
            "flag": flag,
            "morph": morph,
            "review": bool(hit) and not flag and morph == "pass",
        }

    changed = []
    for x, y in zip(ra, rb, strict=True):
        sa, sb = stage(x, va), stage(y, vb)
        if sa != sb:
            changed.append({"name": x["name"], "cls": x["cls"], "mag": x["mag"], "a": sa, "b": sb})
    return {"n": len(ra), "changed": changed}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--a", type=Path, required=True)
    p.add_argument("--b", type=Path, required=True)
    p.add_argument("--a-inj", type=Path)
    p.add_argument("--b-inj", type=Path)
    p.add_argument("--out", type=Path, required=True)
    x = p.parse_args(argv)
    out = {
        "label": "derived; a = refit, b = shared (PLAN §3 tests 2-4)",
        "tolerance_max_abs_dS": TOL_S,
        "maps": maps_diff(x.a, x.b),
        "hits": hits_diff(x.a, x.b),
    }
    if x.a_inj and x.b_inj:
        out["inj_maps"] = maps_diff(x.a_inj, x.b_inj)
        out["inj_hits"] = hits_diff(x.a_inj, x.b_inj)
        out["injections"] = inj_diff(x.a_inj, x.b_inj)
    s_ok = all(
        v["max_abs_diff"] is not None
        and v["max_abs_diff"] <= TOL_S
        and v["finite_mismatch_px"] == 0
        for m in ("maps", "inj_maps")
        if m in out
        for k, v in out[m].items()
        if k.endswith("_s")
    )
    h = [out["hits"]] + ([out["inj_hits"]] if "inj_hits" in out else [])
    h_ok = all(
        not (z["only_a"] or z["only_b"] or z["changed"]) and z["review"][0] == z["review"][1]
        for z in h
    )
    i_ok = not out.get("injections", {}).get("changed")
    out["pass"] = {"maps": s_ok, "hits": h_ok, "injections": i_ok, "all": s_ok and h_ok and i_ok}
    x.out.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out["pass"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
