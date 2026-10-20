"""PLAN.md steps 1 and 3: G0 (finite k = 0 fluxes unchanged under k = 1; e0 against the committed
persistence_run table, e2 against this run's k = 0 table) and the `on_sky` fraction of NEXUS e2
monitored pairs (mosaic S/N >= 10) by exposure class, k = 0 beside k = 1,
with the frozen rules R1-R3.
Usage (repo root): score.py  -> results/dimming/nan_k1/score.json"""

import json
from pathlib import Path

import numpy as np
from astropy.table import Table

OLD, NEW = Path("results/dimming/persistence_run/nexus"), Path("results/dimming/nan_k1/nexus")
EPS_K1, EPS_K0, EPS_GT4 = 0.584, 0.412, 0.961  # PLAN.md


def wilson(k: int, n: int, z: float = 1.96) -> list[float] | None:
    """Same formula as expcount_eff/score.py."""
    if n == 0:
        return None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(float(c - h), 3), round(float(c + h), 3)]


def g0(exp_a: Table, ver_a: Table, exp_b: Table, ver_b: Table) -> dict:
    """Exposures with finite k = 0 flux keep the identical flux; pairs with no NaN k = 0
    exposure keep their verdict; count what k = 1 adds."""
    key = lambda t: [f"{u}/{f}" for u, f in zip(t["uid"], t["file"], strict=True)]  # noqa: E731
    ib = {k: i for i, k in enumerate(key(exp_b))}
    fa = np.asarray(exp_a["flux"], float)
    fin = np.isfinite(fa)
    idx = np.array([ib.get(k, -1) for k in key(exp_a)])
    fb = np.where(idx >= 0, np.asarray(exp_b["flux"], float)[np.maximum(idx, 0)], np.nan)
    vb_of = dict(zip(ver_b["uid"], ver_b["verdict"], strict=True))
    all_finite = {u for u in ver_a["uid"] if np.isfinite(fa[exp_a["uid"] == u]).all()}
    res = {
        "rows_k0": len(exp_a),
        "rows_k1": len(exp_b),
        "missing_in_k1": int((idx < 0).sum()),
        "finite_k0": int(fin.sum()),
        "finite_k0_changed": int((fin & ~(fb == fa)).sum()),
        "nan_k0_finite_k1": int((~fin & np.isfinite(fb)).sum()),
        "pairs_without_nan_k0": len(all_finite),
        "pairs_without_nan_k0_verdict_changed": sum(
            str(vb_of.get(u)) != str(r)
            for u, r in zip(ver_a["uid"], ver_a["verdict"], strict=True)
            if u in all_finite
        ),
    }
    res["pass"] = (
        res["finite_k0_changed"] == 0
        and res["missing_in_k1"] == 0
        and res["pairs_without_nan_k0_verdict_changed"] == 0
    )
    return res


def table(v: Table) -> dict:
    m = v[(v["kind"] == "monitored") & (v["mosaic_snr"] >= 10)]
    out = {}
    for cls, c in (("le4", m["n_exposures"] <= 4), ("gt4", m["n_exposures"] > 4)):
        q = m[c]
        k = int((q["verdict"] == "on_sky").sum())
        out[cls] = {
            "on_sky": k,
            "n": len(q),
            "wilson95": wilson(k, len(q)),
            "persistence": int((q["verdict"] == "persistence").sum()),
            "verdicts": {
                str(a): int(b)
                for a, b in zip(*np.unique(q["verdict"], return_counts=True), strict=True)
            },
        }
    n_mon = len(m)
    n_pers = int((m["verdict"] == "persistence").sum())
    out["monitored"] = {"n": n_mon, "persistence": n_pers, "wilson95": wilson(n_pers, n_mon)}
    return out


def read(d: Path, epoch: int) -> tuple[Table, Table] | None:
    if not (d / "verdicts.ecsv").exists():
        return None
    return Table.read(d / f"exposures_e{epoch}.ecsv"), Table.read(d / "verdicts.ecsv")


def main() -> dict:
    res: dict = {"provenance": "derived", "eps_k1_proposed": EPS_K1}
    e0 = read(NEW / "e0", 0)
    if e0 is not None:
        old_v = Table.read(OLD / "verdicts.ecsv")
        res["G0_e0"] = g0(Table.read(OLD / "exposures_e0.ecsv"), old_v[old_v["epoch"] == 0], *e0)
    k0, k1 = read(NEW / "e2_k0", 2), read(NEW / "e2_k1", 2)
    if k0 is not None and k1 is not None:
        res["G0_e2"] = g0(*k0, *k1)
    if k0 is not None:
        res["e2_k0"] = table(k0[1])
        w = res["e2_k0"]["gt4"]["wilson95"]
        res["gt4_rule"] = "not tested" if w is None else ("holds" if w[1] >= EPS_GT4 else "fails")
    if k1 is None:
        res["decision"] = "R3 (e2 k = 1 unfinished)"
        return res
    t = res["e2_k1"] = table(k1[1])
    le4, mon = t["le4"], t["monitored"]
    g0_ok = res.get("G0_e0", {}).get("pass", False)
    if le4["n"] < 10:
        res["decision"] = f"R3 (n(<= 4) = {le4['n']} < 10)"
    elif (
        g0_ok
        and le4["wilson95"][1] >= EPS_K1
        and le4["wilson95"][0] > EPS_K0
        and (mon["wilson95"] is None or mon["wilson95"][0] <= 0.02)
    ):
        res["decision"] = "R1 (adopt k = 1 in the next PR)"
    else:
        res["decision"] = "R2 (not adopted)"
    return res


if __name__ == "__main__":
    out = main()
    Path("results/dimming/nan_k1/score.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
