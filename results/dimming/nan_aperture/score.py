"""PLAN.md steps 1 and 4: G0 (finite k = 0 fluxes unchanged under k = 2) and the `on_sky` fraction
of NEXUS monitored pairs (mosaic S/N >= 10) by exposure class and split, k = 0 (committed
persistence_run) beside k = 2 (this run), over the epochs this run finished.
Usage (repo root): score.py  -> results/dimming/nan_aperture/score.json"""

import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack

sys.path.insert(0, "scripts")
import persistence_check as pc  # noqa: E402

OLD, NEW = Path("results/dimming/persistence_run/nexus"), Path("results/dimming/nan_aperture/nexus")
TRAIN, HELD = {0, 1, 3, 4}, {5, 6, 7}  # as expcount_eff


def wilson(k: int, n: int, z: float = 1.96) -> list[float] | None:
    """Same formula as expcount_eff/score.py."""
    if n == 0:
        return None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(float(c - h), 3), round(float(c + h), 3)]


def g0(epoch: int) -> dict:
    """Exposures with finite k = 0 flux keep the identical flux and S/N; count what k = 2 adds."""
    a, b = (
        Table.read(OLD / f"exposures_e{epoch}.ecsv"),
        Table.read(NEW / f"e{epoch}" / f"exposures_e{epoch}.ecsv"),
    )
    key = lambda t: [f"{u}/{f}" for u, f in zip(t["uid"], t["file"], strict=True)]  # noqa: E731
    ib = {k: i for i, k in enumerate(key(b))}
    fa = np.asarray(a["flux"], float)
    fin = np.isfinite(fa)
    idx = np.array([ib.get(k, -1) for k in key(a)])
    fb = np.where(idx >= 0, np.asarray(b["flux"], float)[np.maximum(idx, 0)], np.nan)
    nan_b = np.asarray(b["n_ap_nan"])
    return {
        "rows_k0": len(a),
        "rows_k2": len(b),
        "missing_in_k2": int((idx < 0).sum()),
        "finite_k0": int(fin.sum()),
        "finite_k0_changed": int((fin & ~(fb == fa)).sum()),
        "nan_k0_finite_k2": int((~fin & np.isfinite(fb)).sum()),
        "nan_k2_with_1_2_nan_px": int(
            (~np.isfinite(np.asarray(b["flux"], float)) & (nan_b >= 1) & (nan_b <= 2)).sum()
        ),
    }


def table(v: Table, epochs: list[int]) -> dict:
    m = v[(v["kind"] == "monitored") & (v["mosaic_snr"] >= 10)]
    out = {}
    for split, ep in (("train", TRAIN), ("held_out", HELD)):
        s = np.isin(m["epoch"], sorted(ep & set(epochs)))
        out[split] = {}
        for cls, c in (("le4", m["n_exposures"] <= 4), ("gt4", m["n_exposures"] > 4)):
            q = m[s & c]
            k = int((q["verdict"] == "on_sky").sum())
            out[split][cls] = {
                "on_sky": k,
                "n": len(q),
                "wilson95": wilson(k, len(q)),
                "verdicts": {
                    str(a): int(b)
                    for a, b in zip(*np.unique(q["verdict"], return_counts=True), strict=True)
                },
            }
    return out


def k1_approx(epochs: list[int]) -> Table:
    """POST HOC: k = 1 verdicts approximated from the k = 2 tables by dropping exposures whose own
    aperture had 2 NaN pixels (earlier-exposure fills are not undone; no k = 2 row is suspect)."""
    rows = []
    for k in epochs:
        v = Table.read(NEW / f"e{k}" / "verdicts.ecsv")
        e = Table.read(NEW / f"e{k}" / f"exposures_e{k}.ecsv")
        snr = np.where(np.asarray(e["n_ap_nan"]) >= 2, np.nan, np.asarray(e["snr"], float))
        for r in v:
            s = e["uid"] == r["uid"]
            r = dict(zip(r.colnames, r, strict=True))
            r["verdict"] = pc.classify(snr[s], e["suspect"][s], e["n_prior"][s])
            rows.append(r)
    return Table(rows=rows)


def main() -> dict:
    epochs = sorted(int(p.name[1:]) for p in NEW.glob("e*") if (p / "verdicts.ecsv").exists())
    new = vstack(
        [Table.read(NEW / f"e{k}" / "verdicts.ecsv") for k in epochs], metadata_conflicts="silent"
    )
    old = Table.read(OLD / "verdicts.ecsv")
    old = old[np.isin(old["epoch"], epochs)]
    gate = {f"e{k}": g0(k) for k in epochs}
    res = {
        "provenance": "derived",
        "epochs": epochs,
        "G0": gate,
        "G0_pass": all(
            g["finite_k0_changed"] == 0 and g["missing_in_k2"] == 0 for g in gate.values()
        ),
        "k0": table(old, epochs),
        "k2_post_hoc": table(new, epochs),  # G1 failed (g1.json): descriptive only, not adopted
        "k1_post_hoc_approx": table(k1_approx(epochs), epochs),
        "persistence_k2": int((new["verdict"] == "persistence").sum()),
    }
    Path("results/dimming/nan_aperture/score.json").write_text(json.dumps(res, indent=1) + "\n")
    return res


if __name__ == "__main__":
    print(json.dumps(main(), indent=1))
