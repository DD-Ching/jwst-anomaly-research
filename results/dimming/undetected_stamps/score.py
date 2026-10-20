"""PLAN.md step 2: `on_sky` fraction over all and over measurable NEXUS monitored pairs (>= 2
covering exposures with finite aperture flux), by exposure class and split; R1 / R2 via decide().
Usage (repo root): score.py  -> results/dimming/undetected_stamps/score.json"""

import json
from pathlib import Path

import numpy as np
from astropy.table import Table

RUN = Path("results/dimming/persistence_run/nexus")
TRAIN, HELD = {0, 1, 3, 4}, {5, 6, 7}  # as expcount_eff
R1_LOWER = 0.80  # PLAN.md R1 (ASSUMPTION; the D-039 rule-1 level)


def wilson_raw(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Same formula as expcount_eff/score.py (that module runs on import, so it is not imported)."""
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return float(c - h), float(c + h)


def n_finite(per_exp: Table, uid: str) -> int:
    return int(np.isfinite(np.asarray(per_exp["flux"][per_exp["uid"] == uid], float)).sum())


def decide(k_meas: int, n_meas: int) -> str:
    """PLAN.md: R1 if the held-out <= 4 measurable lower bound is >= 0.80, else R2."""
    if n_meas == 0:
        return "not tested"
    return "R1" if wilson_raw(k_meas, n_meas)[0] >= R1_LOWER else "R2"


def stamp_summary() -> dict:
    """Counts from stamps.py's reasons.ecsv (NaN-flux exposures of monitored pairs, e4-e7)."""
    t = Table.read(Path(__file__).parent / "out" / "reasons.ecsv")
    t = t[t["uid"] != "nexus-d30295"]
    out = {}
    for name, q in (("undetected_pairs", t[t["verdict"] == "undetected"]), ("all", t)):
        nd = q[q["reason"] == "aperture_no_data"]
        out[name] = {
            "pairs": len(set(q["uid"].astype(str) + "/" + q["epoch"].astype(str))),
            "exposures": len(q),
            "reason": {
                str(a): int(b)
                for a, b in zip(*np.unique(q["reason"], return_counts=True), strict=True)
            },
            "aperture_no_data_nan_pixels_le2": int((nd["n_ap_nan"] <= 2).sum()),
            "stamp_nan_frac_median": round(float(np.median(q["stamp_nan_frac"])), 4),
        }
    return out


def main() -> dict:
    v = Table.read(RUN / "verdicts.ecsv")
    m = v[(v["kind"] == "monitored") & (v["mosaic_snr"] >= 10)]
    nf = {}
    for k in sorted({int(e) for e in m["epoch"]}):
        e = Table.read(RUN / f"exposures_e{k}.ecsv")
        for u in m["uid"][m["epoch"] == k]:
            nf[(str(u), k)] = n_finite(e, u)
    meas = np.array([nf[(str(r["uid"]), int(r["epoch"]))] >= 2 for r in m])
    out = {}
    for split, ep in (("train", TRAIN), ("held_out", HELD)):
        s_ep = np.isin(m["epoch"], sorted(ep))
        out[split] = {}
        for cls, c in (("le4", m["n_exposures"] <= 4), ("gt4", m["n_exposures"] > 4)):
            row = {}
            for name, sel in (("all", s_ep & c), ("measurable", s_ep & c & meas)):
                q = m[sel]
                k, n = int((q["verdict"] == "on_sky").sum()), len(q)
                row[name] = {
                    "on_sky": k,
                    "n": n,
                    "wilson95": [round(b, 3) for b in wilson_raw(k, n)] if n else None,
                    "verdicts": {
                        str(a): int(b)
                        for a, b in zip(*np.unique(q["verdict"], return_counts=True), strict=True)
                    },
                }
            row["measurable_fraction"] = (
                round(row["measurable"]["n"] / row["all"]["n"], 3) if row["all"]["n"] else None
            )
            out[split][cls] = row
    h = out["held_out"]["le4"]["measurable"]
    out["decision"] = decide(h["on_sky"], h["n"])
    out["stamps"] = stamp_summary()
    out["meta"] = {"provenance": "derived", "mosaic_snr_min": 10, "r1_lower": R1_LOWER}
    Path("results/dimming/undetected_stamps/score.json").write_text(
        json.dumps(out, indent=1) + "\n"
    )
    return out


if __name__ == "__main__":
    print(json.dumps(main(), indent=1))
