"""PLAN.md: one cause per measurable held-out <= 4-exposure NEXUS `inconclusive` pair, from the
committed `persistence_run/nexus` tables (no pixels read). Selection as undetected_stamps/score.py.
Usage (repo root): split.py  -> results/dimming/inconclusive_split/split.json"""

import json
from pathlib import Path

import numpy as np
from astropy.table import Table

RUN = Path("results/dimming/persistence_run/nexus")
HELD = (5, 6, 7)  # as expcount_eff / undetected_stamps
TRAIN = (0, 1, 3, 4)
MIN_SNR = 5.0  # D-039 per-exposure detection rule (ASSUMPTION), as persistence_check.MIN_SNR


SAME_EXPOSURE_S = 60.0  # SW and LW halves of one exposure start < 1 s apart (ASSUMPTION: 60 s)


def n_times(mjd: np.ndarray) -> int:
    """Distinct exposure start times; SW/LW detectors of one exposure count once."""
    t = np.sort(mjd) * 86400.0
    return int(len(t) and 1 + (np.diff(t) > SAME_EXPOSURE_S).sum())


def exposure_flags(e: Table) -> dict:
    """Per-exposure state in the order persistence_check.classify uses it."""
    snr = np.asarray(e["snr"], float)
    det = np.isfinite(snr) & (snr >= MIN_SNR)
    sus = np.asarray(e["suspect"], bool)
    npr = np.asarray(e["n_prior"], int)
    nfail = np.asarray(e["n_prior_failed"], int)
    return {
        "finite": int(np.isfinite(np.asarray(e["flux"], float)).sum()),
        "det": int(det.sum()),
        "clean": int((det & ~sus & (npr > 0)).sum()),
        "det_no_prior_first": int((det & ~sus & (npr == 0) & (nfail == 0)).sum()),
        "det_no_prior_failed": int((det & ~sus & (npr == 0) & (nfail > 0)).sum()),
        "det_suspect": int((det & sus).sum()),
        "finite_below_snr": int((np.isfinite(snr) & (snr < MIN_SNR)).sum()),
        "n_distinct_mjd": n_times(np.asarray(e["mjd"], float)[np.isfinite(snr)]),
        "snr": [round(float(s), 2) for s in snr],
    }


def cause(f: dict) -> str:
    """PLAN.md order: C1 if no-prior detections counted as clean give >= 2; C2 if a finite exposure
    is below S/N 5 (and C1 does not hold); C3 otherwise."""
    if f["clean"] + f["det_no_prior_first"] + f["det_no_prior_failed"] >= 2:
        return "C1_n_prior"
    if f["finite_below_snr"] > 0:
        return "C2_low_snr"
    return "C3_other"


def split(epochs: tuple[int, ...]) -> tuple[list[dict], dict]:
    v = Table.read(RUN / "verdicts.ecsv")
    m = v[
        (v["kind"] == "monitored")
        & (v["mosaic_snr"] >= 10)
        & (v["n_exposures"] <= 4)
        & np.isin(v["epoch"], epochs)
    ]
    rows, counts = [], {}
    for k in epochs:
        e = Table.read(RUN / f"exposures_e{k}.ecsv")
        for r in m[(m["epoch"] == k) & (m["verdict"] == "inconclusive")]:
            f = exposure_flags(e[e["uid"] == r["uid"]])
            if f["finite"] < 2:  # not measurable (undetected_stamps definition)
                continue
            c = cause(f)
            counts[c] = counts.get(c, 0) + 1
            rows.append(
                {
                    "uid": str(r["uid"]),
                    "epoch": k,
                    "n_exposures": int(r["n_exposures"]),
                    "mosaic_snr": round(float(r["mosaic_snr"]), 1),
                    "cause": c,
                    **f,
                }
            )
    return rows, counts


def main() -> dict:
    rows, counts = split(HELD)
    # C1 sub-split: first-of-visit only vs a failed prior read
    sub = {"first_of_visit": 0, "prior_failed": 0}
    for r in rows:
        if r["cause"] == "C1_n_prior":
            sub["prior_failed" if r["det_no_prior_failed"] else "first_of_visit"] += 1
    t_rows, t_counts = split(TRAIN)  # post hoc: train epochs were used to set the D-039 rules
    out = {
        "n_pairs": len(rows),
        "cause_counts": counts,
        "c1_subsplit": sub,
        "single_time_pairs": sum(r["n_distinct_mjd"] == 1 for r in rows),
        "pairs": rows,
        "posthoc_train": {
            "n_pairs": len(t_rows),
            "cause_counts": t_counts,
            "single_time_pairs": sum(r["n_distinct_mjd"] == 1 for r in t_rows),
            "pairs": t_rows,
        },
        "meta": {"provenance": "derived", "min_snr": MIN_SNR, "epochs": list(HELD)},
    }
    Path("results/dimming/inconclusive_split/split.json").write_text(
        json.dumps(out, indent=1) + "\n"
    )
    return out


if __name__ == "__main__":
    o = main()
    print(
        json.dumps(
            {k: o[k] for k in ("n_pairs", "cause_counts", "c1_subsplit", "single_time_pairs")},
            indent=1,
        )
    )
    print("train", {k: v for k, v in o["posthoc_train"].items() if k != "pairs"})
    for r in o["pairs"] + o["posthoc_train"]["pairs"]:
        print(
            r["uid"],
            r["epoch"],
            r["n_exposures"],
            r["mosaic_snr"],
            r["cause"],
            r["snr"],
            "clean",
            r["clean"],
            "first",
            r["det_no_prior_first"],
            "fail",
            r["det_no_prior_failed"],
            "sus",
            r["det_suspect"],
            "t",
            r["n_distinct_mjd"],
        )
