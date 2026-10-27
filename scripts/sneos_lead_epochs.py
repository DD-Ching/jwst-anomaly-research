"""SN Eos leading pair: which public epochs could test it.

Plan: results/transients/sneos_lead/PLAN.md. Rule T per epoch x filter x image: band red of Lya at
z 5.133, 101.1-equivalent date inside the window where 101.1 is known bright (with a +-20 %
Delta-t allowance), and a 101.1 measurement in the same filter.
Metadata only: one MAST query per predicted position.
"""

from __future__ import annotations

import argparse
import json
from datetime import date, timedelta
from pathlib import Path

Z_SN = 5.133
LYA_A = 1216.0
DT_ALLOWANCE = 0.20  # ASSUMPTION (plan rule T2)
WINDOW = (date(2024, 3, 11), date(2025, 9, 1))  # 101.1 known bright (sneos_hst, sneos_ratio)
WINDOW_EXT = (date(2024, 3, 11), date(2026, 9, 28))  # sensitivity variant, not deciding
REF_FILTERS = frozenset(
    {"F110W", "F277W", "F356W", "F444W"}
)  # 101.1 measured (sneos_hst, sneos_ratio)

# Approximate pivot wavelengths in Angstrom (STScI instrument handbooks); only the comparison with
# 1216 * (1 + z) = 7457 A is used, so values need only be good to a few per cent.
PIVOT_A = {
    "F225W": 2359,
    "F275W": 2704,
    "F336W": 3355,
    "F390W": 3921,
    "F435W": 4329,
    "F475W": 4747,
    "FR505N": 5050,
    "F606W": 5887,
    "F625W": 6311,
    "F070W": 7039,
    "F775W": 7693,
    "F814W": 8045,
    "F090W": 9022,
    "F850LP": 9033,
    "F105W": 10552,
    "F110W": 11534,
    "F115W": 11543,
    "F125W": 12486,
    "F140W": 13923,
    "F150W": 15007,
    "F160W": 15369,
    "F150W2": 16592,
    "F200W": 19886,
    "F210M": 20954,
    "F277W": 27617,
    "F300M": 29891,
    "F356W": 35684,
    "F410M": 40822,
    "F444W": 44043,
    "F560W": 56350,
    "F770W": 76393,
    "F1000W": 99531,
    "F2100W": 207950,
}


def rule_band(filt: str, z: float = Z_SN) -> bool:
    """T1: pivot red of Lyman alpha at the SN redshift. Unknown filters raise KeyError."""
    return PIVOT_A[filt] > LYA_A * (1.0 + z)


def parse_filters(field: str) -> list[str]:
    """Science filters in a MAST `filters` value (CLEAR and detection dropped)."""
    toks = [t.strip().upper() for t in field.split(";")]
    return [t for t in toks if t and t not in ("DETECTION", "--", "CLEAR")]


def rule_phase(
    t: date, dt_days: float, window: tuple[date, date] = WINDOW, allowance: float = DT_ALLOWANCE
) -> bool:
    """T2: [t + (1 - a)|dt|, t + (1 + a)|dt|] overlaps the window (dt < 0: image leads 101.1)."""
    lo = t + timedelta(days=(1.0 - allowance) * abs(dt_days))
    hi = t + timedelta(days=(1.0 + allowance) * abs(dt_days))
    return lo <= window[1] and hi >= window[0]


def rule_reference(filt: str) -> bool:
    """T3: 101.1 has a measurement in the same filter."""
    return filt in REF_FILTERS


def classify(t: date, filt: str, dt_days: float) -> dict:
    """Rule T for one row. A filter with no tabulated pivot fails T1 (it cannot pass T3)."""
    t1 = rule_band(filt) if filt in PIVOT_A else False
    t2 = rule_phase(t, dt_days)
    t3 = rule_reference(filt)
    return {
        "band": t1,
        "phase": t2,
        "reference": t3,
        "pass": t1 and t2 and t3,
        "band_known": filt in PIVOT_A,
        "phase_ext": rule_phase(t, dt_days, WINDOW_EXT),
    }


def outcome(rows: list[dict]) -> str:
    return "P" if any(r["pass"] for r in rows) else "U"


def _query(ra: float, dec: float) -> list[tuple[date, str, str, str]]:
    import numpy as np
    from astropy.time import Time
    from astroquery.mast import Observations

    t = Observations.query_criteria(
        coordinates=f"{ra} {dec}",
        radius="0.0003 deg",
        dataproduct_type="image",
        intentType="science",
    )
    out = set()
    for r in t:
        if str(r["obs_collection"]) not in ("HST", "JWST") or not np.isfinite(float(r["t_min"])):
            continue
        filts = parse_filters(str(r["filters"]))
        if not filts:
            continue
        d = date.fromisoformat(Time(float(r["t_min"]), format="mjd").iso[:10])
        for filt in filts:
            out.add((d, str(r["obs_collection"]), str(r["instrument_name"]), filt))
    return sorted(out)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--summary", default="results/transients/sneos/summary.json")
    ap.add_argument("--out", default="results/transients/sneos_lead/epochs.json")
    args = ap.parse_args()
    images = json.loads(Path(args.summary).read_text())["images"]
    lead = [im for im in images if not im["catalogued"] and -2000 < im["dt_days_vs_101.1"] < 0]
    rows = []
    for name, im in zip(
        ("lead_a", "lead_b"), sorted(lead, key=lambda i: i["dt_days_vs_101.1"]), strict=True
    ):
        for d, coll, inst, filt in _query(im["ra"], im["dec"]):
            rows.append(
                {
                    "image": name,
                    "dt_days": im["dt_days_vs_101.1"],
                    "mu": im["mu_ours"],
                    "date": d.isoformat(),
                    "collection": coll,
                    "instrument": inst,
                    "filter": filt,
                    **classify(d, filt, im["dt_days_vs_101.1"]),
                }
            )
    res = {
        "plan": "results/transients/sneos_lead/PLAN.md",
        "window": [w.isoformat() for w in WINDOW],
        "window_ext": [w.isoformat() for w in WINDOW_EXT],
        "dt_allowance": DT_ALLOWANCE,
        "lya_obs_A": LYA_A * (1 + Z_SN),
        "n_rows": len(rows),
        "n_pass": sum(r["pass"] for r in rows),
        "n_unknown_filter": sum(not r["band_known"] for r in rows),
        "n_pass_ext": sum(r["band"] and r["reference"] and r["phase_ext"] for r in rows),
        "outcome": outcome(rows),
        "rows": rows,
    }
    Path(args.out).write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k != "rows"}, indent=1))


if __name__ == "__main__":
    main()
