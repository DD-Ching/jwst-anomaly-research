"""Reproduction gate for a regenerated per-field W3 dimming ``limits.ecsv``.

Frozen with results/dimming/macs0416_regen/PLAN.md before the run: per t_E at the given ρ, the
regenerated ``exposure_source_yr`` must lie within ``rel_tol`` (ASSUMPTION, 10 %) of the reference
values, and the field's forced stage must be calibrated.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from astropy.table import Table

# D-052 MACS0416 exposures (source yr, ρ = 0.1), as rounded in docs/exotic_limits.md and summed
# into the headline
MACS0416_D052 = {0.01: 1.2, 0.03: 4.7, 0.1: 9.9, 0.3: 18.7, 1.0: 53.7, 3.0: 62.6}


def repro_gate(
    lim: Table, reference: dict[float, float], rho: float = 0.1, rel_tol: float = 0.10
) -> dict:
    """Per-t_E comparison; ``passed`` needs every t_E within tolerance AND a calibrated field."""
    rows = []
    for t_e, ref in sorted(reference.items()):
        sel = lim[(abs(lim["t_e_yr"] - t_e) < 1e-9) & (abs(lim["rho"] - rho) < 1e-9)]
        new = float(sel["exposure_source_yr"][0]) if len(sel) else float("nan")
        rel = (new - ref) / ref
        rows.append(
            {
                "t_e_yr": t_e,
                "reference": ref,
                "regenerated": new,
                "rel_diff": rel,
                "within": bool(abs(rel) <= rel_tol),
            }
        )
    calibrated = bool(lim.meta.get("calibrated", False))
    return {
        "rho": rho,
        "rel_tol": rel_tol,
        "calibrated": calibrated,
        "rows": rows,
        "passed": calibrated and all(r["within"] for r in rows),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("limits", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    res = repro_gate(Table.read(args.limits), MACS0416_D052)
    args.out.write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({"passed": res["passed"]}))
    return 0 if res["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
