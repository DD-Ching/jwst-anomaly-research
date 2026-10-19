"""Frozen two-band pair rule for review-list entries (``results/transients/review_list/vet_*``).

Per entry, Δ = flux_B − flux_A and σ = hypot(err_A, err_B) from ``hit_lightcurve.py`` rows of the
pair's own two mosaics in each band. Outcome (rule frozen in that PLAN.md; thresholds ASSUMPTIONs):

- ``untested``: a pair mosaic is missing in either band;
- ``no_change``: |Δ|/σ < 2 in both bands;
- ``change``: |Δ|/σ ≥ 2 in both bands with the same sign;
- ``inconclusive``: anything else (one band only, or opposite signs).

    python scripts/review_pair_rule.py --lc vet/lightcurves.ecsv \\
        --entry h039 jw02736-o001_t001_nircam_clear jw06882-o057_t057_nircam_clear \\
        --out vet/outcome.json
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from astropy.table import Table

BANDS = ("F150W", "F444W")
Z_MIN = 2.0  # ASSUMPTION: the frozen per-band significance threshold


def classify(z: dict[str, float | None]) -> str:
    """Outcome from the per-band z = Δ/σ (``None`` = a pair mosaic is missing)."""
    if any(z.get(b) is None for b in BANDS):
        return "untested"
    vals = [z[b] for b in BANDS]
    if all(abs(v) < Z_MIN for v in vals):
        return "no_change"
    if all(abs(v) >= Z_MIN for v in vals) and len({math.copysign(1, v) for v in vals}) == 1:
        return "change"
    return "inconclusive"


def pair_z(lc: Table, name: str, prefix_a: str, prefix_b: str) -> dict[str, dict | None]:
    """Per band: the pair's fluxes, Δ, σ and z for one entry (``None`` if a mosaic is missing)."""
    out: dict[str, dict | None] = {}
    rows = lc[lc["name"] == name]
    for band in BANDS:
        ra = [r for r in rows if r["obs_id"] == f"{prefix_a}-{band.lower()}"]
        rb = [r for r in rows if r["obs_id"] == f"{prefix_b}-{band.lower()}"]
        if not ra or not rb:
            out[band] = None
            continue
        d = float(rb[0]["flux_njy"] - ra[0]["flux_njy"])
        s = math.hypot(float(ra[0]["err_njy"]), float(rb[0]["err_njy"]))
        out[band] = {
            "flux_a_njy": round(float(ra[0]["flux_njy"]), 2),
            "flux_b_njy": round(float(rb[0]["flux_njy"]), 2),
            "delta_njy": round(d, 2),
            "sigma_njy": round(s, 2),
            "z": round(d / s, 3),
        }
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--lc", type=Path, required=True)
    ap.add_argument(
        "--entry", nargs=3, action="append", required=True, metavar=("NAME", "OBS_A", "OBS_B")
    )
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    lc = Table.read(args.lc)
    res = {}
    for name, a, b in args.entry:
        bands = pair_z(lc, name, a, b)
        z = {k: (v["z"] if v else None) for k, v in bands.items()}
        res[name] = {"pair": [a, b], "bands": bands, "outcome": classify(z)}
        print(name, z, res[name]["outcome"])
    args.out.write_text(json.dumps(res, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
