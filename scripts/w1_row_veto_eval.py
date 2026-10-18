"""K1r, the cost check and the adoption decision of the W1 row/column-segment veto (frozen rule:
results/w1_row_veto/PLAN.md).

    python scripts/w1_row_veto_eval.py --masks <union-mask dir> --runs <inject_shear out dir> \
        --variant A --out results/w1_row_veto/adoption_A.json

``--runs`` holds ``<field>/summary.json`` of the z0.3 screens made with the union masks; the
reference screens are ``results/w1_shear_ncf/summaries/z0.3_<field>.json`` (D-097 masks).
Masks and counts ``derived``; the 75 % threshold is an ASSUMPTION.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w1_streak_mask as wsm  # noqa: E402

NCF = ("macs0416-ncf", "macs1149-ncf", "abell370-ncf", "macs0417-ncf", "macs1423-ncf")
REF = Path(__file__).resolve().parents[1] / "results" / "w1_shear_ncf" / "summaries"
K1R_FIELD = "macs0417-ncf"


def adopt(k1r: bool, cost: dict[str, bool]) -> bool:
    """Frozen rule: adopted only if K1r holds and the cost check passes in every NCF."""
    return bool(k1r) and len(cost) == len(NCF) and all(cost.values())


def evaluate(mask_dir: Path, run_dir: Path, ref_dir: Path = REF) -> dict:
    ref = json.loads((ref_dir / f"z0.3_{K1R_FIELD}.json").read_text())["screen"]
    mask, wcs = wsm.read_mask(mask_dir / f"{K1R_FIELD}_streak_mask.fits.gz")
    ra, dec = ref["s_max_radec"]
    k1r = bool(wsm.veto_rows([ra], [dec], mask, wcs)[0])
    fields = {}
    for f in NCF:
        old = json.loads((ref_dir / f"z0.3_{f}.json").read_text())
        new = json.loads((run_dir / f / "summary.json").read_text())
        if new["field"] != f or "streak_veto" not in new["assumptions"]:
            raise ValueError(f"{run_dir / f}: not a {f} summary run with streak masks")
        fields[f] = {
            "kept_d097": wsm.kept_rows(old),
            "kept_union": wsm.kept_rows(new),
            "kept_frac": wsm.kept_rows(new) / wsm.kept_rows(old),
            "cost_ok": wsm.cost_ok(new, old),
            "union_mask": new["assumptions"]["streak_mask"],
            "screen_union": {
                k: new["screen"][k] for k in ("s_max", "s_max_radec", "p_random_max", "b_max")
            }
            | {"passes_eb_rule": new["screen"]["passes_eb_rule"]},
            "screen_d097": {
                k: old["screen"][k] for k in ("s_max", "p_random_max", "b_max", "passes_eb_rule")
            },
        }
    cost = {f: v["cost_ok"] for f, v in fields.items()}
    return {
        "k1r_field": K1R_FIELD,
        "k1r_radec": [ra, dec],
        "K1r": k1r,
        "cost_ok": cost,
        "adopted": adopt(k1r, cost),
        "fields": fields,
        "rule": "results/w1_row_veto/PLAN.md (K1r, cost check)",
        "provenance": {"masks": "derived", "counts": "derived", "keep_frac": "assumption"},
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--masks", type=Path, required=True)
    ap.add_argument("--runs", type=Path, required=True)
    ap.add_argument("--variant", required=True, choices=("A", "B"))
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    res = {"variant": args.variant, **evaluate(args.masks, args.runs)}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: res[k] for k in ("variant", "K1r", "cost_ok", "adopted")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
