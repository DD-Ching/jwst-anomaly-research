"""Held-out lens-model residuals across fields and models: shared / isolated misses, floors.

Frozen plan: ``results/lens_consistency/heldout_synthesis/PLAN.md``. Reads the stored per-image
held-out tables of D-034 addenda ap, at, au, aw, ax, ay; no data access.

    python scripts/heldout_synthesis.py [--out DIR]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from astropy.table import Table

ROOT = Path("results/lens_consistency")
OUT = ROOT / "heldout_synthesis"
MISS = 3.0  # arcsec, ASSUMPTION (gate H / K threshold)
GOOD = 1.0  # arcsec, ASSUMPTION (family members that fit)

# (field, model, table, primary)
TABLES = [
    ("Abell 370", "MARS", "mars_hff_heldout/heldout_A370.ecsv", True),
    ("MACS0416", "MARS", "mars_hff_heldout/heldout_MACSJ0416.ecsv", True),
    ("MACS1149", "MARS", "mars_m1149_heldout/heldout_MACSJ1149.ecsv", True),
    ("MACS1149", "CATS v4.1", "m1149_hst_heldout/heldout_macs1149-cats.ecsv", True),
    ("MACS1149", "Sharon v4cor", "m1149_hst_heldout/heldout_macs1149-sharon.ecsv", True),
    ("SMACS 0723", "glafic v2", "smacs_relics_heldout/heldout_glafic_v2.ecsv", True),
    (
        "SMACS 0723",
        "zitrin-ltm-gauss v2",
        "smacs_relics_heldout/heldout_zitrin-ltm-gauss_v2.ecsv",
        True,
    ),
    ("SMACS 0723", "lenstool v1", "smacs_relics_heldout/heldout_lenstool_v1.ecsv", False),
    ("Abell 2744", "Bergamini+23", "a2744_heldout/heldout_abell2744-bergamini23.ecsv", True),
    ("Abell 2744", "MARS", "a2744_heldout/heldout_A2744.ecsv", False),
    ("El Gordo", "Caminha+23", "elgordo_heldout/heldout_elgordo-caminha23.ecsv", True),
]


def offsets(t: Table) -> dict[str, tuple[str, float]]:
    """``image_id -> (system, offset)`` for rows inside the map; no prediction is +inf."""
    out = {}
    for r in t:
        if not bool(r["inside"]):
            continue
        d = float(r["dtheta_arcsec"])
        if int(r["n_predicted"]) == 0 or not np.isfinite(d):
            d = np.inf
        out[str(r["image_id"])] = (str(r["system"]), d)
    return out


def model_stats(off: dict[str, tuple[str, float]]) -> dict:
    """Rule 1: n, rms of finite offsets, median / 68 / 95 percentiles (inf counted), > MISS."""
    d = np.array([v for _, v in off.values()], dtype=float)
    if d.size == 0:
        return {"n": 0}
    fin = d[np.isfinite(d)]
    s = np.sort(d)

    def pct(q: float) -> float:
        # nearest-rank percentile, so an inf entry stays inf instead of turning into nan
        return float(s[max(0, int(np.ceil(q / 100 * s.size)) - 1)])

    return {
        "n": int(d.size),
        "n_no_prediction": int((~np.isfinite(d)).sum()),
        "rms_arcsec": float(np.sqrt(np.mean(fin**2))) if fin.size else None,
        "median_arcsec": pct(50),
        "p68_arcsec": pct(68),
        "p95_arcsec": pct(95),
        "frac_over_3": float(np.mean(d > MISS)),
    }


def isolated_misses(off: dict[str, tuple[str, float]]) -> list[str]:
    """Rule 2: offset > MISS while >= 2 other tested images of its family are all <= GOOD."""
    fam: dict[str, list[tuple[str, float]]] = {}
    for iid, (sysid, d) in off.items():
        fam.setdefault(sysid, []).append((iid, d))
    out = []
    for members in fam.values():
        for iid, d in members:
            others = [o for j, o in members if j != iid]
            if d > MISS and len(others) >= 2 and all(o <= GOOD for o in others):
                out.append(iid)
    return sorted(out)


def shared_misses(per_model: dict[str, dict[str, tuple[str, float]]]) -> dict[str, list[str]]:
    """Rule 3: image_id -> models (>= 2) missing it; ``per_model`` is one field's primaries."""
    hit: dict[str, list[str]] = {}
    for m, off in per_model.items():
        for iid, (_, d) in off.items():
            if d > MISS:
                hit.setdefault(iid, []).append(m)
    return {k: sorted(v) for k, v in sorted(hit.items()) if len(v) >= 2}


def flags(shared: dict[str, list[str]], isolated: dict[str, list[str]]) -> list[str]:
    """Rule 4: a shared miss that is isolated in at least one of the models that miss it."""
    return sorted(i for i, ms in shared.items() if any(i in isolated.get(m, []) for m in ms))


def field_floor(stats: dict[str, dict]) -> dict:
    """Rule 5: the 95th percentile of the primary model with the lowest rms."""
    ok = {m: s for m, s in stats.items() if s.get("n") and s.get("rms_arcsec") is not None}
    if not ok:
        return {"model": None, "p95_arcsec": None}
    best = min(ok, key=lambda m: ok[m]["rms_arcsec"])
    return {
        "model": best,
        "p95_arcsec": ok[best]["p95_arcsec"],
        "rms_arcsec": ok[best]["rms_arcsec"],
    }


def _json(x):
    if isinstance(x, float) and not np.isfinite(x):
        return "inf"
    if isinstance(x, dict):
        return {k: _json(v) for k, v in x.items()}
    if isinstance(x, list):
        return [_json(v) for v in x]
    return x


def run(root: Path = ROOT) -> dict:
    fields: dict[str, dict] = {}
    for field, model, rel, primary in TABLES:
        off = offsets(Table.read(root / rel, format="ascii.ecsv"))
        f = fields.setdefault(field, {"models": {}, "_off": {}})
        f["models"][model] = {
            "table": str(root / rel),
            "primary": primary,
            "stats": model_stats(off),
            "isolated_misses": isolated_misses(off),
        }
        if primary:
            f["_off"][model] = off
    for f in fields.values():
        prim = f.pop("_off")
        f["shared_misses"] = shared_misses(prim)
        iso = {m: v["isolated_misses"] for m, v in f["models"].items() if v["primary"]}
        f["flags"] = flags(f["shared_misses"], iso)
        f["floor"] = field_floor({m: f["models"][m]["stats"] for m in prim})
    return {
        "plan": str(OUT / "PLAN.md"),
        "provenance": "offsets: model_prediction of observed positions; thresholds: ASSUMPTION",
        "fields": fields,
        "flags": sorted(f"{k}: {i}" for k, f in fields.items() for i in f["flags"]),
    }


# Post hoc (not in the frozen plan): rule 2 with the model's *seen* images of the family counted as
# members, where the stage stored their offsets (``seen_rows`` in its summary.json).
SEEN = [
    ("Abell 2744", "Bergamini+23", "a2744_heldout", ("models", "abell2744-bergamini23")),
    ("El Gordo", "Caminha+23", "elgordo_heldout", ("all_heldout (reported only)",)),
    ("MACS1149", "CATS v4.1", "m1149_hst_heldout", ("models", "macs1149-cats")),
]


def seen_inclusive_isolated(
    off: dict[str, tuple[str, float]], seen: dict[str, float], family
) -> list[str]:
    """Post hoc rule 2 with seen images (offset = their stored fit residual) as family members."""
    both = dict(off)
    for iid, d in seen.items():
        both.setdefault(iid, (family(iid), float(d)))
    return sorted(i for i in isolated_misses(both) if i in off)


def posthoc(root: Path = ROOT) -> dict:
    from jwst_anomaly.lensmodel import image_family

    rel = {(f, m): r for f, m, r, _ in TABLES}
    out = {"label": "post hoc; not in the frozen plan", "seen_inclusive_isolated": {}}
    for field, model, d, keys in SEEN:
        node = json.loads((root / d / "summary.json").read_text())
        for k in keys:
            node = node[k]
        off = offsets(Table.read(root / rel[(field, model)], format="ascii.ecsv"))
        iso = seen_inclusive_isolated(off, node["seen_rows"], image_family)
        out["seen_inclusive_isolated"][f"{field} / {model}"] = {
            "n_seen": len(node["seen_rows"]),
            "isolated": {i: off[i][1] for i in iso},
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()
    res = run()
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "summary.json").write_text(json.dumps(_json(res), indent=1) + "\n")
    for k, f in res["fields"].items():
        print(f"== {k}: floor {f['floor']}  shared {f['shared_misses']}  flags {f['flags']}")
        for m, v in f["models"].items():
            s = v["stats"]
            print(f"  {m:22s} prim={v['primary']!s:5s} {s}  isolated={v['isolated_misses']}")
    print("FLAGS:", res["flags"])
    ph = posthoc()
    (a.out / "posthoc.json").write_text(json.dumps(_json(ph), indent=1) + "\n")
    print("POST HOC:", ph)


if __name__ == "__main__":
    main()
