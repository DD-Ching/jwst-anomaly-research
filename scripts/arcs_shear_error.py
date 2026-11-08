"""Model shear error at the single-component pooled ``arcs`` anti rows (D-024 addendum 13).

Plan (frozen before the run): ``results/arcs/anti-shear-error/PLAN.md``. Re-classes the deciding
runs' stored rows with ``lens_consistency.orientation_table`` (unchanged) under the RELICS
Lenstool v1 ``range/`` realisations of RXC J0600 and two further RELICS v1 models of WHL0137
(Zitrin LTM-Gauss, WSLAP). Maps are streamed into a cache directory, sha256-recorded and deleted
after use. Predicted axes ``model_prediction``; classes and labels ``derived``.

    REQUESTS_CA_BUNDLE=<proxy bundle + certifi roots> python scripts/arcs_shear_error.py --cache DIR
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table

from jwst_anomaly import lensmodel

_DIR = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("lens_consistency", _DIR / "lens_consistency.py")
lc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(lc)

ROOT = _DIR.parent / "results" / "arcs"
OUT = ROOT / "anti-shear-error"
RELICS = "https://archive.stsci.edu/hlsps/relics/"
COSMO = (70.0, 0.3)

# frozen thresholds (ASSUMPTIONS, PLAN.md)
KAPPA_MAX_REL = 0.05
CONTROL_MIN = 0.70
ROBUST_F = 0.68
SENSITIVE_F = 0.32
RANGE_GATE_MAPS = (0, 49, 99)
N_RANGE = 100

WHL = {
    "row": 2194,
    "run": "whl0137-relics/whl0137-relics-lenstool/arcs.ecsv",
    "z_lens": 0.566,
    "models": {
        "zitrin-ltm-gauss": "whl0137m08/models/zitrin-ltm-gauss/v1/"
        "hlsp_relics_model_model_whl0137m08_zitrin-ltm-gauss_v1",
        "wslap": "whl0137m08/models/wslap/v1/hlsp_relics_model_model_whl0137m08_wslap_v1",
    },
}
RXC = {
    "row": 581,
    "run": "rxc0600-relics/rxc0600-relics-glafic/arcs.ecsv",
    "z_lens": 0.43,
    "range": "rxc0600m20/models/lenstool/range/"
    "hlsp_relics_model_model_rxc0600m20_lenstool-map{:03d}_v1",
}


# ---- rule functions (unit-tested before the run) ----
def map_gate(median_rel_diff: float, maximum: float = KAPPA_MAX_REL) -> bool:
    """Rule 1: our κ against the published κ map, median relative difference <= 0.05."""
    return bool(np.isfinite(median_rel_diff) and median_rel_diff <= maximum)


def retention(classes: list[str]) -> float:
    """Share of control rows still ``aligned``."""
    return sum(c == "aligned" for c in classes) / len(classes) if classes else float("nan")


def control_gate(r: float, minimum: float = CONTROL_MIN) -> bool:
    """Rule 2: control retention >= 0.70."""
    return bool(np.isfinite(r) and r >= minimum)


def models_label(classes: list[str]) -> str:
    """Rule 3 (WHL0137): classes of the row under the valid, gate-passing further models."""
    if not classes:
        return "undetermined"
    k = sum(c == "anti" for c in classes)
    if k == len(classes):
        return "robust"
    if k <= len(classes) / 2:
        return "model_sensitive"
    return "partly"


def range_label(f: float, valid: bool) -> str:
    """Rule 3 (RXC J0600): f = share of realisations in which the row is ``anti``."""
    if not valid or not np.isfinite(f):
        return "undetermined"
    if f >= ROBUST_F:
        return "robust"
    if f <= SENSITIVE_F:
        return "model_sensitive"
    return "partly"


def outcome(labels: dict[str, str]) -> str:
    """Rule 4 over the tested rows."""
    vals = list(labels.values())
    if vals and all(v == "robust" for v in vals):
        return "model error does not remove the tested rows"
    sens = [k for k, v in labels.items() if v == "model_sensitive"]
    if sens:
        return "model error explains " + ", ".join(sens)
    if "undetermined" in vals:
        return "undetermined"
    return "mixed"


# ---- I/O ----
def fetch(url: str, cache: Path) -> tuple[Path, str]:
    import requests

    cache.mkdir(parents=True, exist_ok=True)
    path = cache / url.rsplit("/", 1)[1]
    digest = hashlib.sha256()
    for attempt in range(4):
        try:
            with requests.get(url, stream=True, timeout=120) as r:
                r.raise_for_status()
                with open(path, "wb") as fh:
                    for chunk in r.iter_content(1 << 20):
                        fh.write(chunk)
                        digest.update(chunk)
            return path, digest.hexdigest()
        except requests.RequestException:
            if attempt == 3:
                raise
            digest = hashlib.sha256()
    raise RuntimeError("unreachable")


def build(prefix: str, z_lens: float, cache: Path, with_kappa: bool) -> tuple:
    """(model, sha256 per file, κ check or None); deletes the files after reading."""
    names = ["x-arcsec-deflect", "y-arcsec-deflect"] + (["kappa"] if with_kappa else [])
    files, shas = {}, {}
    for n in names:
        files[n], shas[n] = fetch(f"{RELICS}{prefix}_{n}.fits", cache)
    model = lensmodel.MapLensModel.from_fits(
        files["x-arcsec-deflect"],
        files["y-arcsec-deflect"],
        z_lens,
        FlatLambdaCDM(H0=COSMO[0], Om0=COSMO[1]),
        source=prefix,
    )
    check = lc.map_check(model, files["kappa"], "kappa", None) if with_kappa else None
    for p in files.values():
        p.unlink()
    if check:
        check["map"] = Path(check["map"]).name
    return model, shas, check


def reclass(model, rows: Table) -> Table:
    return lc.orientation_table(model, rows)


def split_rows(run: str, label: int) -> tuple[Table, Table]:
    """(the tested row, the deciding run's aligned strong-shear background control rows)."""
    t = Table.read(ROOT / run)
    bg = t[t["strong_shear"] & (t["population"] == "background")]
    row = bg[bg["label"] == label]
    if len(row) != 1 or row["orientation_class"][0] != "anti":
        raise SystemExit(f"error: row {label} is not a single anti row of {run}")
    return row, bg[bg["orientation_class"] == "aligned"]


def _rec(t: Table) -> dict:
    return {
        "class": str(t["orientation_class"][0]),
        "offset_min": round(float(t["offset_min"][0]), 1),
        "offset_max": round(float(t["offset_max"][0]), 1),
        "g_min": round(float(t["g_min"][0]), 3),
    }


def run_whl(cache: Path) -> dict:
    row, ctrl = split_rows(WHL["run"], WHL["row"])
    per = {}
    for name, prefix in WHL["models"].items():
        model, shas, check = build(prefix, WHL["z_lens"], cache, True)
        valid = map_gate(check["median_rel_diff"])
        r = retention(list(reclass(model, ctrl)["orientation_class"]))
        per[name] = {
            "sha256": shas,
            "kappa_check": check,
            "map_gate": valid,
            "control_n": len(ctrl),
            "control_retention": round(r, 3),
            "control_gate": control_gate(r),
            "row": _rec(reclass(model, row)),
        }
    used = [p["row"]["class"] for p in per.values() if p["map_gate"] and p["control_gate"]]
    return {
        "label": WHL["row"],
        "deciding_run_row": _rec(row),
        "models": per,
        "label_rule3": models_label(used),
    }


def run_rxc(cache: Path) -> dict:
    row, ctrl = split_rows(RXC["run"], RXC["row"])
    real, gates, rets = [], {}, []
    for i in range(N_RANGE):
        with_k = i in RANGE_GATE_MAPS
        model, shas, check = build(RXC["range"].format(i), RXC["z_lens"], cache, with_k)
        if with_k:
            gates[f"map{i:03d}"] = {
                "kappa_check": check,
                "passes": map_gate(check["median_rel_diff"]),
            }
        rets.append(retention(list(reclass(model, ctrl)["orientation_class"])))
        rec = _rec(reclass(model, row))
        rec["map"] = i
        rec["sha256_x"] = shas["x-arcsec-deflect"]
        real.append(rec)
    valid_maps = all(g["passes"] for g in gates.values())
    r_med = float(np.median(rets))
    f = sum(r["class"] == "anti" for r in real) / len(real)
    valid = valid_maps and control_gate(r_med)
    pct = lambda k: [round(float(v), 3) for v in np.percentile([r[k] for r in real], (16, 50, 84))]  # noqa: E731
    return {
        "label": RXC["row"],
        "deciding_run_row": _rec(row),
        "map_gates": gates,
        "control_n": len(ctrl),
        "control_retention_median": round(r_med, 3),
        "control_retention_range": [round(min(rets), 3), round(max(rets), 3)],
        "control_gate": control_gate(r_med),
        "f_anti": f,
        "label_rule3": range_label(f, valid),
        "offset_min_pct": pct("offset_min"),
        "g_min_pct": pct("g_min"),
        "realisations": real,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache", type=Path, required=True)
    args = ap.parse_args(argv)
    whl = run_whl(args.cache)
    rxc = run_rxc(args.cache)
    labels = {"WHL0137 2194": whl["label_rule3"], "RXC J0600 581": rxc["label_rule3"]}
    summary = {
        "plan": "results/arcs/anti-shear-error/PLAN.md",
        "params": {
            "kappa_max_rel": KAPPA_MAX_REL,
            "control_min": CONTROL_MIN,
            "robust_f": ROBUST_F,
            "sensitive_f": SENSITIVE_F,
            "range_gate_maps": list(RANGE_GATE_MAPS),
            "n_range": N_RANGE,
            "cosmology": list(COSMO),
            "orientation_table": "lens_consistency defaults (aligned 30, anti 60, z margin 0.1)",
        },
        "provenance": {"axes": "model_prediction", "classes": "derived"},
        "WHL0137": whl,
        "RXC J0600": rxc,
        "MACS0308 2157": "untested (plan: realisation maps 36 MB each; byte-range route next)",
        "labels": labels,
        "outcome": outcome(labels),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: summary[k] for k in ("labels", "outcome")}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
