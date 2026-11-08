"""Gate A: D-037's repeated-pair and model-dependent-mu functions on Abell 2744, the second field
(plan: results/lens_consistency/d037_rules_a2744/PLAN.md). Offline apart from the catalogued
image list (Bergamini+23 Table 2, VizieR TSV passed with --images)."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import lens_consistency as lc
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

# hand readings (docs/fields/abell2744.md) -> wanted D-037 labels; key (model, system, forced_class)
HAND_A = {
    ("uncover2", "3.1", "absent"): ["model_topology"],
    ("uncover2", "34.1", "absent"): ["model_topology"],
    ("bergamini23", "700.1", "absent"): [],
    ("bergamini23", "18.1", "confused"): [],
    ("bergamini23", "22.1", "confused"): [],
}
N_HAND_ROWS = 6  # 700.1 has two absent rows
# Bergamini model: recorded validate shared_match (docs/fields/abell2744.md; fe.infn.it down)
SHARED_BERGAMINI = ["3.2a", "3.2b", "34.1a", "34.1b", "700.1a", "700.1b"]
# docs-table positions of the shared pairs (route check against Table 2)
DOCS_POS = {
    "3.2a": (3.58921, -30.39382),
    "3.2b": (3.58896, -30.39380),
    "34.1a": (3.59341, -30.41081),
    "34.1b": (3.59380, -30.41069),
    "700.1a": (3.57970, -30.35772),
    "700.1b": (3.57917, -30.35783),
}
ROUTE_POS_ARCSEC = 0.1  # docs vs Table 2 position agreement (ASSUMPTION: docs round to 1e-5 deg)
XM_ROW_MATCH_ARCSEC = 0.05  # same row in crossmodel.ecsv
# post hoc (route check failed as frozen): Table 2 names system 46 "JD1"; the pairing follows from
# each stored sep_image_arcsec (0.39 / 0.76 / 0.28") equalling its Table 2 separation
JD1_ALIAS = {"JD1B": "46.1a", "JD1A": "46.1b", "JD1C": "46.1c"}


def read_vizier_tsv(path: Path, alias: dict[str, str] | None = None) -> Table:
    """Bergamini+23 Table 2 (VizieR asu-tsv: header, units and dashes lines) as an image table."""
    lines = [s for s in path.read_text().splitlines() if s and not s.startswith("#")]
    head = lines[0].split("\t")
    rows = [s.split("\t") for s in lines[3:]]
    col = {h.strip(): i for i, h in enumerate(head)}
    ids = [(alias or {}).get(r[col["ID"]].strip(), r[col["ID"]].strip()) for r in rows]
    out = Table(
        {
            "image_id": ids,
            "system": [i.rstrip("abcdefghijklmnopqrstuvwxyz") for i in ids],
            "ra": [float(r[col["RAJ2000"]]) for r in rows],
            "dec": [float(r[col["DEJ2000"]]) for r in rows],
        }
    )
    out.meta.update(provenance="observed", source=f"{path.name} (VizieR J/ApJ/952/84/table2)")
    return out


def route_check(images: Table, ip_ids: list[str]) -> dict:
    ids = list(np.asarray(images["image_id"]).astype(str))
    missing = sorted(set(ip_ids) - set(ids))
    pos = {}
    for k, (ra, dec) in DOCS_POS.items():
        i = ids.index(k)
        sep = SkyCoord(ra, dec, unit="deg").separation(
            SkyCoord(images["ra"][i], images["dec"][i], unit="deg")
        )
        pos[k] = round(float(sep.arcsec), 3)
    ok = len(images) == 149 and not missing and max(pos.values()) <= ROUTE_POS_ARCSEC
    return {"n_images": len(images), "missing_ids": missing, "docs_pos_sep_arcsec": pos, "ok": ok}


def other_model(xm: Table, model: str, r) -> tuple[float | None, float | None]:
    sel = (np.asarray(xm["model"]).astype(str) == model) & (
        np.asarray(xm["system"]).astype(str) == str(r["system"])
    )
    for x in xm[sel]:
        sep = SkyCoord(r["ra"], r["dec"], unit="deg").separation(
            SkyCoord(x["ra"], x["dec"], unit="deg")
        )
        if sep.arcsec < XM_ROW_MATCH_ARCSEC:
            s, m = float(x["other_sep_arcsec"]), float(x["other_magnification"])
            return (s if np.isfinite(s) else None, m if np.isfinite(m) else None)
    return None, None


def label_rows(model: str, table: Table, images: Table, shared: list[str], xm: Table) -> list[dict]:
    if hasattr(table["forced_class"], "filled"):
        table["forced_class"] = [str(x) for x in table["forced_class"].filled("")]
    topo = lc.repeated_pair_topology(table, images, shared)
    out = []
    for r, t in zip(table, topo, strict=True):
        if r["forced_class"] not in lc.COPY_CLASSES:
            continue
        sep2, mu2 = other_model(xm, model, r)
        md = lc.mu_model_dependent(float(r["magnification"]), mu2, sep2)
        ratio = lc.mu_ratio(float(r["magnification"]), mu2)
        out.append(
            {
                "model": model,
                "system": str(r["system"]),
                "forced_class": str(r["forced_class"]),
                "radec": [round(float(r["ra"]), 6), round(float(r["dec"]), 6)],
                "mu_first": round(float(r["magnification"]), 2),
                "second_sep_arcsec": None if sep2 is None else round(sep2, 2),
                "mu_second_nearest": None if mu2 is None else round(mu2, 2),
                "mu_ratio": None if ratio is None else round(ratio, 2),
                "shared_pair": t,
                "mu_model_dependent": md,
                "labels": (["model_topology"] if t else [])
                + (["mu_model_dependent"] if md else []),
            }
        )
    return out


def gate(rows: list[dict]) -> dict:
    res, n = {}, 0
    for (model, sys_id, cls), want in HAND_A.items():
        got = [
            r["labels"]
            for r in rows
            if (r["model"], r["system"], r["forced_class"]) == (model, sys_id, cls)
        ]
        n += len(got)
        res[f"{model}:{sys_id}:{cls}"] = {
            "want": want,
            "got": got,
            "ok": bool(got) and all(g == want for g in got),
        }
    res["n_hand_rows"] = n
    res["unread_rows"] = [
        f"{r['model']}:{r['system']}:{r['forced_class']} {r['labels']}"
        for r in rows
        if (r["model"], r["system"], r["forced_class"]) not in HAND_A
    ]
    res["passes"] = n == N_HAND_ROWS and all(v["ok"] for v in res.values() if isinstance(v, dict))
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", type=Path, required=True, help="abell2744-uncover2 results directory")
    ap.add_argument("--images", type=Path, required=True, help="VizieR Table 2 TSV")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--post-hoc-alias", action="store_true", help="rename JD1A-C to 46.1b/a/c")
    a = ap.parse_args()
    images = read_vizier_tsv(a.images, JD1_ALIAS if a.post_hoc_alias else None)
    ip = Table.read(a.run / "imageplane_bergamini23.ecsv", format="ascii.ecsv")
    shared_unc = [str(i) for i in ip["image_id"][ip["shared_match"]]]
    route = route_check(images, [str(i) for i in ip["image_id"]])
    summary = {
        "provenance": "derived (forced classes: stored 2026-10-14 run; mu: model_prediction)",
        "inputs": {
            "run": str(a.run),
            "images": str(a.images),
            "images_sha256": hashlib.sha256(a.images.read_bytes()).hexdigest(),
        },
        "assumptions": {
            "pair_radius_arcsec": lc.PAIR_RADIUS_ARCSEC,
            "mu_model_factor": lc.MU_MODEL_FACTOR,
            "mu_second_radius_arcsec": lc.MU_SECOND_RADIUS_ARCSEC,
            "route_pos_arcsec": ROUTE_POS_ARCSEC,
        },
        "post_hoc_alias": JD1_ALIAS if a.post_hoc_alias else None,
        "route_check": route,
        "shared_match": {"uncover2": shared_unc, "bergamini23": SHARED_BERGAMINI},
    }
    if route["ok"]:
        xm = Table.read(a.run / "crossmodel.ecsv", format="ascii.ecsv")
        rows = []
        for model, shared in (("bergamini23", SHARED_BERGAMINI), ("uncover2", shared_unc)):
            t = Table.read(a.run / f"images_predicted_{model}.ecsv", format="ascii.ecsv")
            rows += label_rows(model, t, images, shared, xm)
        summary["flag_rows"] = rows
        summary["gate_A"] = gate(rows)
    a.out.mkdir(parents=True, exist_ok=True)
    name = "summary_post_hoc.json" if a.post_hoc_alias else "summary.json"
    (a.out / name).write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
