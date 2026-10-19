"""Gate V2 (PLAN.md addendum): edge.edge_pair vs per-epoch forced stamp cuts, injection set.
Usage (repo root): gate_v2.py GAIA_ROOT OUT FIELD...
Writes OUT/<field>/pairs.ecsv (one row per proxy-covered injection pair) and, with all three fields,
results/dimming/edge_coverage/gate_v2.json.
"""

import importlib.util
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack

sys.path.insert(0, "scripts")
import dimming_screen as ds  # noqa: E402

from jwst_anomaly import cutouts, pipeline  # noqa: E402

_spec = importlib.util.spec_from_file_location("edge", Path(__file__).with_name("edge.py"))
ed = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ed)

STEP = 32
gaia_root, out_root, fields = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3:]
p = ds.Params()
tabs = []
for name in fields:
    t0 = time.time()
    out = out_root / name
    out.mkdir(parents=True, exist_ok=True)
    field = ds.load_field(ds.DEFAULT_CONFIG, name)
    gaia = Table.read(gaia_root / name / "gaia.ecsv")
    cats, lc, scale, err, res, ordn = ds._screen(field, p, gaia)
    inj = ds.monitored(lc, err, ordn, p) & ~ds.flagged(res)
    idx = np.flatnonzero(inj)
    covered = np.isfinite(lc["flux"][:, :, 0])
    targets = Table(
        {
            "source_uid": [f"{name}-d{i:05d}" for i in idx],
            "ra": np.asarray(lc["ra"])[idx],
            "dec": np.asarray(lc["dec"])[idx],
        }
    )
    epochs = [(k, e) for k, e in enumerate(field["epochs"]) if "F200W" in e["bands"]]

    def one(ke, field=field, targets=targets, out=out):
        k, e = ke
        uri = pipeline.l3_image_uri(field["cloud"], ds.obs_id(e, "F200W"))
        grid, wcs = ed.wht_grid(uri, STEP)
        rel = ed.wht_rel(grid, wcs, targets["ra"], targets["dec"], STEP)
        ct = cutouts.make_cutouts(
            uri, targets, size_arcsec=2.0, extensions=["WHT"], out_dir=out / "stamps" / f"e{k}"
        )
        return k, rel, ct

    with ThreadPoolExecutor(max_workers=4) as pool:
        res_k = list(pool.map(one, epochs))
    rows = []
    for k, rel, ct in res_k:
        m = {str(a): r for a, r in zip(ct["source_uid"], ct, strict=True)}
        for j, i in enumerate(idx):
            if not covered[i, k]:
                continue
            r = m[str(targets["source_uid"][j])]
            rows.append(
                (
                    name,
                    str(targets["source_uid"][j]),
                    k,
                    float(rel[j]),
                    str(r["quality_flag"]),
                    bool(r["on_edge"]),
                    float(r["frac_nan"]),
                    float(r["wht_rel"]),
                )
            )
    t = Table(
        rows=rows,
        names=(
            "field",
            "uid",
            "epoch",
            "grid_wht_rel",
            "quality_flag",
            "on_edge",
            "frac_nan",
            "stamp_wht_rel",
        ),
    )
    t["edge_pair"] = ed.edge_pair(np.ones(len(t), bool), t["grid_wht_rel"])
    t["stamp_edge"] = ed.stamp_edge(
        t["quality_flag"], t["on_edge"], t["frac_nan"], t["stamp_wht_rel"]
    )
    t.meta["provenance"] = "derived from observed F200W i2d WHT / SCI"
    t.write(out / "pairs.ecsv", overwrite=True)
    print(
        name,
        len(t),
        int(t["edge_pair"].sum()),
        int(t["stamp_edge"].sum()),
        f"{time.time() - t0:.0f}s",
        flush=True,
    )
    tabs.append(t)
if set(fields) == {"macs0416", "abell2744", "nexus"}:
    allp = vstack(tabs)
    g = ed.gate_v2(allp["edge_pair"], allp["stamp_edge"])
    g["per_field"] = {
        f: ed.gate_v2(allp["edge_pair"][allp["field"] == f], allp["stamp_edge"][allp["field"] == f])
        for f in fields
    }
    s = json.loads(Path("results/dimming/edge_coverage/summary.json").read_text())
    j = s["joint_exposure_rho0.1"]
    g["decision"] = (
        ed.decision(j["veto_off_masked"], j["veto_off"])
        if g["verdict"] == "passed"
        else f"no_reading_gate_v2_{g['verdict']}"
    )
    g["provenance"] = "derived"
    Path("results/dimming/edge_coverage/gate_v2.json").write_text(json.dumps(g, indent=1) + "\n")
    print(json.dumps(g, indent=1))
