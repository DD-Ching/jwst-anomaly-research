"""D-039 persistence test on the compact single-epoch vanish flags (targets) and on every monitored
source in every epoch it is detected in (known on-sky controls; `on_sky` efficiency vs mosaic S/N).

Usage (repo root): run.py FIELD GAIA_ROOT OUT [targets] [epochs=K,K...]
("targets": the flags only, no monitored controls; used for NEXUS, see RESULT.md "Deviations";
"epochs=": only those epoch indices, to finish an interrupted field; merge with partial.py;
"max_nan=K": fill up to K NaN aperture pixels, results/dimming/nan_aperture/PLAN.md; default 0)
"""

import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack

sys.path.insert(0, "scripts")
import dimming_exposure_check as dec  # noqa: E402
import dimming_screen as ds  # noqa: E402
import persistence_check as pc  # noqa: E402

name, groot, out = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
targets_only = "targets" in sys.argv[4:]
only = next((a[7:] for a in sys.argv[4:] if a.startswith("epochs=")), None)
only = [int(k) for k in only.split(",")] if only else None
max_nan = int(next((a[8:] for a in sys.argv[4:] if a.startswith("max_nan=")), 0))
out.mkdir(parents=True, exist_ok=True)
p = ds.Params()
field = ds.load_field(ds.DEFAULT_CONFIG, name)
if only is not None and (
    len(set(only)) != len(only) or not all(0 <= k < len(field["epochs"]) for k in only)
):
    raise SystemExit(f"epochs= must be distinct indices 0..{len(field['epochs']) - 1}: {only}")
gaia = Table.read(groot / name / "gaia.ecsv")
cats, lc, scale, err, res, ordn = ds._screen(field, p, gaia)
static_veto = ds.catalogue_veto({**ordn, "single_epoch": np.zeros_like(ordn["single_epoch"])})
sel = ds.monitored(lc, err, ordn, p)
n_mon = int(sel.sum())
if targets_only:
    sel = np.zeros_like(sel)
det = lc["det"][:, :, 0]
target = ds.flagged(res) & ordn["single_epoch"] & ~static_veto & ordn["point_like"] & res["vanish"]
rows = []
for k in only if only is not None else range(len(field["epochs"])):  # given order is kept
    e = field["epochs"][k]
    if "F200W" not in e["bands"]:
        continue
    idx = np.flatnonzero((target | sel) & det[:, k])
    if not len(idx):
        continue
    uid = {f"{name}-d{i:05d}": i for i in idx}
    pos = {u: (float(lc["ra"][i]), float(lc["dec"][i])) for u, i in uid.items()}
    files = [f for d in dec.cal_dirs(e["obs"]) for f in pc.list_cal(d) if "_nrc" in f]
    per_exp, summ = pc.run(pos, files, workers=16, max_nan=max_nan)
    per_exp.write(out / f"exposures_e{k}.ecsv", overwrite=True)
    summ["epoch"] = k
    summ["kind"] = ["flag" if target[uid[u]] else "monitored" for u in summ["uid"]]
    summ["mosaic_snr"] = [float(lc["flux"][uid[u], k, 0] / err[uid[u], k, 0]) for u in summ["uid"]]
    summ["n_files"] = len(files)
    summ["n_failed"] = len(summ.meta.get("failed", []))
    summ.meta.pop("failed", None)
    rows.append(summ)
    print(
        name,
        k,
        len(files),
        dict(zip(*np.unique(summ["verdict"], return_counts=True), strict=True)),
        flush=True,
    )
t = (
    vstack(rows, metadata_conflicts="silent")
    if rows
    else Table(names=["uid", "verdict", "kind"], dtype=[str, str, str])
)
t.meta.update(provenance="derived", source=f"results/dimming/persistence_run/run.py {name}")
t.meta["max_nan"] = max_nan
t.write(out / "verdicts.ecsv", overwrite=True)
(out / "counts.json").write_text(
    json.dumps(
        {
            "field": name,
            "n_targets": int(target.sum()),
            "n_monitored": n_mon,
            "targets_only": targets_only,
            "max_nan": max_nan,
        },
        indent=1,
    )
    + "\n"
)
