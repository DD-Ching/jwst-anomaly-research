"""PLAN.md measurement: NEXUS monitored sources detected in the e2 (o014) mosaic with e2 S/N >= 50,
measured by `persistence_check.run` on every `_nrc` `_cal` file of jw05105-o015_t005
(NEXUS-Center-Tile-6), the same file convention as persistence_run/run.py.
Usage (repo root): run_o015.py GAIA_ROOT OUT max_nan=K [snr_min=X]
("snr_min=": post hoc override of the frozen e2 S/N cut, RESULT.md)"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "scripts")
import dimming_exposure_check as dec  # noqa: E402
import dimming_screen as ds  # noqa: E402
import persistence_check as pc  # noqa: E402
from astropy.table import Table  # noqa: E402

OBS, E2, SNR_E2_MIN = "jw05105-o015_t005", 2, 50.0  # PLAN.md: 10 / sqrt(1245 / 31137)

groot, out = Path(sys.argv[1]), Path(sys.argv[2])
max_nan = int(next((a[8:] for a in sys.argv[3:] if a.startswith("max_nan=")), 0))
snr_min = float(next((a[8:] for a in sys.argv[3:] if a.startswith("snr_min=")), SNR_E2_MIN))
out.mkdir(parents=True, exist_ok=True)
p = ds.Params()
field = ds.load_field(ds.DEFAULT_CONFIG, "nexus")
assert field["epochs"][E2]["obs"] == "jw05105-o014_t001"
gaia = Table.read(groot / "nexus" / "gaia.ecsv")
cats, lc, scale, err, res, ordn = ds._screen(field, p, gaia)
snr_e2 = lc["flux"][:, E2, 0] / err[:, E2, 0]
sel = ds.monitored(lc, err, ordn, p) & lc["det"][:, E2, 0] & (snr_e2 >= snr_min)
idx = np.flatnonzero(sel)
if not len(idx):
    raise SystemExit(f"no monitored e2 source with e2 S/N >= {snr_min}")
uid = {f"nexus-d{i:05d}": i for i in idx}
pos = {u: (float(lc["ra"][i]), float(lc["dec"][i])) for u, i in uid.items()}
files = [f for d in dec.cal_dirs(OBS) for f in pc.list_cal(d) if "_nrc" in f]
per_exp, summ = pc.run(pos, files, workers=16, max_nan=max_nan)
per_exp.write(out / "exposures.ecsv", overwrite=True)
summ["e2_snr"] = [float(snr_e2[uid[u]]) for u in summ["uid"]]
summ["n_files"] = len(files)
summ["n_failed"] = len(summ.meta.pop("failed", []))
summ = summ[summ["n_exposures"] > 0]
summ.meta.update(
    provenance="derived", source=f"results/dimming/o015_k1/run_o015.py {OBS}", max_nan=max_nan
)
summ.write(out / "verdicts.ecsv", overwrite=True)
(out / "counts.json").write_text(
    json.dumps(
        {
            "n_candidates": len(idx),
            "n_covered": len(summ),
            "n_files": len(files),
            "snr_min": snr_min,
            "e2_snr_max_monitored": float(np.nanmax(snr_e2[ds.monitored(lc, err, ordn, p)])),
        },
        indent=1,
    )
    + "\n"
)
print(len(idx), len(summ), dict(zip(*np.unique(summ["verdict"], return_counts=True), strict=True)))
