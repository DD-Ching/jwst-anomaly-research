"""Builds vetting.json from the committed exposure tables (run from the repository root)."""

import json
from pathlib import Path

import numpy as np
from astropy.table import Table

here = Path("results/dimming/nexus_vet")
single = json.loads((here / "single_points/exposure_summary.json").read_text())
sp = Table.read(here / "single_points/exposure_fluxes.ecsv")
d7 = Table.read(here / "d73390/exposure_fluxes.ecsv")
calib = Table.read("results/dimming/nexus_calib/survivors_and_flagged_controls.ecsv")


def exposures(uid, k, band):
    s = sp[(sp["uid"] == uid) & (sp["epoch"] == k) & (sp["band"] == band)]
    return [
        {
            "exposure": str(r["exposure"]),
            "flux_ujy": round(float(r["flux_ujy"]), 4),
            "n_jump": int(r["n_jump"]),
            "n_ap": int(r["n_ap"]),
        }
        for r in s
    ]


out = {
    "provenance": "derived",
    "post_hoc": ["chromatic_veto (rule written after the plan; unit-tested)"],
}
for uid, k in [
    ("nexus-d53823", 0),
    ("nexus-d60156", 3),
    ("nexus-d79928", 6),
    ("ctl0276", 1),
    ("ctl0290", 4),
    ("ctl0294", 5),
    ("ctl0298", 7),
]:
    out[f"{uid}|e{k}"] = {
        "forced": single[f"{uid}|e{k}|forced"],
        "F200W_exposures": exposures(uid, k, "F200W"),
    }

per_epoch = {}
for band in ("F200W", "F444W"):
    s = d7[d7["band"] == band]
    med = {
        int(k): float(np.nanmedian(s[s["epoch"] == k]["flux_ujy"])) for k in sorted(set(s["epoch"]))
    }
    med = {k: v for k, v in med.items() if np.isfinite(v)}
    vals = np.array(list(med.values()))
    if vals.min() <= 0:
        raise ValueError(f"{band}: non-positive epoch flux; a max/min ratio is meaningless")
    per_epoch[band] = {
        "epoch_median_cal_ujy_r015": {str(k): round(v, 3) for k, v in med.items()},
        "max_over_min": round(float(vals.max() / vals.min()), 3),
        "min_over_max": round(float(vals.min() / vals.max()), 3),
    }
r = calib[calib["uid"] == "nexus-d73390"][0]
i2d = [float(r[f"F200W_e{i}_fflux"]) for i in range(8)]
fin = [v for v in i2d if np.isfinite(v)]
if min(fin) <= 0:
    raise ValueError("d73390 i2d forced F200W: non-positive epoch flux")
per_epoch["F200W_i2d_forced_max_over_min"] = round(max(fin) / min(fin), 3)
out["nexus-d73390"] = per_epoch

r = calib[calib["uid"] == "nexus-d79956"][0]
out["nexus-d79956"] = {
    "F200W_forced": [
        None if not np.isfinite(r[f"F200W_e{i}_fflux"]) else round(float(r[f"F200W_e{i}_fflux"]), 3)
        for i in range(8)
    ],
    "F444W_forced": [
        None if not np.isfinite(r[f"F444W_e{i}_fflux"]) else round(float(r[f"F444W_e{i}_fflux"]), 3)
        for i in range(8)
    ],
    "cal_e6": exposures("nexus-d79956", 6, "F200W") + exposures("nexus-d79956", 6, "F444W"),
    "simbad_5arcsec": "no object",
    "ned": "not checked (query timed out after 60 s)",
    "tns": "not checked (HTTP 403 from the TNS load balancer for unauthenticated requests)",
}
(here / "vetting.json").write_text(json.dumps(out, indent=1) + "\n")
