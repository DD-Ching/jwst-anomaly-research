"""What-if: per-copy single_epoch veto off (D-039 swap scoping).

Usage (repo root): whatif.py FIELD INFL_F200W INFL_F444W OUT
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "scripts")
import dimming_screen as ds

field_name, i1, i2, out = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), Path(sys.argv[4])
out.mkdir(parents=True, exist_ok=True)
p = ds.Params()
field = ds.load_field(ds.DEFAULT_CONFIG, field_name)
gaia = ds._gaia(field, out)
cats, lc, scale, err, res, ordn = ds._screen(field, p, gaia)
sel = ds.monitored(lc, err, ordn, p)
times = ds._times_yr(field)
static_veto = ds.catalogue_veto({**ordn, "single_epoch": np.zeros_like(ordn["single_epoch"])})
infl = np.array([i1, i2])
covd = np.isfinite(lc["flux"][:, :, 0]).sum(axis=1)
area = ds.area_deg2(lc["ra"], lc["dec"], covd)

# workload: real flags vetoed by single_epoch alone
fl = ds.flagged(res)
other = ds.catalogue_veto({**ordn, "single_epoch": np.zeros_like(ordn["single_epoch"])})
se_only = fl & ordn["single_epoch"] & ~other
work = {
    "flagged": int(fl.sum()),
    "single_epoch_only": int(se_only.sum()),
    "single_epoch_only_vanish": int((se_only & res["vanish"]).sum()),
    "single_epoch_only_point_like": int((se_only & ordn["point_like"]).sum()),
    "single_epoch_only_point_like_vanish": int(
        (se_only & ordn["point_like"] & res["vanish"]).sum()
    ),
}

orig = ds._recovered


def no_se(fi, ei, det, p, inflation):
    r, _ = ds._classify_injected(fi, ei, det, p)
    rec = ds.flagged(r)
    if inflation is not None:
        r2, _ = ds._classify_injected(fi, ei * np.asarray(inflation)[None, None, :], det, p)
        rec &= ds.flagged(r2)
    return rec


lims = {}
for tag, fn in (("veto_on", orig), ("veto_off", no_se)):
    ds._recovered = fn
    eff = ds.efficiency_table(
        lc,
        err,
        static_veto,
        sel,
        times,
        ds.T_E_GRID,
        ds.RHO_GRID,
        ds.DEPTHS,
        p,
        seed=0,
        forced_inflation=infl,
    )
    first = eff[
        (eff["model"] == "w3") & (eff["t_e_yr"] == eff["t_e_yr"][0]) & (eff["rho"] == eff["rho"][0])
    ]
    n = int(np.sum(first["n_sources"]))
    lim = ds.rate_limits(eff, area, n)
    eff.write(out / f"efficiency_{tag}.ecsv", overwrite=True)
    lim.write(out / f"limits_{tag}.ecsv", overwrite=True)
    lims[tag] = lim
ds._recovered = orig
(out / "workload.json").write_text(
    json.dumps({"field": field_name, "n_monitored": int(sel.sum()), **work}, indent=1) + "\n"
)
print(json.dumps(work))
print(lims["veto_on"].colnames)
