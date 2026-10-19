"""Compare veto on/off limits; reproduces RESULT.md tables.

Run from the repo root: gaia_esa.py ROOT FIELD...; whatif.py FIELD INFL1 INFL2 ROOT/FIELD;
compare.py ROOT.
"""

import sys

import numpy as np
from astropy.table import Table

S = sys.argv[1]  # the whatif.py output root
com = {
    "macs0416": "results/dimming/macs0416_regen/limits.ecsv",
    "abell2744": "results/dimming/abell2744_calib/limits.ecsv",
    "nexus": "results/dimming/nexus_calib/limits.ecsv",
}
out = {}
for f, c in com.items():
    on = Table.read(f"{S}/{f}/limits_veto_on.ecsv")
    off = Table.read(f"{S}/{f}/limits_veto_off.ecsv")
    cm = Table.read(c)
    rep = np.allclose(
        np.asarray(on["exposure_source_yr"]),
        np.asarray(cm["exposure_source_yr"]),
        rtol=1e-9,
        equal_nan=True,
    )
    print(f, "reproduces committed:", rep)
    m = on["rho"] == 0.1
    for a, b in zip(on[m], off[m]):
        print(
            f"  tE {a['t_e_yr']:5} exp_on {a['exposure_source_yr']:8.3f} off {b['exposure_source_yr']:8.3f} ratio {b['exposure_source_yr'] / max(a['exposure_source_yr'], 1e-12):.3f}"
        )
    eon = Table.read(f"{S}/{f}/efficiency_veto_on.ecsv")
    eoff = Table.read(f"{S}/{f}/efficiency_veto_off.ecsv")
    for d in (0.2, 0.5, 1.0):
        k = (eon["model"] == "dimming") & (eon["depth"] == d)
        print(
            f"  dimming {d}: rec on {int(eon[k]['n_recovered'].sum())}/{int(eon[k]['n_injected'].sum())} off {int(eoff[k]['n_recovered'].sum())}"
        )
    # joint
    out[f] = (on, off)
tot = {}
for tag, i in (("on", 0), ("off", 1)):
    s = sum(np.asarray(out[f][i][out[f][i]["rho"] == 0.1]["exposure_source_yr"]) for f in out)
    tot[tag] = s
print("joint exposure rho0.1 on", np.round(tot["on"], 2))
print("off", np.round(tot["off"], 2))
print("gain", np.round(tot["off"] / tot["on"], 3))
