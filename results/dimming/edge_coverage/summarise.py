"""Joint exposures per variant and the frozen decision (PLAN.md). Usage (repo root):
summarise.py OUT   (OUT from run.py; writes results/dimming/edge_coverage/summary.json)
"""

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

_spec = importlib.util.spec_from_file_location("edge", Path(__file__).with_name("edge.py"))
ed = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ed)

FIELDS = ("macs0416", "abell2744", "nexus")
VARIANTS = ("veto_on", "veto_off", "veto_on_masked", "veto_off_masked")
COMMITTED = {
    "macs0416": "results/dimming/macs0416_regen/limits.ecsv",
    "abell2744": "results/dimming/abell2744_calib/limits.ecsv",
    "nexus": "results/dimming/nexus_calib/limits.ecsv",
}
root = Path(sys.argv[1])
out = {"fields": {}, "joint_exposure_rho0.1": {}}
t_e = None
joint = {v: 0.0 for v in VARIANTS}
for f in FIELDS:
    s = json.loads((root / f / "edge_summary.json").read_text())
    lim = {v: Table.read(root / f / f"limits_{v}.ecsv") for v in VARIANTS}
    m = lim["veto_on"]["rho"] == 0.1
    t_e = [float(x) for x in lim["veto_on"]["t_e_yr"][m]]
    cm = Table.read(COMMITTED[f])
    s["reproduces_committed_veto_on"] = bool(
        np.allclose(
            np.asarray(lim["veto_on"]["exposure_source_yr"]),
            np.asarray(cm["exposure_source_yr"]),
            rtol=1e-9,
            equal_nan=True,
        )
    )
    s["exposure_rho0.1"] = {
        v: [float(x) for x in lim[v]["exposure_source_yr"][lim[v]["rho"] == 0.1]] for v in VARIANTS
    }
    dim = {}
    for v in VARIANTS:
        e = Table.read(root / f / f"efficiency_{v}.ecsv")
        k = (e["model"] == "dimming") & (e["depth"] == 1.0)
        dim[v] = [int(e[k]["n_recovered"].sum()), int(e[k]["n_injected"].sum())]
    s["dimming_1.0_recovered_injected"] = dim
    for v in VARIANTS:
        joint[v] = joint[v] + np.array(s["exposure_rho0.1"][v])
    out["fields"][f] = s
out["t_e_yr"] = t_e
out["joint_exposure_rho0.1"] = {v: [float(x) for x in joint[v]] for v in VARIANTS}
out["veto_off_masked_over_unmasked"] = [
    float(a / b) for a, b in zip(joint["veto_off_masked"], joint["veto_off"], strict=True)
]
out["gain_veto_off_over_on_unmasked"] = [
    float(a / b) for a, b in zip(joint["veto_off"], joint["veto_on"], strict=True)
]
out["gain_veto_off_over_on_masked"] = [
    float(a / b) for a, b in zip(joint["veto_off_masked"], joint["veto_on_masked"], strict=True)
]
gv = out["fields"]["nexus"]
out["gate_v_passed"] = bool(gv["gate_v_passed"])
out["decision"] = (
    ed.decision(joint["veto_off_masked"], joint["veto_off"])
    if out["gate_v_passed"]
    else "no_reading_gate_v_failed"
)
out["provenance"] = "derived (exposures from simulated injections; weight maps observed)"
Path("results/dimming/edge_coverage/summary.json").write_text(json.dumps(out, indent=1) + "\n")
print(json.dumps({k: v for k, v in out.items() if k != "fields"}, indent=1))
