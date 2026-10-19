"""Joint veto-off exposures of run_stamp.py (post hoc sensitivity). Usage (repo root):
summarise_stamp.py OUT   (OUT from run_stamp.py); writes gate_v2/stamp_sensitivity.json"""

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack

here = Path(__file__).parent
_spec = importlib.util.spec_from_file_location("edge", here / "edge.py")
ed = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ed)
root = Path(sys.argv[1])
joint = {}
for f in ("macs0416", "abell2744", "nexus"):
    for v in ("veto_on", "veto_off", "veto_on_masked", "veto_off_masked"):
        lim = Table.read(root / f / f"limits_{v}.ecsv")
        e = np.asarray(lim["exposure_source_yr"][lim["rho"] == 0.1], float)
        joint[v] = joint.get(v, 0.0) + e
pt = vstack(
    [Table.read(here / f"gate_v2/pairs_{f}.ecsv") for f in ("macs0416", "abell2744", "nexus")]
)
e, r = np.asarray(pt["edge_pair"], bool), np.asarray(pt["stamp_edge"], bool)
out = {
    "post_hoc": True,
    "pairs_grid_and_stamp": int((e & r).sum()),
    "pairs_grid_only": int((e & ~r).sum()),
    "pairs_stamp_only": int((~e & r).sum()),
    "pairs_neither": int((~e & ~r).sum()),
    "mask": "stamp_edge on injection-set pairs, edge_pair elsewhere",
    "joint_exposure_rho0.1": {k: [float(x) for x in v] for k, v in joint.items()},
    "veto_off_masked_over_unmasked": [
        float(x) for x in joint["veto_off_masked"] / joint["veto_off"]
    ],
    "gain_veto_off_over_on_masked": [
        float(x) for x in joint["veto_off_masked"] / joint["veto_on_masked"]
    ],
    "frozen_rule_reading": ed.decision(joint["veto_off_masked"], joint["veto_off"]),
    "provenance": "derived (simulated injections; weight maps observed)",
}
(here / "gate_v2/stamp_sensitivity.json").write_text(json.dumps(out, indent=1) + "\n")
print(json.dumps(out, indent=1))
