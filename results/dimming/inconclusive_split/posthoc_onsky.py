"""Post hoc (not in PLAN.md): how many `on_sky` verdicts rest on clean detections at one exposure
time only (SW + LW of one exposure), per field and kind.
Usage (repo root): posthoc_onsky.py -> results/dimming/inconclusive_split/posthoc_onsky.json"""

import glob
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).parent))
from split import MIN_SNR, n_times  # noqa: E402

RUN = Path("results/dimming/persistence_run")


def main() -> dict:
    out = {}
    for fld in ("nexus", "macs0416", "abell2744"):
        v = Table.read(RUN / fld / "verdicts.ecsv")
        res: dict[str, dict[str, int]] = {}
        for ep in sorted({int(x) for x in v["epoch"]}):
            (path,) = glob.glob(str(RUN / fld / f"exposures_e{ep}.ecsv*"))
            e = Table.read(path, format="ascii.ecsv")
            for r in v[(v["epoch"] == ep) & (v["verdict"] == "on_sky")]:
                q = e[e["uid"] == r["uid"]]
                snr = np.asarray(q["snr"], float)
                det = np.isfinite(snr) & (snr >= MIN_SNR)
                clean = det & ~np.asarray(q["suspect"], bool) & (np.asarray(q["n_prior"]) > 0)
                key = "one_time" if n_times(np.asarray(q["mjd"], float)[clean]) < 2 else "multi"
                d = res.setdefault(str(r["kind"]), {"one_time": 0, "multi": 0})
                d[key] += 1
        out[fld] = res
    out["meta"] = {"provenance": "derived", "posthoc": True}
    Path(__file__).with_name("posthoc_onsky.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


if __name__ == "__main__":
    print(json.dumps(main(), indent=1))
