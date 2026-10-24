"""Post hoc (added after gate K failed; no verdict depends on it): gate K numbers of the stored
mean-source residuals without system 10, and system 10's rows. Reads ``imageplane_zspec*.ecsv``
of ``s1063_wslap.py``; writes ``posthoc.json``.

    python scripts/s1063_wslap_posthoc.py
"""

from __future__ import annotations

import json
from pathlib import Path

import m0417_images as mi
import numpy as np
from astropy.table import Table

OUT = Path("results/lens_consistency/s1063_wslap")
DROP = "10"


def run(out: Path = OUT) -> dict:
    ip = Table.read(out / "imageplane_zspec.ecsv")
    ipw = Table.read(out / "imageplane_zspec_wls.ecsv")
    keep = np.asarray(ip["system"], str) != DROP
    rows = ip[~keep]
    res = {
        "label": "post hoc diagnostic, added after gate K failed; no verdict depends on it",
        "gate_K_numbers_mean_source_without_system_10": mi.gate_k(ip[keep]),
        "gate_K_numbers_wls_without_system_10": mi.gate_k(
            ipw[np.asarray(ipw["system"], str) != DROP]
        ),
        "system_10_mean_source": {
            "z": float(rows["z_used"][0]),
            "dtheta_arcsec": [float(v) for v in rows["dtheta_arcsec"]],
            "n_predicted": int(np.max(rows["n_predicted"])),
            "shared_match": [bool(v) for v in rows["shared_match"]],
        },
    }
    (out / "posthoc.json").write_text(json.dumps(res, indent=1) + "\n")
    return res


if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
