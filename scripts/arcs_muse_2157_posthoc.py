"""Post hoc reading of the stored MUSE 2157 z scans (not in the frozen plan): peak contrast."""

import gzip
import json
from pathlib import Path

import numpy as np

OUT = Path("results/arcs/muse-2157")


def peak_contrast(z: np.ndarray, score: np.ndarray, exclude_dz: float = 0.01) -> dict:
    """Best peak, the highest score farther than ``exclude_dz`` from it, and their ratio."""
    k = int(np.argmax(score))
    far = np.abs(z - z[k]) > exclude_dz
    k2 = int(np.flatnonzero(far)[np.argmax(score[far])])
    return {
        "z_best": float(z[k]),
        "score_best": float(score[k]),
        "z_second": float(z[k2]),
        "score_second": float(score[k2]),
        "ratio": float(score[k] / score[k2]),
    }


def main() -> None:
    s = json.load(gzip.open(OUT / "spectra.json.gz", "rt"))
    z = np.arange(0.0, 1.6, 0.0002)
    rows = {key: peak_contrast(z, np.asarray(s[key]["score"])) for key in ("control", "target")}
    out = {"label": "post hoc (not in PLAN.md)", "exclude_dz": 0.01, "rows": rows}
    (OUT / "posthoc.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
