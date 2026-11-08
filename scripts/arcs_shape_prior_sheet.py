"""Figure: CEERS F200W shape prior vs the Gaussian σ_ε 0.26 null (results/arcs/shape-prior/).

Reads the catalogue pinned in ``data/manifests/reference_sample.ecsv`` into a temporary directory
and deletes it after reading.
"""

from __future__ import annotations

import shutil
import sys
import tempfile
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import arcs_shape_prior as sp  # noqa: E402

from jwst_anomaly.photometry import fetch_catalog  # noqa: E402

SIGMA = 0.26


def main() -> None:
    tmp = Path(tempfile.mkdtemp(prefix="shape-prior-"))
    try:
        sha = sp.ab.manifest_sha(sp.MANIFEST, sp.PRIOR_FILE)
        shapes = sp.lc.load_shapes(fetch_catalog(sp.ab.s3_url(sp.PRIOR_FILE), sha, cache_dir=tmp))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    fig, ax = plt.subplots(figsize=(6, 4))
    bins = np.linspace(0, 0.95, 39)
    for m in (2.0, 6.0):
        e = sp.prior_eps(shapes, m, 10.0)
        ax.hist(e, bins, density=True, histtype="step", label=f"semimajor ≥ {m:g} px, n {len(e)}")
    x = np.linspace(0, 0.95, 200)
    ax.plot(x, x / SIGMA**2 * np.exp(-(x**2) / 2 / SIGMA**2), "k--", label="Gaussian σ_ε 0.26")
    ax.axvline(1 / 3, color="grey", lw=0.8)  # e = 1 - b/a = 0.5
    ax.set(xlabel="|ε| = (a−b)/(a+b)", ylabel="density", yscale="log", ylim=(0.01, 5))
    ax.legend(fontsize=8)
    ax.set_title("CEERS F200W shape prior vs the Gaussian null")
    fig.tight_layout()
    fig.savefig(sp.OUT / "prior.png", dpi=90)


if __name__ == "__main__":
    main()
