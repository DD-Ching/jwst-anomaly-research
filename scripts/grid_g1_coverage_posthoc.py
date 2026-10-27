"""Post hoc split of the MACS J0308 / J2129 grid injections (results/transients/grid_g1_coverage).

Reads each run's ``injections.ecsv`` and ``summary.json`` (``transient_grid.py --coverage``) and
the local mosaics, and writes ``posthoc.json``: recovery at one fraction split by whether the
fraction or the floor set the injected flux and by the smallest epoch-1 exposure count, and per
band the run's floor flux (the smallest injected flux, which the floor sets), the median aperture
errors (epoch 1 on the epoch-2 zero point, as ``band_pass`` does) and the calibrated significance a
floor-sized injection reaches at those medians (no systematic floor). The aperture radius is the
grid default (0.15″) used by both runs. Every number is derived; not in the frozen plan.

    python scripts/grid_g1_coverage_posthoc.py --data data/grid_g1 --runs <dir with m2129/, m0308/>
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import transient_grid as tg  # noqa: E402

RUNS = {
    "m2129": (("F277W", "F444W"), ("jw04903-o004_t004", "jw06882-o040_t043")),
    "m0308": (("F300M", "F410M"), ("jw05293-o005_t001", "jw06882-o015_t015")),
}


def split_recovery(inj: Table, frac: float) -> dict:
    """Recovered fraction by (nominal or floor) × (epoch-1 min exposures 1 or >= 2)."""
    t = inj[np.isclose(inj["frac"], frac)]
    out = {}
    for kind, nom in (("nominal", True), ("floor", False)):
        for cov, sel in (("1", t["e1_n_exp_min"] == 1), (">=2", t["e1_n_exp_min"] >= 2)):
            m = np.asarray(sel & (t["nominal"] == nom))
            rec = np.asarray(t["recovered"])[m]
            out[f"{kind}_e1n_{cov}"] = {
                "n": int(m.sum()),
                "recovered_frac": float(rec.mean()) if m.any() else None,
            }
    return out


def floor_significance(e1: np.ndarray, e2: np.ndarray, noise_scale: float, floor: float) -> float:
    """Calibrated significance of an injected flux ``floor`` at median errors."""
    return float(floor / np.hypot(np.nanmedian(e1), np.nanmedian(e2)) / noise_scale)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--runs", type=Path, required=True)
    ap.add_argument(
        "--out", type=Path, default=Path("results/transients/grid_g1_coverage/posthoc.json")
    )
    a = ap.parse_args(argv)
    res: dict = {"provenance": "derived", "post_hoc": True}
    for name, (bands, stems) in RUNS.items():
        run = a.runs / name
        summary = json.loads((run / "summary.json").read_text())
        grid = Table.read(run / "grid.ecsv")
        ra, dec = np.asarray(grid["ra"]), np.asarray(grid["dec"])
        inj = Table.read(run / "injections.ecsv")
        rec = {"split_frac_1.0": split_recovery(inj, 1.0), "bands": {}}
        for b in bands:
            p1, p2 = (a.data / name / f"{s}_nircam_clear-{b.lower()}_i2d.fits" for s in stems)
            with tg._open(str(p1)) as i1, tg._open(str(p2)) as i2:
                m = tg.measure_band(i1, i2, ra, dec, np.zeros(len(ra), bool), 0.15, None)
            scale = summary["bands"][b]["noise_scale"]
            zp = summary["bands"][b]["zero_point_mag"]
            ok = np.isfinite(m["e1"] + m["e2"])
            e1 = m["e1"][ok] * 10 ** (-0.4 * zp)  # epoch 1 on the epoch-2 scale, as band_pass
            floor = float(np.min(inj[f"{b}_inj_flux"]))  # max(frac * f1, floor) >= floor
            rec["bands"][b] = {
                "floor_flux": floor,
                "median_e1_zp": float(np.median(e1)),
                "median_e2": float(np.median(m["e2"][ok])),
                "noise_scale": scale,
                "floor_sigma_at_median": floor_significance(e1, m["e2"][ok], scale, floor),
            }
        res[name] = rec
    a.out.write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
