"""Gate K of the centre-level axial rule for ``radial`` (results/radial_axial/PLAN.md).

Reads a ``radial_more.py`` output (``<run>/<field>/<model>/radial.json`` and
``radial_sources.ecsv``), and for every convergence peak computes the contributing arcs, their
azimuths about the centre and the axial concentration A (``exotic_screens.axial_concentration``).
Gate K passes when the field's strongest centre has A >= ``exotic_screens.AXIAL_MIN``.

Also writes a post-hoc diagnostic (not part of the frozen rule): ``n_axes``, the number of
distinct axes among the contributors when azimuths within ``AXIS_TOL_DEG`` (mod 180) are merged,
and ``max_on_one_line``, the most contributors within ``LINE_TOL_ARCSEC`` of one straight line
through two of them (collinearity of the arcs themselves, independent of the grid point).

    python scripts/radial_axial.py --run RUN_DIR --field macs0647 --out results/radial_axial
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))

import exotic_screens as es  # noqa: E402
import inject_shear as ish  # noqa: E402

AXIS_TOL_DEG = 10.0  # post-hoc diagnostic only (ASSUMPTION)
LINE_TOL_ARCSEC = 0.5  # post-hoc diagnostic only (ASSUMPTION)


def n_axes(phi_deg, tol: float = AXIS_TOL_DEG) -> int:
    """Distinct axes among azimuths (deg, mod 180): greedy single-linkage on the circle."""
    a = np.sort(np.mod(np.asarray(phi_deg, float), 180.0))
    if a.size == 0:
        return 0
    gaps = np.diff(np.concatenate([a, [a[0] + 180.0]]))
    return max(1, int(np.sum(gaps > tol)))


def on_line(x, y, tol: float = LINE_TOL_ARCSEC) -> int:
    """Most contributors within ``tol`` of one straight line through two of them (post hoc)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    best = min(len(x), 2)
    for i in range(len(x)):
        for j in range(i + 1, len(x)):
            dx, dy = x[j] - x[i], y[j] - y[i]
            norm = np.hypot(dx, dy)
            if norm == 0:
                continue
            perp = np.abs((x - x[i]) * dy - (y - y[i]) * dx) / norm
            best = max(best, int(np.sum(perp <= tol)))
    return best


def gate_k(peaks: list[dict]) -> bool:
    """The frozen gate: the strongest centre is axial (A >= AXIAL_MIN)."""
    top = max(peaks, key=lambda p: p["n_lines"])
    return bool(top["axial_a"] >= es.AXIAL_MIN)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", type=Path, required=True)
    ap.add_argument("--field", required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    spec = ish.field_spec(a.field)
    d = a.run / a.field / spec["model"]
    summ = json.loads((d / "radial.json").read_text())
    cand = Table.read(d / "radial_sources.ecsv")
    model, _, _ = es.lc.load_model(spec["model"])
    es.lc.apply_frame_offset(spec["model"], model)
    x, y = (np.asarray(v, float) for v in model.to_frame(cand["ra"], cand["dec"]))
    pa = np.asarray(cand["pa_obs"], float)
    args = es.radial_defaults()
    tol, max_len = args["line_tol_arcsec"], args["max_len_arcsec"]
    peaks = []
    for p in summ["peaks"]:
        n, a_ax = es.axial_concentration(x, y, pa, p["x"], p["y"], tol, max_len)
        gx, gy = np.array([[p["x"]]]), np.array([[p["y"]]])
        hit = np.array(
            [
                bool(es.line_counts([xi], [yi], [pi], gx, gy, tol, max_len)[0, 0])
                for xi, yi, pi in zip(x, y, pa, strict=True)
            ]
        )
        assert n == hit.sum() == p["n_lines"]
        phi = np.rad2deg(np.arctan2(y[hit] - p["y"], x[hit] - p["x"]))
        peaks.append(
            {
                "x": p["x"],
                "y": p["y"],
                "n_lines": p["n_lines"],
                "p_random": p["p_random"],
                "axial_a": a_ax,
                "axial": bool(a_ax >= es.AXIAL_MIN),
                "labels": [int(v) for v in np.asarray(cand["label"])[hit]],
                "azimuth_deg": [round(float(v), 1) for v in phi],
                "post_hoc_n_axes": n_axes(phi),
                "post_hoc_max_on_one_line": on_line(x[hit], y[hit]),
            }
        )
    out = {
        "field": a.field,
        "model": spec["model"],
        "max_lines": summ["max_lines"],
        "p_random_max": summ["p_random_max"],
        "axial_min": es.AXIAL_MIN,
        "gate_k_pass": gate_k(peaks),
        "post_hoc_axis_tol_deg": AXIS_TOL_DEG,
        "post_hoc_line_tol_arcsec": LINE_TOL_ARCSEC,
        "peaks": peaks,
        "provenance": "derived (shapes observed; azimuths and A derived; thresholds ASSUMPTION)",
    }
    a.out.mkdir(parents=True, exist_ok=True)
    path = a.out / f"gate_k_{a.field}.json"
    path.write_text(json.dumps(out, indent=1) + "\n")
    print(
        json.dumps({k: out[k] for k in ("field", "gate_k_pass")}),
        [
            (
                q["n_lines"],
                round(q["axial_a"], 3),
                q["post_hoc_n_axes"],
                q["post_hoc_max_on_one_line"],
            )
            for q in peaks
        ],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
