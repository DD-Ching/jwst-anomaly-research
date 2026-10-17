"""W1 shear null as a comoving volume density limit over lens redshift (results/w1_volume/PLAN.md).

Reads ``<runs>/z<Z>/<field>/summary.json`` written by ``inject_shear.py --lens-z Z --out <runs>/z<Z>``
and, per mass, integrates the efficiency over the comoving volume behind each field's screened area:

    V_eff = sum_f A_f int eps_f(z) dV_c/(dz dOmega) dz,    n_95 = -ln(0.05) / V_eff

with eps_f linear between the grid redshifts and zero outside them (ASSUMPTION, conservative).
Zero detections (the screen is null in every field). Efficiencies ``derived`` from ``simulated``
injections; the limit is ``derived``; the comparison value (Takahashi & Asada 2013) is a published
limit.

    python scripts/w1_volume.py --runs OUT --fields abell2744 ... --out results/w1_volume/limits.json
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from astropy.cosmology import Planck18

from jwst_anomaly import schema

TA13_H3_MPC3 = 1e-4  # Takahashi & Asada 2013 (arXiv:1303.1301), |M| > 1e12 M_sun, h^3 Mpc^-3
N_FINE = 1401  # integration points between the grid ends


def v_eff_mpc3(z_grid, eff, area_deg2, cosmology=Planck18) -> float:
    """Efficiency-weighted comoving volume (Mpc^3): ``eff`` is (fields, len(z_grid)), linear in z
    between grid points, zero outside; ``area_deg2`` per field."""
    z_grid = np.asarray(z_grid, float)
    eff = np.atleast_2d(np.asarray(eff, float))
    if np.any(np.diff(z_grid) <= 0):
        raise ValueError("z_grid must increase")
    z = np.linspace(z_grid[0], z_grid[-1], N_FINE)
    dvdz = cosmology.differential_comoving_volume(z).to_value("Mpc3 / sr")
    sr = np.asarray(area_deg2, float) * (math.pi / 180.0) ** 2
    weight = sum(a * np.interp(z, z_grid, e) for a, e in zip(sr, eff, strict=True))
    return float(np.trapezoid(weight * dvdz, z))


def limit(v_eff: float, h: float = Planck18.h) -> dict:
    """95 % upper limit on n for zero detections (Poisson)."""
    if v_eff <= 0:
        return {"v_eff_mpc3": v_eff, "n95_mpc3": None, "n95_h3_mpc3": None, "ratio_to_ta13": None}
    n = -math.log(0.05) / v_eff
    n_h = n / h**3
    return {
        "v_eff_mpc3": v_eff,
        "n95_mpc3": n,
        "n95_h3_mpc3": n_h,
        "ratio_to_ta13": n_h / TA13_H3_MPC3,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--runs", type=Path, required=True)
    ap.add_argument("--fields", nargs="+", required=True)
    ap.add_argument("--z", nargs="+", type=float, default=[0.1, 0.3, 0.6, 1.0, 1.5])
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    summ = {
        (z, f): json.loads((args.runs / f"z{z:g}" / f / "summary.json").read_text())
        for z in args.z
        for f in args.fields
    }
    for (z, f), s in summ.items():
        if s["assumptions"].get("w1_lens_z") != z:
            raise SystemExit(f"error: {f} at z {z}: summary has w1_lens_z {s['assumptions']}")
    masses = list(summ[(args.z[0], args.fields[0])]["efficiency"])
    area = [summ[(args.z[0], f)]["screened_area_deg2"] for f in args.fields]
    out = {
        "fields": args.fields,
        "z_grid": args.z,
        "total_area_deg2": float(sum(area)),
        "cosmology": "Planck18",
        "comparison": "Takahashi & Asada 2013 (arXiv:1303.1301): n < 1e-4 h^3 Mpc^-3, |M| > 1e12",
        "provenance": schema.Provenance.DERIVED.value,
        "masses": {},
    }
    for m in masses:
        eff = [[summ[(z, f)]["efficiency"][m]["efficiency"] for z in args.z] for f in args.fields]
        n_rec = [[summ[(z, f)]["efficiency"][m]["n_recovered"] for z in args.z] for f in args.fields]
        out["masses"][m] = limit(v_eff_mpc3(args.z, eff, area)) | {
            "n_recovered_by_z": dict(zip([f"{z:g}" for z in args.z], np.sum(n_rec, 0).tolist())),
            "n_inject_per_field_and_z": summ[(args.z[0], args.fields[0])]["efficiency"][m][
                "n_inject"
            ],
            "theta_e_zs2_arcsec_by_z": {
                f"{z:g}": summ[(z, args.fields[0])]["efficiency"][m]["theta_e_zs2_arcsec"]
                for z in args.z
            },
        }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
