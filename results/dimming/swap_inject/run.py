"""D-039 persistence swap in the injection chain: gates G0 / G1 and swapped limits (PLAN.md).
Usage (repo root): run.py GAIA_ROOT OUT FIELD   (one field per process; then run.py combine OUT)

Writes OUT/<field>/{summary.json, limits_<variant>.ecsv, efficiency_<variant>.ecsv, wht_rel.npz}.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, "scripts")
import dimming_screen as ds  # noqa: E402

COMMITTED = {
    "macs0416": "results/dimming/macs0416_regen",
    "abell2744": "results/dimming/abell2744_calib",
    "nexus": "results/dimming/nexus_calib",
}
WHT = "results/dimming/wht_coverage/runs/{}/limits_{}.ecsv"
VARIANTS = {  # name: (eff_low, eff_high, value for pairs without a committed count)
    "swap": (0.412, 0.961, np.nan),
    "swap_pooled": (0.374, 0.961, np.nan),
    "swap_unknown_low": (0.412, 0.961, 0.374),
}
p = ds.Params()


def gate(a: Table, b: Table, op: str) -> bool:
    x = np.asarray(a["exposure_source_yr"], float)
    y = np.asarray(b["exposure_source_yr"], float)
    tol = 1e-9 * np.maximum(np.abs(x), np.abs(y))
    if op == "eq":
        return bool(np.allclose(x, y, rtol=1e-9, equal_nan=True))
    return bool(np.all(x <= y + tol))


def gain_rho01(lim: Table, base: Table) -> dict:
    """Exposure ratio per t_E at rho = 0.1 (None where the base exposure is 0)."""
    out = {}
    for r, b in zip(lim, base, strict=True):
        if np.isclose(r["rho"], 0.1):
            x, y = float(r["exposure_source_yr"]), float(b["exposure_source_yr"])
            out[f"{r['t_e_yr']:g}"] = round(x / y, 4) if y > 0 else None
    return out


def run_field(gaia_root: Path, out_root: Path, name: str) -> None:
    t0 = time.time()
    out = out_root / name
    out.mkdir(parents=True, exist_ok=True)
    field = ds.load_field(ds.DEFAULT_CONFIG, name)
    gaia = Table.read(gaia_root / name / "gaia.ecsv")
    infl = np.array(
        json.loads(Path(f"{COMMITTED[name]}/inject_summary.json").read_text())["forced_inflation"]
    )
    cats, lc, scale, err, res, ordn = ds._screen(field, p, gaia, out / "wht_rel.npz")
    sel = ds.monitored(lc, err, ordn, p)
    static_veto = ds.catalogue_veto({**ordn, "single_epoch": np.zeros_like(ordn["single_epoch"])})
    covd = np.isfinite(lc["flux"][:, :, 0]).sum(axis=1)
    area = ds.area_deg2(lc["ra"], lc["dec"], covd)
    verdicts = Table.read(f"results/dimming/persistence_run/{name}/verdicts.ecsv")

    def limits(tag, pe):
        eff = ds.efficiency_table(
            lc,
            err,
            static_veto,
            sel,
            ds._times_yr(field),
            ds.T_E_GRID,
            ds.RHO_GRID,
            ds.DEPTHS,
            p,
            seed=0,
            forced_inflation=infl,
            persistence_eff=pe,
        )
        first = eff[
            (eff["model"] == "w3")
            & (eff["t_e_yr"] == eff["t_e_yr"][0])
            & (eff["rho"] == eff["rho"][0])
        ]
        lim = ds.rate_limits(eff, area, int(np.sum(first["n_sources"])))
        lim.meta.update(
            field=name,
            area_deg2=area,
            calibrated=True,
            variant=tag,
            persistence_swap=pe is not None,
        )
        eff.write(out / f"efficiency_{tag}.ecsv", overwrite=True)
        lim.write(out / f"limits_{tag}.ecsv", overwrite=True)
        return lim

    on = limits("veto_on", None)
    ref_on = Table.read(WHT.format(name, "veto_on_on"))
    ref_off = Table.read(WHT.format(name, "veto_off_on"))
    summ = {"field": name, "gate_G0": gate(on, ref_on, "eq")}
    shape = lc["flux"].shape[:2]
    known = np.isfinite(ds.persistence_eff_by_exposures(verdicts, name, shape))
    d = lc["det"][:, :, 0] & sel[:, None]
    summ["monitored_detected_pairs"] = int(d.sum())
    summ["monitored_detected_pairs_with_count"] = int((d & known).sum())
    summ["monitored_detected_pairs_by_epoch_without_count"] = {
        str(k): int((d[:, k] & ~known[:, k]).sum()) for k in range(shape[1])
    }
    if summ["gate_G0"]:
        for tag, (lo, hi, unk) in VARIANTS.items():
            pe = ds.persistence_eff_by_exposures(verdicts, name, shape, eff_low=lo, eff_high=hi)
            pe[~known] = unk
            lim = limits(tag, pe)
            summ[tag] = {
                "gate_G1_lower": gate(on, lim, "le"),
                "gate_G1_upper": gate(lim, ref_off, "le"),
                "gain_vs_veto_on_rho0.1": gain_rho01(lim, on),
            }
        summ["veto_off_gain_rho0.1"] = gain_rho01(ref_off, on)
    summ["wall_s"] = round(time.time() - t0)
    (out / "summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    print(json.dumps(summ), flush=True)


def combine(out_root: Path) -> None:
    summ = {}
    for tag in ("veto_on", *VARIANTS):
        per = {f: Table.read(out_root / f / f"limits_{tag}.ecsv") for f in COMMITTED}
        c = ds.combine_limits(per, calibrated_only=True)
        c.write(out_root / f"limits_combined_{tag}.ecsv", overwrite=True)
        summ[tag] = c
    off = ds.combine_limits(
        {f: Table.read(WHT.format(f, "veto_off_on")) for f in COMMITTED}, calibrated_only=False
    )
    res = {
        "fields": summ["veto_on"].meta["fields"],
        "gate_G0_combined": gate(
            summ["veto_on"], Table.read("results/dimming/wht_coverage/limits_combined.ecsv"), "eq"
        ),
    }
    for tag in VARIANTS:
        res[tag] = {"gain_vs_veto_on_rho0.1": gain_rho01(summ[tag], summ["veto_on"])}
    res["veto_off_gain_rho0.1"] = gain_rho01(off, summ["veto_on"])
    res["unknown_low_over_swap_rho0.1"] = gain_rho01(summ["swap_unknown_low"], summ["swap"])
    (out_root / "summary.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "combine":
        combine(Path(sys.argv[2]))
    else:
        run_field(Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3])
