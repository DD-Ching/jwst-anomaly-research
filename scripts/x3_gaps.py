"""X3 G0: W3 positive controls and false alarms on unselected Gaia epoch photometry (GAPS dev).

Unit file: docs/state/units/X3.md. The Gaia Andromeda Photometric Survey (DR3,
`in_andromeda_survey`) publishes G/BP/RP epoch photometry for all 1,257,319 sources around M31,
with no variability selection. The split is frozen in the unit file: dev = `source_id % 10 == 0`;
nothing else is fetched here.

Subcommands (raw light curves under ``$JWST_ANOMALY_DATA/raw/gaia_dr3_gaps/``, derived tables under
``$JWST_ANOMALY_DATA/derived/x3/``, small tracked results under ``results/x3/``):

- ``ids``: dev source list with G <= ``--gmax`` from the Gaia archive (TAP), cached.
- ``fetch``: a reproducible random dev subsample (``--n``, ``--seed``) of epoch photometry.
- ``g0``: on every fetched light curve, (a) the v0 screen on the light curve as observed (false
  alarms), (b) one W3 injection and (c) one PSPL injection, cells drawn as frozen in X3 G0.
- ``summary``: efficiencies per cell, flag fractions, and the G0 gate verdict.

Screen v0 (frozen in X3 G0): dBIC(N1neg - best ordinary) < ``w3_microlensing.P.flag_dbic`` from the
D-057 fitter, same starts for injected and real light curves. Labels: light curves `observed`,
injections `simulated`, fits and rates `derived`. A flag is never a candidate before vetting.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from astropy.table import Table

from jwst_anomaly import gaia_mulens as gm
from jwst_anomaly import paths, schema

_DIR = Path(__file__).resolve().parent


def _load(name: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, _DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


w3 = _load("w3_microlensing")
wg = _load("w3_gaia")
wg.w3 = w3  # w3_gaia loads its own copy; one class object is needed for pickling to workers

# The frozen dev split is MOD(source_id, 10) = 0. Sync TAP caps results at 50,000 rows and anonymous
# jobs at 60 s, so the GAPS cone (5.5 deg around M31) is queried per HEALPix level-5 pixel through
# its source_id range (source_id = HEALPix level-12 nested index * 2^35), which the archive indexes.
DEV_QUERY = (
    "SELECT source_id, ra, dec, phot_g_mean_mag, bp_rp FROM gaiadr3.gaia_source "
    "WHERE source_id BETWEEN {lo} AND {hi} AND in_andromeda_survey = 'True' "
    "AND MOD(source_id, 10) = 0 AND phot_g_mean_mag <= {gmax}"
)
M31 = (10.6847, 41.2690)
GAPS_RADIUS_DEG = 5.6
HPX_LEVEL = 5
# Injection cells frozen in docs/state/units/X3.md (G0)
TE_CELLS = (10.0, 30.0, 100.0, 300.0)
RHO_CELLS = (0.01, 0.1)
MIN_POINTS = 20  # light curves with fewer accepted G epochs are skipped (ASSUMPTION)


def raw_dir() -> Path:
    d = paths.data_root() / "raw" / "gaia_dr3_gaps"
    d.mkdir(parents=True, exist_ok=True)
    return d


def out_dir() -> Path:
    d = paths.data_root() / "derived" / "x3"
    d.mkdir(parents=True, exist_ok=True)
    return d


def results_dir() -> Path:
    d = paths.repo_root() / "results" / "x3"
    d.mkdir(parents=True, exist_ok=True)
    return d


class GapsDev(gm.GaiaDR3Microlensing):
    """DataLink epoch photometry for GAPS dev sources (reuses the D-061 adapter's fetch)."""

    name = "gaia-dr3-gaps-dev"

    def __init__(self):
        super().__init__(cache=raw_dir())


TAP_ASYNC = "https://gea.esac.esa.int/tap-server/tap/async"


def tap_async_csv(query: str, poll_s: float = 5.0, max_s: float = 1800.0) -> str:
    """Asynchronous TAP job (sync jobs time out on GAPS-wide scans); returns the CSV result."""
    body = urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "PHASE": "RUN", "QUERY": query}
    ).encode()
    with urllib.request.urlopen(TAP_ASYNC, data=body, timeout=120) as resp:  # noqa: S310
        job = resp.geturl()
    t1 = time.time()
    while time.time() - t1 < max_s:
        with urllib.request.urlopen(job + "/phase", timeout=60) as resp:  # noqa: S310
            phase = resp.read().decode().strip()
        if phase == "COMPLETED":
            with urllib.request.urlopen(job + "/results/result", timeout=600) as resp:  # noqa: S310
                return resp.read().decode()
        if phase in ("ERROR", "ABORTED"):
            with urllib.request.urlopen(job, timeout=60) as resp:  # noqa: S310
                msg = ET.tostring(ET.fromstring(resp.read()), encoding="unicode")[-500:]
            raise RuntimeError(f"TAP job {phase}: {msg}")
        time.sleep(poll_s)
    raise TimeoutError(f"TAP job {job} not finished after {max_s} s")


def gaps_pixels() -> list[int]:
    """Level-5 nested HEALPix pixels overlapping the GAPS cone."""
    import astropy.units as u
    from astropy.coordinates import Latitude, Longitude
    from astropy_healpix import HEALPix

    hp = HEALPix(nside=2**HPX_LEVEL, order="nested")
    pix = hp.cone_search_lonlat(
        Longitude(M31[0] * u.deg), Latitude(M31[1] * u.deg), GAPS_RADIUS_DEG * u.deg
    )
    return sorted(int(p) for p in pix)


def ids_table(gmax: float) -> Table:
    path = out_dir() / f"dev_ids_g{gmax:g}.csv"
    if not path.exists():
        shift = 35 + 2 * (12 - HPX_LEVEL)
        parts, header = [], None
        for p in gaps_pixels():
            q = DEV_QUERY.format(lo=p << shift, hi=((p + 1) << shift) - 1, gmax=gmax)
            lines = gm._tap_csv(q).strip().splitlines()
            if len(lines) - 1 >= 50000:
                raise RuntimeError(f"pixel {p} hit the 50,000-row TAP cap")
            header = lines[0]
            parts.extend(lines[1:])
        gm._write_atomic(path, ("\n".join([header, *parts]) + "\n").encode())
    t = Table.read(path, format="ascii.csv")
    t.meta.update(
        provenance=schema.Provenance.OBSERVED.value,
        source="Gaia archive TAP per level-5 HEALPix source_id range: " + DEV_QUERY,
    )
    return t


def subsample(gmax: float, n: int, seed: int) -> Table:
    t = ids_table(gmax)
    rng = np.random.default_rng(seed)
    return t[np.sort(rng.choice(len(t), size=min(n, len(t)), replace=False))]


def light_curve(survey: GapsDev, row) -> w3.LightCurve | None:
    sid = str(row["source_id"])
    try:
        tab = survey.light_curve(sid)
    except KeyError:
        return None
    if len(tab) < MIN_POINTS:
        return None
    return wg.gaia_light_curve(tab, sid, (float(row["ra"]), float(row["dec"])))


def screen(lc: w3.LightCurve) -> dict:
    """v0 screen: the D-057 fitter on one light curve, starts as in D-061's injections."""
    k = np.convolve(lc.f, np.ones(3) / 3, mode="same")
    t0g = float(lc.t[int(np.argmax(k))])
    res = w3.fit_event(lc, t0g, 30.0, 0.3, models=("PSPL", "PAR", "N1neg"), te_grid=(10.0, 100.0))
    best = min(res[m]["bic"] for m in w3.ORDINARY if m in res)
    dbic = float(res["N1neg"]["bic"] - best)
    return {"dbic_N1neg": dbic, "flagged": bool(dbic < w3.P.flag_dbic)}


def _job(args):
    kind, lc, prm = args
    row = {"kind": kind, "event_id": lc.event_id, "n_points": int(len(lc.t)), **prm}
    try:
        if kind != "real":
            base = float(np.median(lc.f))
            fs = prm["fs_frac"] * base
            if kind == "W3":
                lc = w3.inject_w3(lc, base, fs, prm["t0"], prm["tE"], prm["u0"], prm["rho"])
            else:
                lc = w3.inject_w3(lc, base, fs, prm["t0"], prm["tE"], prm["u0"], 0.0, 1.0, 1)
        row.update(screen(lc))
        row["error"] = ""
    except Exception as exc:  # noqa: BLE001
        row.update({"dbic_N1neg": np.nan, "flagged": False, "error": repr(exc)[:200]})
    return row


def run_g0(gmax: float, n: int, seed: int, procs: int) -> Path:
    survey = GapsDev()
    sample = subsample(gmax, n, seed)
    rng = np.random.default_rng(seed + 1)
    jobs = []
    for row in sample:
        lc = light_curve(survey, row)
        if lc is None:
            continue
        base = {
            "G_mean": float(row["phot_g_mean_mag"]),
            "t_first": float(lc.t.min()),
            "t_last": float(lc.t.max()),
        }
        jobs.append(
            (
                "real",
                lc,
                {
                    **base,
                    "tE": np.nan,
                    "rho": np.nan,
                    "u0": np.nan,
                    "t0": np.nan,
                    "fs_frac": np.nan,
                },
            )
        )
        for kind in ("W3", "PSPL"):
            prm = {
                **base,
                "tE": float(TE_CELLS[int(rng.integers(len(TE_CELLS)))]),
                "rho": float(RHO_CELLS[int(rng.integers(len(RHO_CELLS)))]) if kind == "W3" else 0.0,
                "u0": float(rng.uniform(0.0, 2.0 if kind == "W3" else 1.0)),
                "t0": float(rng.uniform(lc.t.min(), lc.t.max())),
                "fs_frac": float(rng.uniform(0.5, 1.0)),
            }
            jobs.append((kind, lc, prm))
    t1 = time.time()
    with Pool(procs) as pool:
        rows = list(pool.imap_unordered(_job, jobs, 4))
    tab = Table(rows=rows, names=list(rows[0].keys()))
    tab.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source="scripts/x3_gaps.py g0: v0 screen (D-057 fitter) on GAPS dev light curves "
        "(observed) and on W3/PSPL injections (simulated) into them",
        gmax=gmax,
        n_requested=n,
        seed=seed,
        flag_dbic=w3.P.flag_dbic,
        wall_time_s=round(time.time() - t1, 1),
    )
    path = out_dir() / f"g0_seed{seed}.ecsv"
    tab.write(path, overwrite=True)
    w3.write_ecsv_gz(tab, results_dir() / f"g0_seed{seed}.ecsv.gz")
    print(f"wrote {path}: {len(tab)} rows in {time.time() - t1:.0f} s")
    return path


def _col(tab: Table, name: str, fill):
    c = tab[name]
    return np.asarray(c.filled(fill) if hasattr(c, "filled") else c)


def summarise(tab: Table) -> dict:
    tab = tab.copy()
    for name, fill in (("error", ""), ("flagged", False), ("kind", "")):
        tab[name] = _col(tab, name, fill)
    ok = tab["error"] == ""
    out = {"n_rows": int(len(tab)), "n_errors": int((~ok).sum())}
    real = ok & (tab["kind"] == "real")
    out["n_real"] = int(real.sum())
    out["flag_fraction_real"] = float(np.mean(tab["flagged"][real])) if real.any() else None
    pspl = ok & (tab["kind"] == "PSPL")
    out["pspl_flagged_as_w3"] = float(np.mean(tab["flagged"][pspl])) if pspl.any() else None
    cells = {}
    for te in TE_CELLS:
        for rho in RHO_CELLS:
            m = ok & (tab["kind"] == "W3") & (tab["tE"] == te) & (tab["rho"] == rho)
            cells[f"tE{te:g}_rho{rho:g}"] = {
                "n": int(m.sum()),
                "eff": float(np.mean(tab["flagged"][m])) if m.any() else None,
            }
    out["w3_cells"] = cells
    w = ok & (tab["kind"] == "W3")
    effs = [
        float(np.mean(tab["flagged"][w & (tab["tE"] == te)]))
        for te in (30.0, 100.0, 300.0)
        if (w & (tab["tE"] == te)).any()
    ]
    out["mean_eff_te30_300"] = float(np.mean(effs)) if effs else None
    e_all = [c["eff"] for c in cells.values() if c["eff"] is not None]
    if out["mean_eff_te30_300"] is None or out["flag_fraction_real"] is None:
        verdict = "incomplete"
    elif out["mean_eff_te30_300"] >= 0.20 and out["flag_fraction_real"] <= 0.01:
        verdict = "continue to G1"
    elif e_all and max(e_all) < 0.05:
        verdict = "stop: Gaia cadence cannot test W3"
    else:
        verdict = "one recorded screen revision allowed, re-gate on fresh dev"
    out["gate"] = verdict
    out["gate_definition"] = "docs/state/units/X3.md G0 (frozen 2026-10-14 before data)"
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("cmd", choices=("ids", "fetch", "g0", "summary"))
    ap.add_argument("--gmax", type=float, default=19.0)
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=314)
    ap.add_argument("--procs", type=int, default=4)
    a = ap.parse_args(argv)
    if a.cmd == "ids":
        t = ids_table(a.gmax)
        print(f"dev sources with G <= {a.gmax}: {len(t)}")
    elif a.cmd == "fetch":
        ids = [str(s) for s in subsample(a.gmax, a.n, a.seed)["source_id"]]
        t1 = time.time()
        n = GapsDev().prefetch(ids)
        print(f"fetched {n} new light curves of {len(ids)} in {time.time() - t1:.0f} s")
    elif a.cmd == "g0":
        run_g0(a.gmax, a.n, a.seed, a.procs)
    else:
        tab = Table.read(out_dir() / f"g0_seed{a.seed}.ecsv")
        s = summarise(tab)
        (results_dir() / f"g0_summary_seed{a.seed}.json").write_text(json.dumps(s, indent=1) + "\n")
        print(json.dumps(s, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
