"""W3 in the MOA-II 9-year release, pilot field gb22: pre-screen, fits, vetting, limit (D-062).

The release holds every Cut-0 variable object (difference-image detections of positive *or
negative* PSF profiles; ``jwst_anomaly.moa``), before any bump or PSPL cut, so a W3 event (the
source flux drops toward zero inside an umbra between two caustic spikes; ``exotic_sim``) can be in
it. Fitting all 18,599 gb22 light curves with every model is too slow, so:

- ``prescreen``: a W3-shaped matched statistic on every light curve (``deficit_scan``): the most
  significant box of width W = 1 … 300 d that is below both its flanks and the median flux, in
  units of the light curve's own spread of that statistic (red noise), and no second such deficit
  elsewhere (``prescreen_pass``). Thresholds are set from the real distribution and the injections.
- ``fit``: the passes are fitted with the ordinary (PSPL, FSPL, PAR) and exotic (N1neg, E2pos,
  E2neg) models of ``scripts/w3_microlensing.py`` on one trajectory, blend flux free (difference
  flux), extra exotic starts on the deficit. Flag: ΔBIC < −10 (ASSUMPTION, as D-057).
- ``vet``: flags through ordinary explanations, cheapest first (``vet_one``).
- ``inject``: W3 events (``exotic_sim``, ``simulated``) added to real gb22 light curves of quiet
  objects, through a light-curve-level Cut-0 emulation, the pre-screen, the fit and the vetting.
- ``limit``: 95 % upper limit on the W3 rate per monitored star per year for gb22.

Outputs go to ``$JWST_ANOMALY_DATA/derived/w3_moa/`` (never in git). Exotic physics is a
hypothesis: a flag is an anomaly to vet, never a discovery.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from dataclasses import asdict, dataclass, replace
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack

sys.path.insert(0, str(Path(__file__).resolve().parent))

import w3_microlensing as w3  # noqa: E402

from jwst_anomaly import exotic_sim as es  # noqa: E402
from jwst_anomaly import moa, paths, schema  # noqa: E402

FIELD = "gb22"
REPEAT_GAP = 4.0  # ASSUMPTION: boxes this many widths from the deficit count as a second dip
SCAN_KEYS = ("z_min", "s_min", "z_min2", "spread", "width", "t_lo", "t_hi", "err_scale")


@dataclass(frozen=True)
class Params:
    """Every threshold is an ASSUMPTION unless marked as a published cut."""

    widths: tuple = (1.0, 3.0, 10.0, 30.0, 100.0, 300.0)  # deficit boxes (days)
    n_min: int = 3  # epochs in a box
    prescreen_z: float = 10.0  # pass: z_min < −prescreen_z (self-normalised notch) …
    prescreen_s: float = 5.0  # … and S_min < −prescreen_s (white-noise S/N of the same box) …
    prescreen_repeat: float = 8.0  # … and no second deficit z_min2 < −prescreen_repeat
    flag_dbic: float = -10.0  # as D-057
    # Cut-0 (Koshimoto et al. 2023, Table 2; published): S/N > 2.7 on the difference image,
    # N_continue,8 ≥ 3 detections each within 8 days of the previous one
    cut0_snr: float = 2.7
    cut0_n: int = 3
    cut0_gap_days: float = 8.0
    max_points_fit: int = 6000


P = Params()


def out_dir() -> Path:
    d = paths.data_root() / "derived" / "w3_moa"
    d.mkdir(parents=True, exist_ok=True)
    return d


# ----------------------------------------------------------------------------- pre-screen


def robust_scale(f: np.ndarray, sf: np.ndarray) -> tuple[float, float]:
    """Median flux and the error scale from point-to-point scatter (1.4826 MAD of successive
    differences of f/σ, / √2, ≥ 1): insensitive to slow trends and to the event itself."""
    b = float(np.median(f))
    if f.size < 3:
        return b, 1.0
    d = np.diff(f) / np.hypot(sf[1:], sf[:-1])
    s = 1.4826 * float(np.median(np.abs(d - np.median(d))))
    return b, max(s, 1.0)


def nights(t: np.ndarray) -> np.ndarray:
    """Night number (MOA, New Zealand: nights split at 00:00 UT = JD fraction 0.5)."""
    return np.floor(np.asarray(t) - 0.5).astype(np.int64)


def spread_of(x: np.ndarray, valid: np.ndarray) -> float:
    """Robust spread (1.4826 MAD, ≥ 1) of a box statistic over the valid boxes."""
    if valid.sum() < 10:
        return 1.0
    v = x[valid]
    return max(1.4826 * float(np.median(np.abs(v - np.median(v)))), 1.0)


def deficit_scan(t, f, sf, widths=None, n_min=None) -> dict:
    """Most significant local flux deficit (``derived``): a notch filter.

    For every box [t_i, t_i + W) the weighted mean flux is compared with the weighted mean of the
    flanks [t_i − W, t_i) ∪ [t_i + W, t_i + 2W): S = (m_box − m_flank) / σ, errors multiplied by
    the point-to-point scale. Slow trends cancel; a W3 umbra (deficit between two brightenings)
    gives S ≪ 0. Boxes need ``n_min`` epochs on ≥ 2 nights, flanks ``n_min`` epochs. Returns S_min
    with its box, the most significant deficit in a box apart from it (``s_min2``: repeated dips
    are eclipses or pulsations) and the most significant excess S_max.
    """
    widths = P.widths if widths is None else widths
    n_min = P.n_min if n_min is None else n_min
    t = np.asarray(t, float)
    f = np.asarray(f, float)
    sf = np.asarray(sf, float)
    out = {"z_min": 0.0, "s_min": 0.0, "z_max": 0.0, "z_min2": 0.0, "spread": 1.0,
           "width": np.nan, "t_lo": np.nan, "t_hi": np.nan, "n_box": 0, "depth": np.nan,
           "depth_err": np.nan}  # fmt: skip
    if t.size < 2 * n_min:
        out.update(base=float(np.median(f)) if f.size else np.nan, err_scale=1.0, n_points=t.size)
        return out
    b, s = robust_scale(f, sf)
    w = 1.0 / (sf * s) ** 2
    r = f - b
    cw = np.concatenate([[0.0], np.cumsum(w)])
    cwr = np.concatenate([[0.0], np.cumsum(w * r)])
    nt = nights(t)
    cnew = np.concatenate([[0], np.cumsum(np.diff(nt) != 0)])  # night changes up to each epoch
    i = np.arange(t.size)
    per_width = []
    for width in widths:
        j = np.searchsorted(t, t + width, side="left")
        a = np.searchsorted(t, t - width, side="left")
        c = np.searchsorted(t, t + 2 * width, side="left")
        n_in = j - i
        n_fl = (i - a) + (c - j)
        jm = np.maximum(j - 1, i)
        n_nights = 1 + cnew[jm] - cnew[i]
        sw_in = cw[j] - cw[i]
        sw_fl = (cw[i] - cw[a]) + (cw[c] - cw[j])
        with np.errstate(invalid="ignore", divide="ignore"):
            m_in = (cwr[j] - cwr[i]) / sw_in
            m_fl = ((cwr[i] - cwr[a]) + (cwr[c] - cwr[j])) / sw_fl
            st = (m_in - m_fl) / np.sqrt(1.0 / sw_in + 1.0 / sw_fl)
        valid = (n_in >= n_min) & (n_nights >= 2) & (n_fl >= n_min) & np.isfinite(st)
        st = np.where(valid, st, 0.0)
        with np.errstate(invalid="ignore", divide="ignore"):
            sg = np.where(valid, m_in * np.sqrt(sw_in), 0.0)  # box mean vs the median flux

        # red noise: S of real light curves is wider than N(0, 1); normalise by its own spread
        spread = spread_of(st, valid)
        # below its flanks AND below the median flux: the weaker of the two significances
        z = np.maximum(st / spread, sg / spread_of(sg, valid))
        k = int(np.argmin(z))
        per_width.append((width, z, j))
        if z[k] < out["z_min"]:
            out.update(
                z_min=float(z[k]),
                spread=spread,
                s_min=float(min(st[k], sg[k])),
                width=float(width),
                t_lo=float(t[k]),
                t_hi=float(t[jm[k]]),
                n_box=int(n_in[k]),
                depth=float(m_in[k] - m_fl[k]),
                depth_err=float(np.sqrt(1.0 / sw_in[k] + 1.0 / sw_fl[k])),
            )
        out["z_max"] = max(out["z_max"], float(z.max()))
    if np.isfinite(out["width"]):
        lo, hi, wd = out["t_lo"], out["t_hi"], out["width"]
        for width, z, j in per_width:
            t_end = t[np.maximum(j - 1, 0)]
            gap = REPEAT_GAP * max(wd, width)  # a long umbra spans several boxes of one width
            apart = (t_end < lo - gap) | (t > hi + gap)
            if apart.any():
                out["z_min2"] = min(out["z_min2"], float(z[apart].min()))
    out.update(base=b, err_scale=s, n_points=int(t.size))
    return out


def cut0_emulated(t, signal, sf, snr=None, n_req=None, gap=None) -> bool:
    """Light-curve-level Cut-0: ≥ ``n_req`` epochs with |signal| / σ > ``snr`` (either sign), each
    within ``gap`` days of the previous one. ``signal`` is the injected difference flux without
    noise: the emulation asks whether DAOFIND could have seen the event itself (ASSUMPTION: the
    light-curve S/N equals the image S/N; spurious-detection filters are not emulated)."""
    snr = P.cut0_snr if snr is None else snr
    n_req = P.cut0_n if n_req is None else n_req
    gap = P.cut0_gap_days if gap is None else gap
    det = np.asarray(t)[np.abs(np.asarray(signal)) / np.asarray(sf) > snr]
    if det.size < n_req:
        return False
    run = best = 1
    for d in np.diff(det):
        run = run + 1 if d <= gap else 1
        best = max(best, run)
    return best >= n_req


def _prescreen_worker(job):
    ids, tar_path = job
    field = moa.MoaField(FIELD, tar_path=Path(tar_path))
    rows = []
    with open(tar_path, "rb") as fh:
        idx = field.index()
        for eid in ids:
            try:
                field._index = idx
                cols = field.read_member(eid, fh)
                t, f, sf, kind = moa.select_flux(cols)
                ok = np.isfinite(f) & np.isfinite(sf) & (sf > 0)
                t, f, sf = t[ok], f[ok], sf[ok]
                row = {"event_id": eid, "flux_kind": kind, **deficit_scan(t, f, sf)}
                row["error"] = ""
            except Exception as exc:  # noqa: BLE001 — one bad file must not stop the pass
                row = {"event_id": eid, "error": repr(exc)[:200]}
            rows.append(row)
    return rows


def rows_to_table(rows: list[dict], fill=None) -> Table:
    fill = fill or {}
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    cols = {}
    for k in keys:
        vals = [r.get(k, fill.get(k, np.nan)) for r in rows]
        if any(isinstance(v, str) for v in vals):
            vals = [
                "" if (v is None or (isinstance(v, float) and np.isnan(v))) else str(v)
                for v in vals
            ]
        cols[k] = vals
    return Table(cols)


def run_prescreen(procs: int) -> Path:
    field = moa.MoaField(FIELD)
    idx = field.index()
    ids = [k for k, _ in sorted(idx.items(), key=lambda kv: kv[1][0])]
    tar = str(field.tar_path())
    chunks = [(ids[i : i + 250], tar) for i in range(0, len(ids), 250)]
    t1 = time.time()
    rows = []
    with Pool(procs) as pool:
        for part in pool.imap_unordered(_prescreen_worker, chunks):
            rows.extend(part)
    tab = rows_to_table(rows, {"error": "", "flux_kind": ""})
    tab.sort("event_id")
    tab.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source=f"scripts/w3_moa.py prescreen: deficit_scan on every {field.name} light curve",
        params=json.dumps(asdict(P)),
        wall_time_s=round(time.time() - t1, 1),
    )
    path = out_dir() / f"prescreen_{FIELD}.ecsv"
    tab.write(path, overwrite=True)
    ok = np.asarray(tab["error"], str) == ""
    s = np.asarray(tab["s_min"])[ok]
    z = np.asarray(tab["z_min"])[ok]
    print(f"wrote {path}: {len(tab)} light curves ({(~ok).sum()} errors), {time.time() - t1:.0f} s")
    for thr in (5, 6, 8, 10, 12, 15, 20, 30):
        print(
            f"  z_min < -{thr}: {(z < -thr).sum()} (and S_min < -{P.prescreen_s}: "
            f"{((z < -thr) & (s < -P.prescreen_s)).sum()})"
        )
    print(f"passes: {len(passes(tab))}")
    return path


# ----------------------------------------------------------------------------- fit


def to_lightcurve(t, f, sf, ra=None, dec=None, event_id="", err_scale=1.0) -> w3.LightCurve:
    """Fitter light curve in difference-flux counts: blend flux unbounded (``f_min = inf``),
    errors multiplied by the robust scale of the pre-screen (``err_scale``; ASSUMPTION)."""
    return w3.LightCurve(t, f, np.asarray(sf) * err_scale, ra, dec, event_id, f_min=np.inf)


def deficit_starts(model: str, t_lo: float, t_hi: float) -> list:
    """Starts that put the umbra of a repulsive lens (n, ε < 0) on the pre-screen deficit."""
    n, sign = w3.EXOTIC[model]
    if sign != -1:
        return []
    bc = es.caustic_beta(n)
    tc = 0.5 * (t_lo + t_hi)
    half = max(0.5 * (t_hi - t_lo), 0.3)
    out = []
    for grow in (1.0, 1.5, 3.0):
        for u0 in bc * np.array([0.2, 0.6, 0.9, 0.98]):
            te = grow * half / math.sqrt(bc * bc - u0 * u0)
            if w3.P.te_bounds[0] < te < w3.P.te_bounds[1]:
                for lr in (-2.5, -1.2):
                    out.append((tc, math.log10(te), u0, lr))
    return out


def fit_moa(lc: w3.LightCurve, scan: dict, models=None) -> dict:
    """``w3_microlensing.fit_event`` (PSPL bump started on the brightest 3-epoch mean) plus
    repulsive-lens starts on the deficit found by the pre-screen; the better optimum is kept."""
    models = models or ("PSPL", "FSPL", "PAR", *w3.EXOTIC)
    k = np.convolve(lc.f, np.ones(3) / 3, mode="same")
    t0g = float(lc.t[int(np.argmax(k))])
    res = w3.fit_event(lc, t0g, 10.0, 0.3, models=models, te_grid=(3.0, 30.0, 100.0))
    if np.isfinite(scan.get("t_lo", np.nan)):
        for m in w3.EXOTIC:
            if m not in models:
                continue
            st = deficit_starts(m, scan["t_lo"], scan["t_hi"])
            if st:
                r = w3.optimise(m, lc, st, n_best=3)
                if r["bic"] < res[m]["bic"]:
                    res[m] = r
    return res


def _fit_worker(job):
    ev, t, f, sf, scan = job
    t1 = time.time()
    lc = to_lightcurve(t, f, sf, ev["ra"], ev["dec"], ev["event_id"], scan["err_scale"])
    try:
        res = fit_moa(lc, scan, models=SCREEN_MODELS)
        row = w3.summarise(ev["event_id"], res, lc.t.size)
        row["error"] = ""
    except Exception as exc:  # noqa: BLE001
        row = {"event_id": ev["event_id"], "n_points": lc.t.size, "error": repr(exc)[:200]}
    row.update({k: scan[k] for k in SCAN_KEYS})
    row["seconds"] = time.time() - t1
    return row


AUX = ("fwhm", "airmass", "sky")  # observing conditions per epoch (vetting regressors)


def load_arrays(field: moa.MoaField, ids, aux: bool = False) -> dict:
    """``event_id -> (t, f, sf)`` (plus a dict of ``AUX`` columns when ``aux``) for the requested
    light curves, included epochs only (one pass over the tar)."""
    out = {}
    for eid, cols in field.iter_light_curves(ids):
        t, f, sf, _kind = moa.select_flux(cols)
        inc = cols["included"] & np.isfinite(cols["flux"]) & np.isfinite(cols["flux_err"])
        ok = np.isfinite(f) & np.isfinite(sf) & (sf > 0)
        if aux:
            ax = {k: cols[k][inc][ok] for k in AUX if k in cols}
            out[eid] = (t[ok], f[ok], sf[ok], ax)
        else:
            out[eid] = (t[ok], f[ok], sf[ok])
    return out


def no_error(tab: Table) -> np.ndarray:
    """Rows without an error (ECSV reads empty strings back as masked)."""
    col = tab["error"]
    vals = col.filled("") if hasattr(col, "filled") else col
    return np.asarray(vals, str) == ""


def prescreen_pass(z_min, s_min, z_min2):
    """The pre-screen: one significant local deficit (z and white-noise S) and no second one."""
    z, s, z2 = (np.asarray(x, float) for x in (z_min, s_min, z_min2))
    return (z < -P.prescreen_z) & (s < -P.prescreen_s) & (z2 > -P.prescreen_repeat)


def passes(pre: Table) -> Table:
    return pre[no_error(pre) & prescreen_pass(pre["z_min"], pre["s_min"], pre["z_min2"])]


def results_dir() -> Path:
    """Tracked chunk fit tables (as D-059): an ephemeral session fits one chunk, `merge-chunks`
    joins them into the table `vet` reads."""
    d = paths.repo_root() / "results" / "w3_moa"
    d.mkdir(parents=True, exist_ok=True)
    return d


def chunk_name(chunk: tuple[int, int]) -> str:
    return f"fits_{FIELD}_chunk{chunk[0] + 1}of{chunk[1]}.ecsv"


def run_fit(procs: int, limit: int | None = None, chunk: tuple[int, int] | None = None) -> Path:
    """Fit the pre-screen passes; ``chunk=(k, n)`` fits every n-th pass from k (0-based)."""
    field = moa.MoaField(FIELD)
    pre = Table.read(out_dir() / f"prescreen_{FIELD}.ecsv")
    sel = passes(pre)
    if chunk is not None:
        sel = sel[chunk[0] :: chunk[1]]
    if limit:
        sel = sel[:limit]
    ev = field.events()
    pos = {r["event_id"]: (float(r["ra"]), float(r["dec"])) for r in ev}
    arrays = load_arrays(field, sel["event_id"])
    jobs = []
    for r in sel:
        eid = str(r["event_id"])
        t, f, sf = arrays[eid]
        scan = {k: float(r[k]) for k in SCAN_KEYS}
        d = {"event_id": eid, "ra": pos[eid][0], "dec": pos[eid][1]}
        jobs.append((d, t, f, sf, scan))
    print(f"fitting {len(jobs)} pre-screen passes of {len(pre)}", flush=True)
    t1 = time.time()
    rows = []
    with Pool(procs) as pool:
        for i, row in enumerate(pool.imap_unordered(_fit_worker, jobs, 1)):
            rows.append(row)
            if (i + 1) % 25 == 0:
                print(f"{i + 1}/{len(jobs)} fitted, {time.time() - t1:.0f} s", flush=True)
    tab = rows_to_table(rows, {"error": ""})
    tab.sort("event_id")
    tab.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source=f"scripts/w3_moa.py fit: w3_microlensing models on {field.name} pre-screen passes",
        params=json.dumps(asdict(P)),
        fit_params=json.dumps(asdict(w3.P)),
        wall_time_s=round(time.time() - t1, 1),
        n_prescreen=len(pre),
        n_passes=len(passes(pre)),
        chunk="" if chunk is None else f"{chunk[0] + 1}/{chunk[1]}",
    )
    path = out_dir() / f"fits_{FIELD}.ecsv"
    tab.write(path, overwrite=True)
    if chunk is not None and limit is None:  # complete chunks only: what `merge-chunks` joins
        w3.write_ecsv_gz(tab, results_dir() / f"{chunk_name(chunk)}.gz")
    ok = no_error(tab)
    flags = ok & (np.asarray(tab["dbic_min"], float) < P.flag_dbic)
    print(
        f"wrote {path}: {len(tab)} fits, {(~ok).sum()} errors, {flags.sum()} flags, "
        f"{time.time() - t1:.0f} s"
    )
    return path


def merge_chunks(n: int) -> Path:
    """Join the tracked chunk tables 1..n of n into the table `vet` reads. Refused unless every
    chunk is present, was fitted with the current fit ``Params`` and holds exactly its own passes
    of the current pre-screen (which is deterministic and recomputed in each session)."""
    pre = Table.read(out_dir() / f"prescreen_{FIELD}.ecsv")
    ids = list(passes(pre)["event_id"])
    fit_params = json.dumps(asdict(w3.P))
    parts = []
    for k in range(n):
        path = results_dir() / f"{chunk_name((k, n))}.gz"
        if not path.exists():
            raise SystemExit(f"missing chunk {k + 1}/{n}: {path}")
        tab = Table.read(path, format="ascii.ecsv")
        if tab.meta.get("fit_params") != fit_params:
            raise SystemExit(f"chunk {k + 1}/{n} was fitted with other Params; refit it")
        if sorted(tab["event_id"]) != sorted(ids[k::n]):
            raise SystemExit(f"chunk {k + 1}/{n} does not hold exactly its pre-screen passes")
        parts.append(tab)
    tab = vstack(parts, metadata_conflicts="silent")
    tab.sort("event_id")
    tab.meta["chunk"] = ""
    tab.meta["source"] = f"scripts/w3_moa.py merge-chunks: {n} chunk fit tables of {FIELD}"
    tab.meta["wall_time_s"] = round(sum(float(t.meta.get("wall_time_s", 0)) for t in parts), 1)
    path = out_dir() / f"fits_{FIELD}.ecsv"
    tab.write(path, overwrite=True)
    print(f"wrote {path}: {len(tab)} fits from {n} chunks")
    return path


# ----------------------------------------------------------------------------- vetting

REPEAT_S = 6.0  # ASSUMPTION: a second deficit this significant outside the feature = variable star
NEIGHBOUR_PX = 12.0  # ASSUMPTION: Cut-0 objects within 12 px (7″, ~3.5 seeing FWHM) share flux
BASELINE_CHI2 = 2.0  # ASSUMPTION (as D-057): χ²/dof of a constant outside the feature
NEIGHBOUR_S = 5.0  # ASSUMPTION: |S| of a neighbour's notch over the same window = shared feature


def feature_window(r: dict, model: str) -> tuple[float, float]:
    """Time span of an exotic fit's distinguishing feature: the umbra crossing of a repulsive lens
    (half-duration t_E √(β_c² − u0²)) or the bump of an attractive one, plus spike margins."""
    n, sign = w3.EXOTIC[model]
    te, u0 = r["tE"], r["u0"]
    if sign == -1:
        bc = es.caustic_beta(n)
        half = te * math.sqrt(max(bc * bc - u0 * u0, 0.25))
    else:
        half = te * max(u0, 0.5)
    pad = 0.3 * te + 2.0
    return r["t0"] - half - pad, r["t0"] + half + pad


def notch_in_window(t, f, sf, lo: float, hi: float) -> float:
    """S of the box [lo, hi] against flanks of the same length (``deficit_scan``'s statistic)."""
    t, f, sf = (np.asarray(x, float) for x in (t, f, sf))
    if t.size < 2 * P.n_min:
        return 0.0
    _b, s = robust_scale(f, sf)
    w = 1.0 / (sf * s) ** 2
    width = max(hi - lo, 1.0)
    box = (t >= lo) & (t <= hi)
    fl = ((t >= lo - width) & (t < lo)) | ((t > hi) & (t <= hi + width))
    if box.sum() < P.n_min or fl.sum() < P.n_min:
        return 0.0
    m_in = np.sum(w[box] * f[box]) / np.sum(w[box])
    m_fl = np.sum(w[fl] * f[fl]) / np.sum(w[fl])
    return float((m_in - m_fl) / math.sqrt(1.0 / np.sum(w[box]) + 1.0 / np.sum(w[fl])))


def with_regressors(lc: w3.LightCurve, aux: dict) -> w3.LightCurve:
    """Copy whose baseline has a constant plus linear terms in seeing, airmass and sky (the
    observing-condition systematics of difference photometry)."""
    rows = [np.ones_like(lc.t)]
    for k in AUX:
        v = np.asarray(aux.get(k, []), float)
        if v.size == lc.t.size and np.all(np.isfinite(v)) and np.std(v) > 0:
            rows.append((v - v.mean()) / v.std())
    out = w3.LightCurve(lc.t, lc.f, lc.sf, lc.ra, lc.dec, lc.event_id, lc.f_min)
    out.seasons = np.vstack(rows)
    out.n_extra = out.seasons.shape[0] - 1
    return out


def _start(model: str, r: dict) -> tuple:
    x = (r["t0"], math.log10(r["tE"]), r["u0"])
    if model in ("FSPL", *w3.EXOTIC):
        x = x + (math.log10(r["rho"]),)
    if model == "PAR":
        x = x + (r["pi_E_N"], r["pi_E_E"])
    return x


def refit_both(lcx: w3.LightCurve, res: dict, ex: str) -> float:
    """ΔBIC (exotic − best ordinary) on a modified light curve, each family restarted at its fit."""
    best_o = np.inf
    for m in w3.ORDINARY:
        if m in res:
            r = w3.optimise(m, lcx, [_start(m, res[m])], t0_par=res[m].get("t0_par"))
            best_o = min(best_o, r["bic"])
    e = w3.optimise(ex, lcx, [_start(ex, res[ex])])
    return float(e["bic"] - best_o)


def res_from_row(row) -> dict:
    """Fit results per model from one row of the ``fit`` table (``w3_microlensing.summarise``)."""
    cols = row.colnames
    res = {}
    n = int(row["n_points"])
    for m in (*w3.ORDINARY, *w3.EXOTIC):
        key = f"{m}_bic"
        if key not in cols or not np.isfinite(float(row[key])):
            continue
        r = {"model": m}
        for k in ("chi2", "dof", "bic", "t0", "tE", "u0", "rho", "pi_E_N", "pi_E_E", "fs", "fb"):
            c = f"{m}_{k}"
            if c in cols and np.isfinite(float(row[c])):
                r[k] = float(row[c])
        r["k"] = n - int(r["dof"])
        if m == "PAR":
            r["t0_par"] = round(r["t0"], 1)
        res[m] = r
    return res


SCREEN_MODELS = ("PSPL", "FSPL", *w3.EXOTIC)  # PAR only in vetting (it can only remove flags)


def vet_one(job) -> dict:
    """Ordinary explanations for one flag, cheapest first; the first failed test ends the vetting
    (``tests`` lists what ran; ``survives`` only if all ran and passed; ``complete`` then says the
    binary-lens fit ran). ``res``: the screen fits (``None``: fit here). ``neigh``: [(id, t, f,
    sf)] of nearby Cut-0 objects. The VSX/Gaia match is added by ``run_vet``."""
    ev, t, f, sf, ax, scan, neigh, binary_lens, res = job
    lc = to_lightcurve(t, f, sf, ev["ra"], ev["dec"], ev["event_id"], scan["err_scale"])
    out = {"event_id": ev["event_id"], "tests": [], "complete": False, "survives": False}

    def add(name, ok, note):
        out["tests"].append((name, bool(ok), note))
        return bool(ok)

    if res is None:
        res = fit_moa(lc, scan, models=SCREEN_MODELS)

    def summary():
        ex = min(w3.EXOTIC, key=lambda m: res[m]["bic"])
        ordinary = {m: res[m] for m in w3.ORDINARY if m in res}
        best_o = min(ordinary, key=lambda m: ordinary[m]["bic"])
        return ex, ordinary, best_o, res[ex]["bic"] - ordinary[best_o]["bic"]

    ex, ordinary, best_o, d0 = summary()
    out.update(exotic=ex, dbic_all=float(d0), best_ordinary=best_o)
    lo, hi = feature_window(res[ex], ex)
    out["feature_window"] = (lo, hi)

    def record():
        out["res"] = {m: {k: v for k, v in r.items() if k != "model"} for m, r in res.items()}
        return out

    if not add("screen_flag", d0 < P.flag_dbic, f"ΔBIC {d0:.1f} ({ex} vs {best_o})"):
        return record()
    # --- no fitting: shape and neighbourhood
    outside = (lc.t < lo) | (lc.t > hi)
    s2 = deficit_scan(lc.t[outside], lc.f[outside], lc.sf[outside])["z_min"]
    out["z_min_outside"] = s2
    if not add("repeated_deficit", s2 > -REPEAT_S, f"z_min outside the feature {s2:.1f}"):
        return record()
    far = (lc.t < lo - 60.0) | (lc.t > hi + 60.0)
    if far.sum() > 20:
        fo, wo = lc.f[far], lc.w[far]
        chi_out = float(np.sum((fo - np.sum(fo * wo) / np.sum(wo)) ** 2 * wo) / (far.sum() - 1))
    else:
        chi_out = np.nan
    out["chi2_out"] = chi_out
    if not add(
        "variable_baseline" if np.isfinite(chi_out) else "variable_baseline_untestable",
        not (chi_out > BASELINE_CHI2),
        f"baseline χ²/dof {chi_out:.2f} ({int(far.sum())} pts, errors × point-to-point scale)",
    ):
        return record()
    s_n = [(nid, notch_in_window(nt, nf, nsf, lo, hi)) for nid, nt, nf, nsf in neigh]
    hit = [(nid, round(v, 1)) for nid, v in s_n if abs(v) > NEIGHBOUR_S]
    out["neighbours"] = s_n
    if not add("neighbour_shares_feature", not hit, f"{len(s_n)} neighbours; |S| > 5: {hit}"):
        return record()
    # --- refits from the screen optimum: systematics models
    bad = w3.isolated_outliers(lc, w3.model_flux(best_o, lc, ordinary[best_o])) | (
        w3.isolated_outliers(lc, w3.model_flux(ex, lc, res[ex]))
    )
    lc2 = lc.subset(~bad)
    scale = max(ordinary[best_o]["chi2"] / ordinary[best_o]["dof"], 1.0)
    k_pen = (res[ex]["k"] - ordinary[best_o]["k"]) * math.log(lc2.t.size)
    d1 = (refit_both(lc2, res, ex) - k_pen) / scale + k_pen
    if not add(
        "robust_errors_outliers",
        d1 < P.flag_dbic,
        f"{int(bad.sum())} isolated outliers removed, errors ×{math.sqrt(scale):.2f}; "
        f"ΔBIC {d1:.1f}",
    ):
        return record()
    d2 = refit_both(lc.with_season_offsets(), res, ex)
    if not add("season_offsets", d2 < P.flag_dbic, f"free baseline per season; ΔBIC {d2:.1f}"):
        return record()
    d2t = refit_both(lc.with_season_offsets(trend=True), res, ex)
    if not add("season_trends", d2t < P.flag_dbic, f"offset + drift per season; ΔBIC {d2t:.1f}"):
        return record()
    d2c = refit_both(with_regressors(lc, ax), res, ex)
    if not add(
        "observing_conditions", d2c < P.flag_dbic, f"baseline linear in {AUX}; ΔBIC {d2c:.1f}"
    ):
        return record()
    # --- every ordinary single-lens model, then binary source and lens
    ps = res["PSPL"]
    lt = math.log10(ps["tE"])
    fs_ = w3.optimise(
        "FSPL", lc, [(ps["t0"], lt, ps["u0"], math.log10(r)) for r in (0.003, 0.03, 0.3)]
    )
    if "FSPL" not in res or fs_["bic"] < res["FSPL"]["bic"]:
        res["FSPL"] = fs_
    if w3.have_mm():
        res["PAR"] = w3.optimise("PAR", lc, w3.par_starts(ps), t0_par=round(ps["t0"], 1))
    ex, ordinary, best_o, d0 = summary()
    out.update(exotic=ex, dbic_all=float(d0), best_ordinary=best_o)
    if not add("refit_all_ordinary", d0 < P.flag_dbic, f"ΔBIC {d0:.1f} vs {best_o} (FSPL, PAR)"):
        return record()
    bs = w3.fit_binary_source(lc, ps)
    d3 = res[ex]["bic"] - min(ordinary[best_o]["bic"], bs["bic"])
    if not add("binary_source", d3 < P.flag_dbic, f"ΔBIC {d3:.1f}"):
        return record()
    if binary_lens and w3.have_mm():
        bl = w3.fit_binary_lens(lc, ps)
        d4 = res[ex]["bic"] - min(ordinary[best_o]["bic"], bs["bic"], bl["bic"])
        if not add("binary_lens", d4 < P.flag_dbic, f"ΔBIC {d4:.1f}"):
            return record()
    # --- D-057 revet tests: feature sampled, jackknife, two unrelated events
    cov = w3.feature_coverage(lc, ordinary[best_o], res[ex], best_o, ex)
    out["feature"] = cov
    ok_cov = cov["n_epochs"] >= 3 and cov["dchi2_in"] < P.flag_dbic
    if not add(
        "exotic_feature_sampled",
        ok_cov,
        f"{cov['n_epochs']} epochs inside, Δχ² in {cov['dchi2_in']:.1f} / out "
        f"{cov['dchi2_out']:.1f}",
    ):
        return record()
    n_drop = w3.jackknife_n_drop(cov["n_epochs"])
    if n_drop:
        jk = w3.jackknife_worst_epochs(lc, ordinary[best_o], res[ex], best_o, ex, n_drop)
        if not add("jackknife_epochs", jk[-1] < P.flag_dbic, f"ΔBIC {[round(v, 1) for v in jk]}"):
            return record()
    d2e, _two = w3.two_events_dbic(lc, ordinary["PSPL"], res[ex], "PSPL", ex)
    if d2e is not None and not add("two_unrelated_events", d2e < P.flag_dbic, f"ΔBIC {d2e:.1f}"):
        return record()
    out["complete"] = bool(binary_lens and w3.have_mm())
    out["survives"] = True
    return record()


def _vet_worker(job):
    try:
        return vet_one(job)
    except Exception as exc:  # noqa: BLE001 — a failed vetting keeps the flag open
        return {
            "event_id": job[0]["event_id"],
            "tests": [("vet_error", True, repr(exc)[:200])],
            "survives": True,
            "complete": False,
        }


def neighbours_of(ev: Table, eid: str, radius_px: float = NEIGHBOUR_PX) -> list[str]:
    """Cut-0 objects of the same chip and subframe within ``radius_px`` pixels."""
    r = ev[ev["event_id"] == eid][0]
    same = (ev["chip"] == r["chip"]) & (ev["subframe"] == r["subframe"])
    d = np.hypot(ev["x"] - r["x"], ev["y"] - r["y"])
    return [str(x) for x in ev["event_id"][same & (d < radius_px) & (ev["event_id"] != eid)]]


def run_vet(procs: int) -> Path:
    field = moa.MoaField(FIELD)
    fits = Table.read(out_dir() / f"fits_{FIELD}.ecsv")
    ok = no_error(fits)
    flags = fits[ok & (np.asarray(fits["dbic_min"], float) < P.flag_dbic)]
    ev = field.events()
    pos = {r["event_id"]: r for r in ev}
    neigh = {str(e): neighbours_of(ev, str(e)) for e in flags["event_id"]}
    ids = set(map(str, flags["event_id"])) | {n for v in neigh.values() for n in v}
    arrays = load_arrays(field, ids, aux=True)
    jobs = []
    for r in flags:
        eid = str(r["event_id"])
        t, f, sf, ax = arrays[eid]
        scan = {k: float(r[k]) for k in SCAN_KEYS}
        nb = [(n, *arrays[n][:3]) for n in neigh[eid] if n in arrays]
        d = {"event_id": eid, "ra": float(pos[eid]["ra"]), "dec": float(pos[eid]["dec"])}
        jobs.append((d, t, f, sf, ax, scan, nb, True, res_from_row(r)))
    print(f"vetting {len(jobs)} flags", flush=True)
    t1 = time.time()
    out = []
    with Pool(procs) as pool:
        for i, o in enumerate(pool.imap_unordered(_vet_worker, jobs, 1)):
            out.append(o)
            if (i + 1) % 10 == 0:
                print(f"{i + 1}/{len(jobs)} vetted, {time.time() - t1:.0f} s", flush=True)
    alive = [o for o in out if o.get("survives")]
    var = (
        w3.variable_catalogue_matches(
            [float(pos[o["event_id"]]["ra"]) for o in alive],
            [float(pos[o["event_id"]]["dec"]) for o in alive],
            radius_arcsec=2.0,
        )
        if alive
        else {}
    )
    failed = [c for c, idx in var.items() if any(isinstance(j, str) for j in idx)]
    for i, o in enumerate(alive):
        hits = [c for c, idx in var.items() if i in idx]
        o["tests"].append(("variable_catalogues", not hits, f"matches: {hits or 'none'} (2″)"))
        if failed:
            o["tests"].append(("variable_catalogues_failed", True, f"queries failed: {failed}"))
            o["complete"] = False
        o["survives"] = all(ok for _, ok, _ in o["tests"])
    path = out_dir() / f"vetting_{FIELD}.json"
    rec = {
        "provenance": "derived",
        "n_fit": len(fits),
        "n_flags": len(flags),
        "wall_time_s": time.time() - t1,
        "variable_xmatch": var,
        "flags": sorted(out, key=lambda o: o["event_id"]),
    }
    path.write_text(json.dumps(rec, indent=1, default=float))
    print(f"wrote {path}; survivors: {[o['event_id'] for o in out if o.get('survives')]}")
    print(vetting_funnel(rec))
    return path


def vetting_funnel(vet: dict) -> list[tuple[str, int]]:
    """Flags still alive after each test, in the order the tests ran."""
    order = []
    for o in vet["flags"]:
        for name, _ok, _note in o["tests"]:
            if name not in order:
                order.append(name)
    alive = {o["event_id"] for o in vet["flags"]}
    out = [("flags", len(alive))]
    for name in order:
        for o in vet["flags"]:
            for nm, ok, _ in o["tests"]:
                if nm == name and not ok:
                    alive.discard(o["event_id"])
        out.append((name, len(alive)))
    return out


def contact_sheet(path_png: Path, ids=None, max_panels: int = 24) -> None:
    """Light curves of vetted flags around the feature with the best ordinary and exotic models."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    vet = json.loads((out_dir() / f"vetting_{FIELD}.json").read_text())["flags"]
    vet = [o for o in vet if "res" in o and (ids is None or o["event_id"] in ids)]
    vet = sorted(vet, key=lambda o: (not o.get("survives"), o["dbic_all"]))[:max_panels]
    if not vet:
        return
    arrays = load_arrays(moa.MoaField(FIELD), [o["event_id"] for o in vet])
    ncol = 3
    nrow = math.ceil(len(vet) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 3.3 * nrow), squeeze=False)
    for ax, o in zip(axes.ravel(), vet, strict=False):
        t, f, sf = arrays[o["event_id"]]
        lo, hi = o["feature_window"]
        span = max(hi - lo, 20.0)
        sel = (t > lo - span) & (t < hi + span)
        ax.errorbar(t[sel] - lo, f[sel], sf[sel], fmt=".", ms=2, color="0.3", lw=0.4)
        tt = np.linspace(lo - span, hi + span, 3000)
        fine = w3.LightCurve(tt, np.ones_like(tt), np.ones_like(tt), None, None, f_min=np.inf)
        for m, c in ((o["best_ordinary"], "k"), (o["exotic"], "tab:red")):
            mm, r = (m, o["res"][m]) if m != "PAR" else ("PSPL", o["res"]["PSPL"])
            ax.plot(tt - lo, w3.model_flux(mm, fine, r), color=c, lw=0.8, label=m)
        failed = [n for n, ok, _ in o["tests"] if not ok]
        ax.set_title(
            f"{o['event_id']} {o['exotic']} ΔBIC {o['dbic_all']:.0f}\n"
            f"fails: {', '.join(failed[:3]) or 'none'}",
            fontsize=7,
        )
        ax.legend(fontsize=6)
        ax.set_xlabel(f"HJD − {lo:.1f}", fontsize=7)
        ax.set_ylabel("difference flux (counts)", fontsize=7)
        ax.tick_params(labelsize=6)
    for ax in axes.ravel()[len(vet) :]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(path_png, dpi=100)


# ----------------------------------------------------------------------------- injection

INJ_TE = (3.0, 10.0, 30.0, 100.0, 300.0)  # days (ASSUMPTION: grid, as D-057)
INJ_RHO = (0.01, 0.1)
INJ_IS = (14.2, 21.4)  # source MOA-Red magnitudes, uniform (Koshimoto et al. 2023's range)
LF_SLOPE = 0.319  # dex/mag: Gaia DR3 RP counts in gb22, RP 13–17.5 (ASSUMPTION: extends to 21.4)


def quiet_carriers(pre: Table, n: int, seed: int) -> list[str]:
    """Light curves with no significant notch either way (``|z| < 4``) and ≥ 1,000 epochs: stand-ins
    for the difference light curve of a constant star (ASSUMPTION; constant stars are not Cut-0
    objects, so their own light curves are not in the release)."""
    ok = (
        no_error(pre)
        & (np.asarray(pre["z_min"], float) > -4.0)
        & (np.asarray(pre["z_max"], float) < 4.0)
        & (np.asarray(pre["s_min"], float) > -P.prescreen_s)
        & (np.asarray(pre["n_points"]) >= 1000)
    )
    ids = np.asarray(pre["event_id"][ok], str)
    rng = np.random.default_rng(seed)
    return list(rng.choice(ids, size=min(n, ids.size), replace=False))


def injected_signal(t, kind: str, fs: float, prm: dict) -> np.ndarray:
    """Difference flux F_s (A − 1) of an injected event (``simulated``; ``exotic_sim``)."""
    n, sign, rho = (1.0, -1, prm["rho"]) if kind == "W3" else (1.0, 1, 0.0)
    inj = es.inject_light_curve(t, fs, prm["t0"], prm["tE"], prm["u0"], n, sign, rho, 1.0)
    return np.asarray(inj["flux"], float) - fs


def _inject_worker(job):
    kind, eid, t, f, sf, ax, prm = job
    row = {"kind": kind, "event_id": eid, **prm}
    try:
        chip = moa.parse_event_id(eid)[1]
        fs = float(moa.mag_to_counts(prm["Is"], chip))
        sig = injected_signal(t, kind, fs, prm)
        f2 = f + sig
        row["cut0"] = cut0_emulated(t, sig, sf)
        scan = deficit_scan(t, f2, sf)
        row.update(z_min=scan["z_min"], s_min=scan["s_min"], z_min2=scan["z_min2"])
        row["prescreen"] = bool(prescreen_pass(scan["z_min"], scan["s_min"], scan["z_min2"]))
        row.update(flagged=False, survives=False, failed_test="")
        if row["cut0"] and row["prescreen"] and not prm.get("prescreen_only"):
            t1 = time.time()
            d = {"event_id": eid, "ra": prm["ra"], "dec": prm["dec"]}
            o = vet_one((d, t, f2, sf, ax, {k: scan[k] for k in SCAN_KEYS}, [], False, None))
            row["flagged"] = bool(o["tests"][0][1])
            row["survives"] = bool(o["survives"])
            row["failed_test"] = next((nm for nm, ok, _ in o["tests"] if not ok), "")
            row["dbic_all"] = o["dbic_all"]
            row["exotic"] = o["exotic"]
            row["seconds"] = time.time() - t1
        row["recovered"] = bool(row["cut0"] and row["prescreen"] and row["survives"])
        row["error"] = ""
    except Exception as exc:  # noqa: BLE001
        row["error"] = repr(exc)[:200]
    return row


def run_inject(
    procs: int, per_cell: int, per_ctrl: int, seed: int = 60, prescreen_only: bool = False
) -> Path:
    """Injection-recovery through the whole chain: Cut-0 emulation, pre-screen, fit, vetting
    (binary lens and the VSX/Gaia match excepted). W3: n = 1, ε < 0, u0 ~ U[0, 2); PSPL controls:
    u0 ~ U[0, 1) (they measure how often an ordinary event ends as a W3 survivor)."""
    field = moa.MoaField(FIELD)
    pre = Table.read(out_dir() / f"prescreen_{FIELD}.ecsv")
    carriers = quiet_carriers(pre, 300, seed)
    arrays = load_arrays(field, carriers, aux=True)
    ev = field.events()
    pos = {r["event_id"]: (float(r["ra"]), float(r["dec"])) for r in ev}
    rng = np.random.default_rng(seed)
    jobs = []
    for te in INJ_TE:
        for kind, rhos, n_per, umax in (
            ("W3", INJ_RHO, per_cell, 2.0),
            ("PSPL", (0.0,), per_ctrl, 1.0),
        ):
            for rho in rhos:
                for _ in range(n_per):
                    eid = carriers[int(rng.integers(len(carriers)))]
                    t, f, sf, ax = arrays[eid]
                    prm = {
                        "tE": te,
                        "rho": rho,
                        "u0": float(rng.uniform(0.0, umax)),
                        "t0": float(rng.uniform(moa.T_START, moa.T_END)),
                        "Is": float(rng.uniform(*INJ_IS)),
                        "ra": pos[eid][0],
                        "dec": pos[eid][1],
                    }
                    if prescreen_only:
                        prm["prescreen_only"] = True
                    jobs.append((kind, eid, t, f, sf, ax, prm))
    print(f"{len(jobs)} injections on {len(arrays)} carriers", flush=True)
    t1 = time.time()
    rows = []
    with Pool(procs) as pool:
        for i, row in enumerate(pool.imap_unordered(_inject_worker, jobs, 1)):
            rows.append(row)
            if (i + 1) % 50 == 0:
                print(f"{i + 1}/{len(jobs)} injections, {time.time() - t1:.0f} s", flush=True)
    tab = rows_to_table(rows, {"error": "", "failed_test": "", "exotic": "", "kind": ""})
    tab.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source=(
            "scripts/w3_moa.py inject: exotic_sim.inject_light_curve (simulated) added to real "
            f"{field.name} difference light curves of quiet Cut-0 objects"
        ),
        params=json.dumps(asdict(P)),
        wall_time_s=round(time.time() - t1, 1),
        seed=seed,
    )
    path = out_dir() / f"injections_{FIELD}{'_prescreen' if prescreen_only else ''}.ecsv"
    tab.write(path, overwrite=True)
    print(f"wrote {path}: {len(tab)} rows, {time.time() - t1:.0f} s")
    return path


# ----------------------------------------------------------------------------- limit


def lf_weights(mag) -> np.ndarray:
    """Luminosity-function weights ∝ 10^(LF_SLOPE · I) (ASSUMPTION, see ``LF_SLOPE``)."""
    return 10.0 ** (LF_SLOPE * (np.asarray(mag, float) - INJ_IS[1]))


def lf_fraction_injected() -> float:
    """Fraction of the 10 ≤ I ≤ 21.4 stars (N_s) inside the injected range 14.2–21.4."""
    a = LF_SLOPE * math.log(10.0)

    def integral(lo, hi):
        return (math.exp(a * hi) - math.exp(a * lo)) / a

    return integral(*INJ_IS) / integral(10.0, INJ_IS[1])


def run_limit() -> Path:
    """95 % limit per monitored star per year; needs a complete null vetting (zero survivors)."""
    vet = json.loads((out_dir() / f"vetting_{FIELD}.json").read_text())
    open_flags = [o["event_id"] for o in vet["flags"] if o.get("survives")]
    if open_flags:
        raise SystemExit(f"no zero-event limit: flags survive: {open_flags}")
    inj = Table.read(out_dir() / f"injections_{FIELD}.ecsv")
    inj = inj[no_error(inj)]
    n_s, n_lo, n_hi = moa.star_count_estimate(moa.CUT0_PER_FIELD[22])
    years = (moa.T_END - moa.T_START) / 365.25
    frac = lf_fraction_injected()
    kind = np.asarray(inj["kind"], str)
    rows = []
    for te in INJ_TE:
        ctrl = inj[(kind == "PSPL") & (np.asarray(inj["tE"]) == te)]
        for rho in INJ_RHO:
            w = inj[
                (kind == "W3") & (np.asarray(inj["tE"]) == te) & (np.asarray(inj["rho"]) == rho)
            ]
            rec = np.asarray(w["recovered"], bool)
            wt = lf_weights(w["Is"])
            eff_inj = float(np.sum(wt * rec) / np.sum(wt)) if len(w) else np.nan
            eff = eff_inj * frac  # stars brighter than 14.2 counted with efficiency 0
            lim = 3.0 / (n_s * years * eff) if eff > 0 else np.inf
            rows.append(
                {
                    "tE_days": te,
                    "rho": rho,
                    "mass_msun_model": (te / w3.einstein_time_days(1.0)[0]) ** 2,
                    "n_inj": len(w),
                    "p_cut0": float(np.mean(w["cut0"])),
                    "p_prescreen": float(
                        np.mean(np.asarray(w["cut0"]) & np.asarray(w["prescreen"]))
                    ),
                    "p_flag": float(np.mean(np.asarray(w["flagged"], bool))),
                    "p_recovered": float(rec.mean()),
                    "p_recovered_bright": float(rec[np.asarray(w["Is"]) < 19.0].mean())
                    if (np.asarray(w["Is"]) < 19.0).any()
                    else np.nan,
                    "eff_per_star": eff,
                    "rate95_per_star_yr": lim,
                    "rate95_conservative": 3.0 / (n_lo * years * eff) if eff > 0 else np.inf,
                    "n_ctrl": len(ctrl),
                    "ctrl_false_w3": int(np.sum(np.asarray(ctrl["recovered"], bool))),
                }
            )
    tab = Table(rows)
    tab.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source="scripts/w3_moa.py limit: injection-recovery through the full gb22 chain",
        assumptions=(
            f"N_s = {n_s:.3g} (range {n_lo:.3g}–{n_hi:.3g}) monitored stars, 10 ≤ I ≤ 21.4, from "
            "N_s per Cut-0 object of Nunota et al. 2024's 20 fields; T = "
            f"{years:.2f} yr; LF ∝ 10^({LF_SLOPE} I) (Gaia DR3 RP); injections I_s ~ "
            "U[14.2, 21.4], "
            "MOA-Red ≈ I; carrier noise unchanged by the source (sky/blend dominated); Cut-0 "
            "emulated at light-curve level; 95 % Poisson for 0 events = 3.0; mass: n = 1, "
            "D_L = 4 kpc, D_S = 8 kpc, μ_rel = 5 mas/yr (model_prediction)"
        ),
    )
    path = out_dir() / f"limits_{FIELD}.ecsv"
    tab.write(path, overwrite=True)
    tab.pprint(max_width=250, max_lines=50)
    return path


LAST_MODIFIED = {  # HTTP Last-Modified of the pinned files, read 2026-10-08 (the release version)
    "metadata.ipac.tar.gz": "2023-10-13",
    "gb22.tar": "2023-10-10",
}


def write_manifest() -> Path:
    """Pinned MOA-II files as a tracked manifest (URL, sha256, size, retrieval date, release)."""
    names = [u.rsplit("/", 1)[-1] for u in moa.FILES]
    tab = Table(
        {
            "url": list(moa.FILES),
            "sha256": [v[0] for v in moa.FILES.values()],
            "size_bytes": [v[1] for v in moa.FILES.values()],
            "retrieved_utc": ["2026-10-08"] * len(names),
            "last_modified": [LAST_MODIFIED[n] for n in names],
            "release": ["MOA-II 9-year bulge release (2006-2014), NASA Exoplanet Archive"]
            * len(names),
        }
    )
    tab.meta.update(
        provenance=schema.Provenance.OBSERVED.value,
        source="jwst_anomaly.moa.FILES (MOA-II 9-year release; D-062)",
    )
    path = paths.manifests_dir() / "moa_ii.ecsv"
    tab.write(path, overwrite=True)
    return path


# ----------------------------------------------------------------------------- main


# The fitter's bounds for MOA (ASSUMPTION): t_E ≤ 1,000 d (longer "events" in a 3,146-d light curve
# are baseline trends, and an exact finite-source integral on every epoch made one fit take 200 s)
# and ρ ≤ 0.3 (the injections use ρ ≤ 0.1). Applied by the CLI only, so tests and
# w3_microlensing keep their defaults.
MOA_FIT_BOUNDS = {"te_bounds": (0.1, 1000.0), "log_rho_bounds": (-3.5, -0.5)}


def use_moa_fit_params() -> None:
    w3.P = replace(w3.P, **MOA_FIT_BOUNDS)


def main(argv=None) -> int:
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    use_moa_fit_params()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prescreen", help="deficit statistic on every light curve of the field")
    p.add_argument("--procs", type=int, default=4)
    f = sub.add_parser("fit", help="ordinary and exotic fits of the pre-screen passes")
    f.add_argument("--procs", type=int, default=4)
    f.add_argument("--limit", type=int, default=None)
    f.add_argument("--chunk", default=None, help="K/N: fit every N-th pass from the K-th")
    m = sub.add_parser("merge-chunks", help="join the tracked chunk fit tables 1..N")
    m.add_argument("--n", type=int, required=True)
    v = sub.add_parser("vet", help="vet the flags of the fit stage")
    v.add_argument("--procs", type=int, default=4)
    c = sub.add_parser("sheet", help="contact sheet of the vetted flags")
    c.add_argument("--out", type=Path, default=None)
    c.add_argument("--survivors-only", action="store_true")
    i = sub.add_parser("inject", help="injection-recovery through the whole chain")
    i.add_argument("--procs", type=int, default=4)
    i.add_argument("--per-cell", type=int, default=30)
    i.add_argument("--per-ctrl", type=int, default=20)
    i.add_argument("--prescreen-only", action="store_true", help="Cut-0 and pre-screen only")
    sub.add_parser("limit", help="95 %% rate limit per monitored star per year")
    sub.add_parser("manifest", help="write data/manifests/moa_ii.ecsv")
    a = ap.parse_args(argv)
    if a.cmd == "prescreen":
        run_prescreen(min(a.procs, 4))
    elif a.cmd == "fit":
        run_fit(min(a.procs, 4), a.limit, w3.parse_chunk(a.chunk))
    elif a.cmd == "merge-chunks":
        merge_chunks(a.n)
    elif a.cmd == "vet":
        run_vet(min(a.procs, 4))
    elif a.cmd == "inject":
        run_inject(min(a.procs, 4), a.per_cell, a.per_ctrl, prescreen_only=a.prescreen_only)
    elif a.cmd == "limit":
        run_limit()
    elif a.cmd == "manifest":
        print(write_manifest())
    elif a.cmd == "sheet":
        vet = json.loads((out_dir() / f"vetting_{FIELD}.json").read_text())
        ids = (
            [o["event_id"] for o in vet["flags"] if o.get("survives")] if a.survivors_only else None
        )
        png = a.out or out_dir() / f"contact_sheet_{FIELD}.png"
        contact_sheet(png, ids)
        print(png, vetting_funnel(vet))
    return 0


if __name__ == "__main__":
    sys.exit(main())
