"""W3 in the MOA-II 9-year release, any field gb1 … gb22: pre-screen, fits, vetting, limit.

D-062 (method, gb22 pilot) and D-068 (streaming, calibration, all fields). The release holds every
Cut-0 variable object (difference-image detections of positive *or negative* PSF profiles;
``jwst_anomaly.moa``), before any bump or PSPL cut, so a W3 event (the source flux drops toward
zero inside an umbra between two caustic spikes; ``exotic_sim``) can be in it. Fitting every light
curve with every model is too slow, so (``--field gbF``, default gb22):

- ``prescreen``: streams the field tar with concurrent HTTP range reads (``moa_stream``; the tar is
  never stored) into a process pool computing a W3-shaped matched statistic on every light curve
  (``deficit_scan``): the most significant box of width W = 1 … 300 d that is below both its
  flanks and the median flux, in units of the light curve's own spread of that statistic (red
  noise), and no second such deficit elsewhere (``prescreen_pass``). One tracked table per 4 GiB
  of tar in ``results/w3_moa/prescreen/`` (resumable); ``merge-prescreen`` joins them.
- ``fit``: the passes are fitted with the ordinary (PSPL, FSPL, PAR) and exotic (N1neg, E2pos,
  E2neg) models of ``scripts/w3_microlensing.py`` on one trajectory, blend flux free (difference
  flux), extra exotic starts on the deficit. Flag: ΔBIC < −10 (ASSUMPTION, as D-057).
- ``vet``: flags through ordinary explanations, cheapest first (``vet_one``); the variable-baseline
  threshold is calibrated on the field's quiet light curves (``calibrate_baseline``).
- ``inject``: W3 events (``exotic_sim``, ``simulated``) added to real light curves of the field's
  quiet objects, through a light-curve-level Cut-0 emulation, the pre-screen, the fit and vetting.
- ``limit``: 95 % upper limit on the W3 rate per monitored star per year for the field;
  ``combine``: the same over every field with a tracked limit table.
- ``run-field``: every stage in order, resumable.

Outputs go to ``$JWST_ANOMALY_DATA/derived/w3_moa/``, except the small tracked tables in
``results/w3_moa/``. Exotic physics is a hypothesis: a flag is an anomaly to vet, never a discovery.
"""

from __future__ import annotations

import argparse
import functools
import json
import math
import os
import sys
import time
import zlib
from dataclasses import asdict, dataclass, replace
from multiprocessing import Pool
from pathlib import Path

# one thread per process: the stages run process pools sized to the cores (never oversubscribe);
# set before numpy loads its BLAS
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402
from astropy.table import Table, vstack  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))

import w3_microlensing as w3  # noqa: E402

from jwst_anomaly import exotic_sim as es  # noqa: E402
from jwst_anomaly import moa, moa_stream, paths, schema  # noqa: E402

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


def _median(x: np.ndarray) -> float:
    """``np.median`` of a 1-d array without its wrapper overhead (bit-identical: the mean of the
    two middle values is their sum / 2, as ``np.mean`` computes it; NaN or empty: ``np.median``)."""
    n = x.size
    if n == 0 or np.isnan(x).any():
        return float(np.median(x))
    lo, hi = np.partition(x, ((n - 1) // 2, n // 2))[[(n - 1) // 2, n // 2]]
    return float(lo) if n % 2 else float((lo + hi) / 2.0)


def robust_scale(f: np.ndarray, sf: np.ndarray) -> tuple[float, float]:
    """Median flux and the error scale from point-to-point scatter (1.4826 MAD of successive
    differences of f/σ, / √2, ≥ 1): insensitive to slow trends and to the event itself."""
    b = _median(f)
    if f.size < 3:
        return b, 1.0
    d = np.diff(f) / np.hypot(sf[1:], sf[:-1])
    s = 1.4826 * _median(np.abs(d - _median(d)))
    return b, max(s, 1.0)


def nights(t: np.ndarray) -> np.ndarray:
    """Night number (MOA, New Zealand: nights split at 00:00 UT = JD fraction 0.5)."""
    return np.floor(np.asarray(t) - 0.5).astype(np.int64)


def spread_of(x: np.ndarray, valid: np.ndarray) -> float:
    """Robust spread (1.4826 MAD, ≥ 1) of a box statistic over the valid boxes."""
    if valid.sum() < 10:
        return 1.0
    v = x[valid]
    return max(1.4826 * _median(np.abs(v - _median(v))), 1.0)


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
                s_min=float(max(st[k], sg[k])),  # the weaker of the two, as z
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


PRESCREEN_COLUMNS = ("HJD", "flux", "cor_flux", "flux_err", "included")


def baseline_chi2(f, sf, err_scale: float) -> float:
    """χ²/dof of a constant (weighted mean) with errors × the point-to-point scale: the statistic
    of the variable-baseline vetting test, here on a whole light curve (``derived``)."""
    f = np.asarray(f, float)
    if f.size < 3:
        return float("nan")
    w = 1.0 / (np.asarray(sf, float) * err_scale) ** 2
    m = np.sum(w * f) / np.sum(w)
    return float(np.sum((f - m) ** 2 * w) / (f.size - 1))


def scan_member(eid: str, raw: bytes) -> dict:
    """Pre-screen row of one light-curve member (gzipped IPAC bytes)."""
    try:
        cols = moa.parse_lightcurve(raw, columns=PRESCREEN_COLUMNS)
        t, f, sf, kind = moa.select_flux(cols)
        ok = np.isfinite(f) & np.isfinite(sf) & (sf > 0)
        t, f, sf = t[ok], f[ok], sf[ok]
        row = {"event_id": eid, "flux_kind": kind, **deficit_scan(t, f, sf)}
        row["chi2_const"] = baseline_chi2(f, sf, row["err_scale"])
        row["error"] = ""
    except Exception as exc:  # noqa: BLE001 — one bad file must not stop the pass
        row = {"event_id": eid, "error": repr(exc)[:200]}
    return row


def _scan_batch(batch):
    """Worker: ``[(event_id, offset_data, size, gz_bytes)]`` -> pre-screen rows with offsets."""
    out = []
    for eid, off, size, raw in batch:
        row = scan_member(eid, raw)
        row["offset"], row["size"] = int(off), int(size)
        out.append(row)
    return out


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


CHUNK_BYTES = (
    4 * 2**30
)  # tar bytes per tracked pre-screen chunk (~1–2 min of download, ~1 min CPU/GB)
SEG_BYTES = 64 * 2**20  # one HTTP range read
QUIET_TRACK_MOD = 4  # ASSUMPTION: 1 in 4 quiet light curves (by crc32 of the id) is tracked
TRACK_COLUMNS = (
    "event_id",
    "flux_kind",
    "z_min",
    "s_min",
    "z_max",
    "z_min2",
    "spread",
    "width",
    "t_lo",
    "t_hi",
    "err_scale",
    "n_points",
    "chi2_const",
    "offset",
    "size",
    "error",
)


def field_number(field: str | None = None) -> int:
    return int((field or FIELD).removeprefix("gb"))


def n_chunks(field: str | None = None) -> int:
    return math.ceil(moa.TAR_BYTES[field_number(field)] / CHUNK_BYTES)


def local_tar(field: str | None = None) -> Path | None:
    """The whole tar under ``raw/moa/`` if a session downloaded it (only gb22 is pinned)."""
    url = moa.tar_url(field_number(field))
    if url not in moa.FILES:
        return None
    p = paths.data_root() / "raw" / "moa" / f"{moa.FILES[url][0][:12]}_{url.rsplit('/', 1)[-1]}"
    return p if p.exists() else None


def reader_for(field: str | None = None, source: str = "auto") -> moa_stream.RangeReader:
    """``source``: ``auto`` (local tar if present, else HTTP), ``local`` or ``http``."""
    p = None if source == "http" else local_tar(field)
    if source == "local" and p is None:
        raise SystemExit(f"no local tar for {field or FIELD}")
    if p is not None:
        return moa_stream.RangeReader(path=p)
    f = field_number(field)
    return moa_stream.RangeReader(url=moa.tar_url(f), total=moa.TAR_BYTES[f])


def is_quiet(tab: Table) -> np.ndarray:
    """No significant notch either way (|z| < 4), S > −5 and ≥ 1,000 epochs: stand-ins for the
    difference light curve of a constant star (ASSUMPTION; see ``quiet_carriers``)."""
    return (
        no_error(tab)
        & (np.asarray(tab["z_min"], float) > -4.0)
        & (np.asarray(tab["z_max"], float) < 4.0)
        & (np.asarray(tab["s_min"], float) > -P.prescreen_s)
        & (np.asarray(tab["n_points"], float) >= 1000)
    )


def tracked_mask(tab: Table) -> np.ndarray:
    """Rows kept in the tracked chunk tables: every error, every deficit in the shared-epoch
    population (z < −``COINC_Z``, a superset of the pre-screen passes) and a deterministic 1 in
    ``QUIET_TRACK_MOD`` sample of the quiet light curves (injection carriers and the baseline
    calibration)."""
    ids = np.asarray(tab["event_id"], str)
    sample = np.array([zlib.crc32(e.encode()) % QUIET_TRACK_MOD == 0 for e in ids], bool)
    z = np.asarray(tab["z_min"], float) if "z_min" in tab.colnames else np.zeros(len(tab))
    return ~no_error(tab) | (z < -COINC_Z) | (is_quiet(tab) & sample)


def prescreen_chunk_path(field: str, k: int, n: int) -> Path:
    d = results_dir() / "prescreen"
    d.mkdir(parents=True, exist_ok=True)
    return d / f"{field}_c{k + 1:03d}of{n:03d}.ecsv.gz"


def _chunk_done(field: str, k: int, n: int) -> bool:
    p = prescreen_chunk_path(field, k, n)
    if not p.exists():
        return False
    tab = Table.read(p, format="ascii.ecsv")
    return (
        tab.meta.get("params") == json.dumps(asdict(P))
        and tab.meta.get("chunk_bytes") == CHUNK_BYTES
        and set(TRACK_COLUMNS) <= set(tab.colnames)
    )


def _write_prescreen_chunk(field: str, k: int, n: int, rows: list, stats: dict) -> Path:
    tab = rows_to_table(rows, {"error": "", "flux_kind": ""})
    for c in TRACK_COLUMNS:
        if c not in tab.colnames:
            tab[c] = [""] * len(tab) if c in ("flux_kind", "error") else np.full(len(tab), np.nan)
    tab.sort("event_id")
    ok = no_error(tab)
    shape = ok & prescreen_pass(tab["z_min"], tab["s_min"], tab["z_min2"])
    meta = {
        "provenance": schema.Provenance.DERIVED.value,
        "source": (
            f"scripts/w3_moa.py stream-prescreen: deficit_scan on every light curve of "
            f"{moa.tar_url(field_number(field))} bytes [{k * CHUNK_BYTES}, "
            f"{min((k + 1) * CHUNK_BYTES, moa.TAR_BYTES[field_number(field)])})"
        ),
        "params": json.dumps(asdict(P)),
        "chunk_bytes": CHUNK_BYTES,
        "chunk": f"{k + 1}/{n}",
        "n_members": len(tab),
        "n_errors": int((~ok).sum()),
        "n_shape": int(shape.sum()),
        "n_quiet": int(is_quiet(tab).sum()),
        "quiet_track_mod": QUIET_TRACK_MOD,
        "quiet_chi2_hist": chi2_histogram(np.asarray(tab["chi2_const"], float)[is_quiet(tab)]),
        **stats,
    }
    full = out_dir() / "prescreen_full"
    full.mkdir(exist_ok=True)
    tab.meta.update(meta)
    tab.write(full / f"{field}_c{k + 1:03d}of{n:03d}.ecsv", overwrite=True)  # untracked
    small = tab[tracked_mask(tab)][list(TRACK_COLUMNS)]
    small.meta.update(meta)
    path = prescreen_chunk_path(field, k, n)
    w3.write_ecsv_gz(small, path)
    return path


def _stream_worker(wid: int, segs: list, total: int, url, path, pin, prefetch: int, queue) -> None:
    """One pre-screen process: reads its own byte ranges (``prefetch`` threads download the next
    segments while this one is scanned) and puts ``(segment, rows, MB, retries)`` on ``queue``;
    the parent holds no light-curve bytes. ``None`` marks the end, ``("error", text)`` a failure."""
    from collections import deque
    from concurrent.futures import ThreadPoolExecutor

    moa_stream.limit_heap_growth()
    try:
        reader = moa_stream.RangeReader(url=url, path=path, total=pin)
        with ThreadPoolExecutor(prefetch) as ex:
            futs = deque()
            nxt = 0
            while nxt < len(segs) and len(futs) < prefetch:
                k, a, b = segs[nxt]
                futs.append(ex.submit(moa_stream.read_segment_hashed, reader, a, b, total))
                nxt += 1
            for s in segs:
                members, sha = futs.popleft().result()
                if nxt < len(segs):
                    k, a, b = segs[nxt]
                    futs.append(ex.submit(moa_stream.read_segment_hashed, reader, a, b, total))
                    nxt += 1
                rows = _scan_batch(members)
                del members
                queue.put((s, rows, reader.n_retries, sha))
        queue.put(None)
    except Exception as exc:  # noqa: BLE001 — reported to the parent, which stops the run
        queue.put(("error", f"worker {wid}: {exc!r}"[:500]))


def run_stream_prescreen(field: str, procs: int, conns: int, chunks=None, source: str = "auto",
                         log_every: float = 30.0) -> list[Path]:  # fmt: skip
    """Pre-screen a whole field tar without storing it. ``procs`` processes each read their own
    64 MB byte ranges (segments dealt round-robin) with ``conns / procs`` prefetch threads, so
    downloads overlap the scan and no process holds more than a few segments. One tracked table
    per ``CHUNK_BYTES`` of tar; finished chunks are skipped (resumable)."""
    import multiprocessing as mp
    import queue as queue_mod
    from collections import defaultdict

    n, total = n_chunks(field), moa.TAR_BYTES[field_number(field)]
    todo = [k for k in (range(n) if chunks is None else chunks) if not _chunk_done(field, k, n)]
    if not todo:
        print(f"{field}: all {n} pre-screen chunks done")
        return []
    reader = reader_for(field, source)
    segs = [
        (k, a, b)
        for k in todo
        for a, b in moa_stream.segments(
            k * CHUNK_BYTES, min((k + 1) * CHUNK_BYTES, total), SEG_BYTES
        )
    ]
    seg_left = {k: sum(1 for s in segs if s[0] == k) for k in todo}
    seg_pos = {s: i for i, s in enumerate(segs)}
    seg_sha: dict = {}
    prefetch = max(1, conns // procs)
    print(f"{field}: {len(todo)} of {n} chunks, {len(segs)} segments of {SEG_BYTES >> 20} MB, "
          f"{procs} processes × {prefetch} prefetching connections, source "
          f"{reader.url or reader.path}", flush=True)  # fmt: skip
    ctx = mp.get_context("fork") if "fork" in mp.get_all_start_methods() else mp.get_context()
    queue = ctx.Queue(maxsize=4 * procs)
    workers = [
        ctx.Process(
            target=_stream_worker,
            args=(i, segs[i::procs], total, reader.url, reader.path, reader.total, prefetch, queue),
            daemon=True,
        )
        for i in range(procs)
    ]
    for w in workers:
        w.start()
    rows = defaultdict(list)
    t_chunk = {k: None for k in todo}
    meter_chunk, bytes_chunk = {}, defaultdict(int)
    meter = moa_stream.CpuMeter()
    retries = [0] * procs
    written = []
    t0 = t_log = time.time()
    for k in todo[:1]:
        t_chunk[k], meter_chunk[k] = t0, moa_stream.CpuMeter()
    n_items = n_items_log = bytes_log = 0
    finished = 0
    try:
        while finished < procs:
            try:
                msg = queue.get(timeout=60)
            except queue_mod.Empty:
                dead = [w.exitcode for w in workers if not w.is_alive() and w.exitcode]
                if dead:
                    raise SystemExit(
                        f"{field} pre-screen: a worker died (exit codes {dead})"
                    ) from None
                continue
            if msg is None:
                finished += 1
                continue
            if msg[0] == "error":
                raise SystemExit(f"{field} pre-screen failed: {msg[1]}")
            s, out, nret, sha = msg
            seg_sha[s] = sha
            k = s[0]
            if t_chunk[k] is None:  # a later chunk started (rough: its first finished segment)
                t_chunk[k], meter_chunk[k] = time.time(), moa_stream.CpuMeter()
            retries[seg_pos[s] % procs] = nret
            rows[k].extend(out)
            n_items += len(out)
            bytes_chunk[k] += s[2] - s[1]
            bytes_log += s[2] - s[1]
            seg_left[k] -= 1
            if seg_left[k] == 0:
                wall = time.time() - t_chunk[k]
                busy = meter_chunk[k].busy()
                st = {
                    "wall_time_s": round(wall, 1),
                    "items_per_s": round(len(rows[k]) / wall, 1),
                    "mb_per_s": round(bytes_chunk[k] / 1e6 / wall, 1),
                    "cpu_percent": None if busy is None else round(busy, 1),
                    "procs": procs,
                    "conns": procs * prefetch,
                    "segment_bytes": SEG_BYTES,
                    "segment_sha256": [seg_sha[x] for x in segs if x[0] == k],
                }
                written.append(_write_prescreen_chunk(field, k, n, rows.pop(k), st))
                shown = {k2: v for k2, v in st.items() if k2 != "segment_sha256"}
                print(f"{field} chunk {k + 1}/{n}: {shown}", flush=True)
            now = time.time()
            if now - t_log >= log_every:
                busy = meter.busy()
                print(
                    f"  {n_items} light curves, {(n_items - n_items_log) / (now - t_log):.0f}/s, "
                    f"{bytes_log / 1e6 / (now - t_log):.0f} MB/s, CPU "
                    f"{busy if busy is None else round(busy)} %, retries {sum(retries)}",
                    flush=True,
                )
                t_log, n_items_log, bytes_log = now, n_items, 0
    finally:
        for w in workers:
            if w.is_alive() and finished < procs:
                w.terminate()
            w.join()
    print(f"{field}: {n_items} light curves in {time.time() - t0:.0f} s", flush=True)
    return written


def merge_prescreen(field: str) -> Path:
    """Join the tracked pre-screen chunks of a field into the table the later stages read.
    Refused unless every chunk exists with the current ``Params`` and the member count equals the
    field's Cut-0 count in the metadata (no member lost or duplicated at a range boundary)."""
    n = n_chunks(field)
    parts = []
    for k in range(n):
        p = prescreen_chunk_path(field, k, n)
        if not p.exists():
            raise SystemExit(f"missing pre-screen chunk {k + 1}/{n} of {field}: {p}")
        tab = Table.read(p, format="ascii.ecsv")
        if tab.meta.get("params") != json.dumps(asdict(P)):
            raise SystemExit(f"pre-screen chunk {k + 1}/{n} of {field} has other Params; rerun")
        parts.append(tab)
    tab = vstack(parts, metadata_conflicts="silent")
    tab.sort("event_id")
    n_members = sum(int(t.meta["n_members"]) for t in parts)
    expected = moa.CUT0_PER_FIELD[field_number(field)]
    if n_members != expected:
        raise SystemExit(f"{field}: {n_members} light curves streamed, metadata has {expected}")
    if len(set(tab["event_id"])) != len(tab):
        raise SystemExit(f"{field}: duplicate light curves across chunks")
    keys = ("n_members", "n_errors", "n_shape", "n_quiet", "wall_time_s")
    tab.meta = {
        "provenance": schema.Provenance.DERIVED.value,
        "source": f"scripts/w3_moa.py merge-prescreen: {n} tracked chunk tables of {field}",
        "params": json.dumps(asdict(P)),
        "tracked_subset": True,
        **{k: sum(float(t.meta.get(k, 0)) for t in parts) for k in keys},
        "quiet_track_mod": QUIET_TRACK_MOD,
        "quiet_chi2_hist": np.sum([t.meta["quiet_chi2_hist"] for t in parts], axis=0).tolist(),
        "range_sha256": field_range_digest(field) or "",
    }
    tab.meta["n_members"] = int(tab.meta["n_members"])
    path = out_dir() / f"prescreen_{field}.ecsv"
    tab.write(path, overwrite=True)
    print(
        f"wrote {path}: {len(tab)} tracked rows of {n_members} light curves; "
        f"{int(tab.meta['n_shape'])} shape passes, {len(passes(tab))} not at a shared epoch"
    )
    return path


def read_prescreen(field: str | None = None) -> Table:
    return Table.read(out_dir() / f"prescreen_{field or FIELD}.ecsv")


def n_light_curves(pre: Table) -> int:
    """Light curves screened (a merged table holds a tracked subset; its meta has the count)."""
    return int(pre.meta.get("n_members", len(pre)))


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


@functools.lru_cache(maxsize=2)
def _local_index(path: Path) -> dict[str, tuple[int, int]]:
    """Member offsets of a local field tar (header walk only)."""
    return moa.MoaField(tar_path=path).index()


def load_arrays(ids, aux: bool = False, pre: Table | None = None, field: str | None = None) -> dict:
    """``event_id -> (t, f, sf)`` (plus a dict of ``AUX`` columns when ``aux``), included epochs
    only. Members with a tar offset in the pre-screen table are read by byte range (local tar or
    HTTP); others (vetting neighbours) from the per-object archive files."""
    pre = read_prescreen(field) if pre is None else pre
    loc = {}
    if "offset" in pre.colnames:
        loc = {
            str(e): (int(o), int(s))
            for e, o, s in zip(pre["event_id"], pre["offset"], pre["size"], strict=True)
        }
    ids = list(dict.fromkeys(map(str, ids)))
    reader = reader_for(field)
    if reader.path is not None and any(e not in loc for e in ids):
        loc = {**_local_index(reader.path), **loc}  # the pinned local tar has every member
    raws = moa_stream.fetch_members(reader, [(e, *loc[e]) for e in ids if e in loc])
    missing = [e for e in ids if e not in loc]
    if missing:
        raws.update(moa_stream.fetch_objects(missing))
    want = (*PRESCREEN_COLUMNS, *AUX) if aux else PRESCREEN_COLUMNS
    out = {}
    for eid in ids:
        cols = moa.parse_lightcurve(raws[eid], columns=want)
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
    """Light curves passing the shape cuts (``prescreen_pass``) whose deficit epoch is not shared
    by improbably many other objects of the field or chip (``epoch_artefact``)."""
    sel = pre[no_error(pre) & prescreen_pass(pre["z_min"], pre["s_min"], pre["z_min2"])]
    pop, chips = deficit_population(pre), chip_populations(pre)
    sh = [epoch_artefact(pop, chips, r, str(r["event_id"])) for r in sel]
    for k in ("field_n", "field_p", "chip_n", "chip_p"):
        sel[f"shared_{k}"] = [x[k] for x in sh]
    keep = np.array([x["p"] >= COINC_P for x in sh], bool)
    return sel[keep]


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
    pre = read_prescreen()
    sel = passes(pre)
    if chunk is not None:
        sel = sel[chunk[0] :: chunk[1]]
    if limit:
        sel = sel[:limit]
    ev = field.events()
    pos = {r["event_id"]: (float(r["ra"]), float(r["dec"])) for r in ev}
    arrays = load_arrays(sel["event_id"], pre=pre)
    jobs = []
    for r in sel:
        eid = str(r["event_id"])
        t, f, sf = arrays[eid]
        scan = {k: float(r[k]) for k in SCAN_KEYS}
        d = {"event_id": eid, "ra": pos[eid][0], "dec": pos[eid][1]}
        jobs.append((d, t, f, sf, scan))
    print(f"fitting {len(jobs)} pre-screen passes of {n_light_curves(pre)}", flush=True)
    t1 = time.time()
    rows = []
    with _pool(procs) as pool:
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
        n_prescreen=n_light_curves(pre),
        n_passes=len(passes(pre)),
        chunk="" if chunk is None else f"{chunk[0] + 1}/{chunk[1]}",
    )
    # a chunk run never replaces the (possibly merged) table that `vet` reads
    path = out_dir() / (f"fits_{FIELD}.ecsv" if chunk is None else chunk_name(chunk))
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


def fit_chunk_problem(tab: Table, ids: list, k: int, n: int) -> str | None:
    """Why the tracked fit chunk ``k`` of ``n`` cannot be used (None: it can): fitted with other
    ``Params``, or not exactly its share ``ids[k::n]`` of the current pre-screen passes."""
    if tab.meta.get("fit_params") != json.dumps(asdict(w3.P)) or tab.meta.get(
        "params"
    ) != json.dumps(asdict(P)):
        return "was fitted with other Params; refit it"
    if sorted(tab["event_id"]) != sorted(ids[k::n]):
        return "does not hold exactly its pre-screen passes"
    return None


def merge_chunks(n: int) -> Path:
    """Join the tracked chunk tables 1..n of n into the table `vet` reads. Refused unless every
    chunk is present, was fitted with the current pre-screen and fit ``Params`` and holds exactly
    its own passes of the current pre-screen (deterministic, recomputed in each session)."""
    pre = read_prescreen()
    ids = list(passes(pre)["event_id"])
    parts = []
    for k in range(n):
        path = results_dir() / f"{chunk_name((k, n))}.gz"
        if not path.exists():
            raise SystemExit(f"missing chunk {k + 1}/{n}: {path}")
        tab = Table.read(path, format="ascii.ecsv")
        problem = fit_chunk_problem(tab, ids, k, n)
        if problem:
            raise SystemExit(f"chunk {k + 1}/{n} {problem}")
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
BASELINE_CHI2 = 2.0  # D-057/D-062 fixed threshold; used only when no field calibration is set
BASELINE_Q = 0.95  # ASSUMPTION (D-068): threshold = this quantile of the field's quiet χ²/dof
NEIGHBOUR_S = 5.0  # ASSUMPTION: |S| of a neighbour's notch over the same window = shared feature
MIN_FEATURE_NIGHTS = 3  # ASSUMPTION: nights with epochs inside the exotic feature
COINC_Z = 5.0  # deficits with z_min < −5 form the population for the shared-epoch test
COINC_P = 1e-3  # ASSUMPTION: Poisson probability below which a shared epoch is a frame systematic

DOMAIN_U0_MAX = 2.0  # the injection and limit domain: umbra crossings with u0 < 2 (D-068)
FS_MAX_FACTOR = 3.0  # ASSUMPTION: fitted source flux ≤ 3 × the object's reference flux
FS_REF_DEFAULT_MAG = 14.2  # ASSUMPTION: no reference magnitude → the bright end of the injections
REF_MATCH_ARCSEC = 1.0  # ASSUMPTION: Gaia DR3 counterpart radius for the reference flux
BRACKET_MIN_EPOCHS = 20  # ASSUMPTION: baseline epochs required before ingress and after egress

_POP: tuple | None = None  # (field, per-chip) deficit populations (passed to the workers)
_BASELINE: float | None = None  # calibrated variable-baseline threshold of the field


def deficit_population(pre: Table, chip: int | None = None) -> dict:
    """Centres of every significant deficit box of the field (or of one chip), per box width and
    pooled over widths (key ``"all"``) (``derived``)."""
    ok = no_error(pre) & (np.asarray(pre["z_min"], float) < -COINC_Z)
    ids = np.asarray(pre["event_id"], str)
    if chip is not None:
        ok &= np.array([moa.parse_event_id(e)[1] == chip for e in ids])
    w = np.asarray(pre["width"], float)[ok]
    c = 0.5 * (np.asarray(pre["t_lo"], float) + np.asarray(pre["t_hi"], float))[ok]
    ids = ids[ok]
    out = {}
    for key in [*np.unique(w).tolist(), "all"]:
        m = np.ones(w.size, bool) if key == "all" else w == key
        order = np.argsort(c[m])
        out[key] = (c[m][order], ids[m][order])
    return out


def calibrate_baseline(pre: Table, q: float | None = None) -> dict:
    """Variable-baseline threshold of a field (``derived``, D-068): the ``BASELINE_Q`` quantile
    of the whole-light-curve χ²/dof about a constant (errors × the point-to-point scale, the
    vetting statistic) over the field's quiet light curves, the injection carriers. Difference
    photometry has red noise, so a fixed χ²/dof > 2 removed 35 % of quiet gb22 carriers; the
    quantile removes a fixed small fraction of them and every baseline noisier than that."""
    q = BASELINE_Q if q is None else q
    hist = pre.meta.get("quiet_chi2_hist")
    if hist is not None:  # every quiet light curve of the field (the rows are a sample)
        h = np.asarray(hist, float)
        n = int(h.sum())
        cdf = np.concatenate([[0.0], np.cumsum(h)]) / max(n, 1)
        thr, med = (float(10 ** np.interp(x, cdf, np.log10(CHI2_BINS))) for x in (q, 0.5))
        above2 = float(1.0 - np.interp(math.log10(2.0), np.log10(CHI2_BINS), cdf))
    else:
        c = np.asarray(pre["chi2_const"], float)[is_quiet(pre)]
        c = c[np.isfinite(c)]
        n = int(c.size)
        thr, med, above2 = (
            (float(np.quantile(c, q)), float(np.median(c)), float(np.mean(c > 2.0)))
            if n
            else (np.nan,) * 3
        )
    if n < 50:
        raise SystemExit(f"{FIELD}: {n} quiet light curves; too few to calibrate the baseline")
    return {"threshold": thr, "q": q, "n_quiet": n, "frac_above_2": above2, "median": med}


CHI2_BINS = np.logspace(-1.0, 3.0, 801)  # χ²/dof histogram edges (0.5 % wide) for the calibration


def chi2_histogram(c: np.ndarray) -> list[int]:
    """Counts of ``c`` in ``CHI2_BINS`` (values outside are clipped into the end bins)."""
    c = np.clip(c[np.isfinite(c)], CHI2_BINS[0], CHI2_BINS[-1] * 0.999999)
    return np.histogram(c, CHI2_BINS)[0].astype(int).tolist()


def baseline_threshold() -> float:
    return BASELINE_CHI2 if _BASELINE is None else _BASELINE


def set_field_context(pre: Table) -> None:
    """Shared-epoch populations and the baseline calibration of the field, for vet and inject."""
    global _POP, _BASELINE
    _POP = (deficit_population(pre), chip_populations(pre))
    _BASELINE = calibrate_baseline(pre)["threshold"] if "chi2_const" in pre.colnames else None


def chip_populations(pre: Table) -> dict:
    return {chip: deficit_population(pre, chip) for chip in range(1, 11)}


def shared_epoch(
    pop: dict, width: float, centre: float, self_id: str = "", pooled: bool = False
) -> dict:
    """How many other objects have a deficit box of the same width (``pooled``: any width) centred
    within max(W/2, 1 d) of this one, against the Poisson expectation for boxes spread uniformly
    over the survey. Many unrelated stars dimming together is an artefact of the frame or of the
    chip (bad seeing, clouds, subtraction, a detector problem), not lensing of one star."""
    key = "all" if pooled else width
    if key not in pop:
        return {"n": 0, "expected": 0.0, "p": 1.0}
    cs, ids = pop[key]
    tol = max(width / 2.0, 1.0)
    lo, hi = np.searchsorted(cs, centre - tol), np.searchsorted(cs, centre + tol)
    n = int(np.sum(ids[lo:hi] != self_id))
    lam = max(len(cs) - 1, 0) * 2.0 * tol / (moa.T_END - moa.T_START)
    from scipy.stats import poisson

    return {"n": n, "expected": float(lam), "p": float(poisson.sf(n - 1, lam)) if n else 1.0}


def epoch_artefact(field_pop: dict, chip_pops: dict, scan: dict, event_id: str) -> dict:
    """Shared-epoch probabilities in the whole field (same width) and on the object's chip (any
    width); the smaller one decides."""
    c = 0.5 * (scan["t_lo"] + scan["t_hi"])
    fld = shared_epoch(field_pop, scan["width"], c, event_id)
    chip = moa.parse_event_id(event_id)[1]
    chp = shared_epoch(chip_pops.get(chip, {}), scan["width"], c, event_id, pooled=True)
    return {"field_n": fld["n"], "field_p": fld["p"], "chip_n": chp["n"], "chip_p": chp["p"],
            "p": min(fld["p"], chp["p"])}  # fmt: skip


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


def reference_flux(ev: dict) -> tuple[float, str]:
    """Counts of the object's reference flux for the source-flux bound of ``exotic_in_domain``:
    DoPHOT magnitude from the release where present, else the Gaia DR3 RP magnitude of a
    counterpart within ``REF_MATCH_ARCSEC`` (``ev["ref_mag"]``; MOA-Red ≈ RP, ASSUMPTION), else
    ``FS_REF_DEFAULT_MAG`` (the brightest injected source)."""
    chip = moa.parse_event_id(str(ev["event_id"]))[1] if "-R-" in str(ev["event_id"]) else 1
    for key, label in (("dophot_mag", "DoPHOT"), ("ref_mag", "Gaia DR3 RP")):
        m = ev.get(key)
        if m is not None and np.isfinite(float(m)):
            return float(moa.mag_to_counts(float(m), chip)), f"{label} {float(m):.2f}"
    return float(moa.mag_to_counts(FS_REF_DEFAULT_MAG, chip)), f"default {FS_REF_DEFAULT_MAG}"


def exotic_domain(res: dict, ev: dict) -> dict:
    """Per exotic model: inside the domain the limit is defined on (u0 < ``DOMAIN_U0_MAX``) and with
    a physical source flux (0 < f_s ≤ ``FS_MAX_FACTOR`` × the reference flux)? A far-field fit
    (u0 ≫ 1) with cancelling giant f_s and f_b can mimic any smooth dip (gb20-R-4-0-49379)."""
    f_ref, label = reference_flux(ev)
    out = {}
    for m in w3.EXOTIC:
        if m not in res:
            continue
        u0, fs = float(res[m].get("u0", np.nan)), float(res[m].get("fs", np.nan))
        ok = bool(u0 < DOMAIN_U0_MAX and 0.0 < fs <= FS_MAX_FACTOR * f_ref)
        out[m] = (ok, f"u0 {u0:.2f}, f_s {fs:.3g} vs {FS_MAX_FACTOR:g} × {f_ref:.3g} ({label})")
    return out


def gauss_dip(t, tc: float, sigma: float) -> np.ndarray:
    return np.exp(-0.5 * ((np.asarray(t, float) - tc) / sigma) ** 2)


def fit_smooth_dip(lc: w3.LightCurve, lo: float, hi: float) -> dict:
    """Ordinary smooth dimming: baseline − A × Gaussian(tc, σ), A ≥ 0 and baseline from a weighted
    linear solve; Nelder–Mead over (tc, log σ) (slow red-giant variability, a slow subtraction
    residual). k = 4."""
    from scipy.optimize import minimize

    span = float(lc.t.max() - lc.t.min())

    def chi2(x):
        tc, ls = x
        if not (-1.0 <= ls <= math.log10(span)):
            return 1e30
        g = gauss_dip(lc.t, tc, 10.0**ls)
        coef, c2 = w3.linear_fluxes_n(-g[None, :], lc.f, lc.w)
        return c2 if coef[0] >= 0 else 1e30

    tc0, w0 = 0.5 * (lo + hi), max(hi - lo, 0.5)
    starts = [
        (tc0 + dt * w0, math.log10(w0 * g))
        for dt in (-0.25, 0.0, 0.25)
        for g in (0.1, 0.25, 0.5, 1.0)
    ]
    vals = [chi2(np.array(x)) for x in starts]
    best = None
    for i in np.argsort(vals)[:3]:
        r = minimize(chi2, np.array(starts[i]), method="Nelder-Mead",
                     options={"maxfev": 1000, "xatol": 1e-5})  # fmt: skip
        if best is None or r.fun < best.fun:
            best = r
    k = 4
    return {"chi2": float(best.fun), "k": k, "bic": float(best.fun) + k * math.log(lc.t.size),
            "tc": float(best.x[0]), "sigma": float(10.0 ** best.x[1])}  # fmt: skip


def bracketing(t, lo: float, hi: float) -> tuple[int, int]:
    """Epochs before the exotic feature's ingress and after its egress, inside the data."""
    t = np.asarray(t, float)
    return int(np.sum(t < lo)), int(np.sum(t > hi))


def step_ramp(t, ts: float) -> np.ndarray:
    """Design rows of the ordinary step model: H(t − ts) and (t − ts) H(t − ts) (days)."""
    x = np.asarray(t, float) - ts
    h = (x >= 0).astype(float)
    return np.vstack([h, x * h])


def fit_step_ramp(lc: w3.LightCurve) -> dict:
    """Ordinary step: baseline, a change of level at ts and a linear ramp after it (a reference or
    photometric-scale change, a secular change of the star; gb19-R-4-4-31159). Linear in the
    three fluxes; ts scanned over every 3rd night, then every epoch near the best. k = 4."""
    t = np.asarray(lc.t, float)
    nt = nights(t)
    first = np.r_[0, np.flatnonzero(np.diff(nt)) + 1]  # first epoch of each night
    cand = t[first[1::3]] if first.size > 3 else t[first]

    def chi2(ts):
        if not (t[0] < ts <= t[-1]):
            return np.inf, None
        coef, c2 = w3.linear_fluxes_n(step_ramp(t, ts), lc.f, lc.w)
        return c2, coef

    vals = [chi2(x)[0] for x in cand]
    best = float(cand[int(np.argmin(vals))])
    near = t[(t > best - 10.0) & (t < best + 10.0)]
    vals2 = [chi2(x)[0] for x in near] if near.size else [np.inf]
    if near.size and min(vals2) < min(vals):
        best = float(near[int(np.argmin(vals2))])
    c2, coef = chi2(best)
    k = 4
    return {"chi2": float(c2), "k": k, "bic": float(c2) + k * math.log(t.size), "ts": best,
            "step": float(coef[0]), "ramp_per_day": float(coef[1])}  # fmt: skip


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

    allowed = list(w3.EXOTIC)

    def summary():
        ex = min(allowed, key=lambda m: res[m]["bic"])
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
    # --- the fit must lie in the limit's domain with a physical source flux; a flag is re-judged
    # on its in-domain exotic models only
    dom = exotic_domain(res, ev)
    out["domain"] = {m: v[1] for m, v in dom.items()}
    allowed = [m for m in w3.EXOTIC if dom.get(m, (False,))[0]]
    if not allowed:
        add("exotic_in_domain", False, f"no in-domain exotic fit: {out['domain']}")
        return record()
    ex, ordinary, best_o, d0 = summary()
    out.update(exotic=ex, dbic_all=float(d0), best_ordinary=best_o)
    lo, hi = feature_window(res[ex], ex)
    out["feature_window"] = (lo, hi)
    if not add("exotic_in_domain", d0 < P.flag_dbic, f"{ex} ΔBIC {d0:.1f}; {dom[ex][1]}"):
        return record()
    # --- the feature must be bracketed by baseline: a one-sided step or secular change is not an
    # umbra crossing (gb19-R-4-4-31159)
    n_bef, n_aft = bracketing(lc.t, lo, hi)
    out["bracket"] = (n_bef, n_aft)
    if not add(
        "feature_bracketed",
        min(n_bef, n_aft) >= BRACKET_MIN_EPOCHS,
        f"{n_bef} epochs before ingress, {n_aft} after egress (need {BRACKET_MIN_EPOCHS} each)",
    ):
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
        not (chi_out > baseline_threshold()),
        f"baseline χ²/dof {chi_out:.2f} > {baseline_threshold():.2f}? ({int(far.sum())} pts, "
        "errors × point-to-point scale)",
    ):
        return record()
    s_n = [(nid, notch_in_window(nt, nf, nsf, lo, hi)) for nid, nt, nf, nsf in neigh]
    hit = [(nid, round(v, 1)) for nid, v in s_n if abs(v) > NEIGHBOUR_S]
    out["neighbours"] = s_n
    if not add("neighbour_shares_feature", not hit, f"{len(s_n)} neighbours; |S| > 5: {hit}"):
        return record()
    # --- an eclipse or occultation: a trapezoidal dip with no caustic spikes
    ecl = fit_eclipse(lc, lo, hi)
    out["eclipse"] = ecl
    d_ecl = res[ex]["bic"] - min(ordinary[best_o]["bic"], ecl["bic"])
    if not add("eclipse_dip", d_ecl < P.flag_dbic, f"ΔBIC {d_ecl:.1f} vs a trapezoidal dip"):
        return record()
    # --- slow smooth dimming (red-giant variability, a slow subtraction residual)
    sd = fit_smooth_dip(lc, lo, hi)
    out["smooth_dip"] = sd
    d_sd = res[ex]["bic"] - min(ordinary[best_o]["bic"], sd["bic"])
    if not add(
        "smooth_dip",
        d_sd < P.flag_dbic,
        f"ΔBIC {d_sd:.1f} vs a Gaussian dip (σ {sd['sigma']:.1f} d)",
    ):
        return record()
    st = fit_step_ramp(lc)
    out["step_ramp"] = st
    d_st = res[ex]["bic"] - min(ordinary[best_o]["bic"], st["bic"])
    if not add(
        "step_ramp",
        d_st < P.flag_dbic,
        f"ΔBIC {d_st:.1f} vs a step at {st['ts']:.1f} ({st['step']:.0f} counts, "
        f"ramp {st['ramp_per_day']:.2f}/d)",
    ):
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
    fe = w3.model_flux(ex, lc, res[ex])
    fo = w3.model_flux(best_o, lc, ordinary[best_o])
    cov["n_nights"] = int(np.unique(nights(lc.t[np.abs(fe - fo) > 3.0 * np.median(lc.sf)])).size)
    out["feature"] = cov
    # MOA takes several exposures a night: the D-057 "≥ 3 epochs" becomes ≥ 3 nights (ASSUMPTION)
    ok_cov = cov["n_nights"] >= MIN_FEATURE_NIGHTS and cov["dchi2_in"] < P.flag_dbic
    if not add(
        "exotic_feature_sampled",
        ok_cov,
        f"{cov['n_epochs']} epochs on {cov['n_nights']} nights inside, Δχ² in "
        f"{cov['dchi2_in']:.1f} / out {cov['dchi2_out']:.1f}",
    ):
        return record()
    n_drop = w3.jackknife_n_drop(cov["n_epochs"])
    if n_drop:
        jk = w3.jackknife_worst_epochs(lc, ordinary[best_o], res[ex], best_o, ex, n_drop)
        if not add("jackknife_epochs", jk[-1] < P.flag_dbic, f"ΔBIC {[round(v, 1) for v in jk]}"):
            return record()
    jn = jackknife_nights(lc, res, ex, best_o, ecl)
    out["jackknife_nights"] = jn
    if jn["n_dropped"] and not add(
        "jackknife_nights",
        jn["dbic"] < P.flag_dbic,
        f"ΔBIC {jn['dbic']:.1f} after dropping the {jn['n_dropped']} most influential night(s) "
        f"of {jn['n_nights']} in the feature",
    ):
        return record()
    d2e, _two = w3.two_events_dbic(lc, ordinary["PSPL"], res[ex], "PSPL", ex)
    if d2e is not None and not add("two_unrelated_events", d2e < P.flag_dbic, f"ΔBIC {d2e:.1f}"):
        return record()
    out["complete"] = bool(binary_lens and w3.have_mm())
    out["survives"] = True
    return record()


def trapezoid(t, tc: float, dur: float, fin: float) -> np.ndarray:
    """Dip shape in [0, 1]: flat bottom for |t − tc| < dur/2 (1 − fin), linear ingress/egress."""
    x = np.abs(np.asarray(t, float) - tc)
    half = 0.5 * dur
    ramp = max(fin, 1e-3) * half
    return np.clip((half - x) / ramp, 0.0, 1.0)


def fit_eclipse(lc: w3.LightCurve, lo: float, hi: float) -> dict:
    """Ordinary dip model: baseline − depth × trapezoid (eclipse, occultation, dipper), depth ≥ 0
    and baseline from a weighted linear solve; Nelder–Mead over (tc, log duration, ingress)."""
    from scipy.optimize import minimize

    def chi2(x):
        tc, ld, fin = x
        if not (0.0 <= fin <= 1.0 and -1.5 <= ld <= 3.5):
            return 1e30
        s = trapezoid(lc.t, tc, 10.0**ld, fin)
        coef, c2 = w3.linear_fluxes_n(-s[None, :], lc.f, lc.w)
        return c2 if coef[0] >= 0 else 1e30

    tc0, dur0 = 0.5 * (lo + hi), max(hi - lo, 0.2)
    starts = [
        (tc0 + dt * dur0, math.log10(dur0 * g), fin)
        for dt in (-0.2, 0.0, 0.2)
        for g in (0.3, 0.6, 1.0)
        for fin in (0.05, 0.3, 0.7)
    ]
    vals = [chi2(np.array(s)) for s in starts]
    best = None
    for i in np.argsort(vals)[:3]:
        r = minimize(
            chi2, np.array(starts[i]), method="Nelder-Mead", options={"maxfev": 1500, "xatol": 1e-5}
        )
        if best is None or r.fun < best.fun:
            best = r
    n = lc.t.size
    k = 3 + 2
    tc, ld, fin = best.x
    return {"chi2": float(best.fun), "k": k, "bic": float(best.fun) + k * math.log(n),
            "tc": float(tc), "duration": float(10.0**ld), "ingress": float(fin)}  # fmt: skip


def eclipse_flux(lc: w3.LightCurve, ecl: dict) -> np.ndarray:
    """Model flux of a ``fit_eclipse`` result on ``lc`` (baseline and depth re-solved)."""
    sh = trapezoid(lc.t, ecl["tc"], ecl["duration"], ecl["ingress"])
    coef, _ = w3.linear_fluxes_n(-sh[None, :], lc.f, lc.w)
    return coef[1] - coef[0] * sh


def jackknife_nights(lc: w3.LightCurve, res: dict, ex: str, best_o: str, ecl: dict) -> dict:
    """Drop whole nights (MOA takes several exposures a night, so one bad night is several bad
    epochs) in order of their contribution to the exotic preference, at most 2 and never fewer than
    2 nights left inside the feature (ASSUMPTION). The nights are ranked twice, against the best
    single-lens model and against the eclipse (trapezoid) model, because the preference over each
    can rest on different nights (gb7-R-8-6-94052: one bright night at the predicted ingress spike
    carried the whole preference over the trapezoid). Returns the larger ΔBIC (exotic − best
    ordinary, single lens or eclipse, refitted on the subset) of the two rankings."""
    fe = w3.model_flux(ex, lc, res[ex])
    fo = w3.model_flux(best_o, lc, res[best_o])
    nt = nights(lc.t)
    inside = np.abs(fe - fo) > 3.0 * np.median(lc.sf)
    feat_nights = np.unique(nt[inside])
    n_drop = int(min(2, max(0, feat_nights.size - 2)))
    out = {"n_nights": int(feat_nights.size), "n_dropped": n_drop, "dbic": np.nan}
    if not n_drop:
        return out
    lo, hi = float(lc.t[inside].min()), float(lc.t[inside].max())
    fecl = eclipse_flux(lc, ecl)
    best, dropped = -np.inf, []
    for alt in (fo, fecl):
        dchi = (lc.f - fe) ** 2 * lc.w - (lc.f - alt) ** 2 * lc.w
        per = {k: float(dchi[nt == k].sum()) for k in np.unique(nt)}
        worst = sorted(per, key=per.get)[:n_drop]
        sub = lc.subset(~np.isin(nt, worst))
        d = refit_both(sub, res, ex)
        e2 = fit_eclipse(sub, lo, hi)
        e_bic = w3.optimise(ex, sub, [_start(ex, res[ex])])["bic"]
        dd = float(max(d, e_bic - e2["bic"]))
        if dd > best:
            best, dropped = dd, [int(x) for x in worst]
    out["dbic"] = best
    out["dropped"] = dropped
    return out


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


def gaia_rp_mags(ra, dec) -> tuple[np.ndarray, str]:
    """Gaia DR3 RP magnitude of the nearest counterpart within ``REF_MATCH_ARCSEC`` (NaN: none), one
    batched CDS XMatch (``observed``); a failed query gives NaN everywhere (then the default
    reference flux, the lenient bound) and says so."""
    out = np.full(len(ra), np.nan)
    if not len(ra):
        return out, "no flags"
    import astropy.units as u
    from astroquery.xmatch import XMatch

    pos = Table({"ra": np.asarray(ra, float), "dec": np.asarray(dec, float)})
    pos["idx"] = np.arange(len(pos))
    try:
        m = XMatch.query(cat1=pos, cat2="vizier:I/355/gaiadr3",
                         max_distance=REF_MATCH_ARCSEC * u.arcsec, colRA1="ra",
                         colDec1="dec")  # fmt: skip
    except Exception as exc:  # noqa: BLE001
        return out, f"query failed: {exc!r}"[:200]
    m.sort("angDist")
    for row in m[::-1]:  # nearest written last
        if np.isfinite(float(row["RPmag"])):
            out[int(row["idx"])] = float(row["RPmag"])
    return (
        out,
        f"vizier:I/355/gaiadr3 within {REF_MATCH_ARCSEC}″: {int(np.isfinite(out).sum())} matched",
    )


def neighbours_of(ev: Table, eid: str, radius_px: float = NEIGHBOUR_PX) -> list[str]:
    """Cut-0 objects of the same chip and subframe within ``radius_px`` pixels."""
    r = ev[ev["event_id"] == eid][0]
    same = (ev["chip"] == r["chip"]) & (ev["subframe"] == r["subframe"])
    d = np.hypot(ev["x"] - r["x"], ev["y"] - r["y"])
    return [str(x) for x in ev["event_id"][same & (d < radius_px) & (ev["event_id"] != eid)]]


def check_complete_fits(fits: Table) -> None:
    """Vet only a fit table of every pre-screen pass: a chunk table (``fit --chunk`` writes one
    before `merge-chunks`) or a ``--limit`` run would turn a partial screen into a null result."""
    if fits.meta.get("chunk", ""):
        raise SystemExit(f"fit table holds chunk {fits.meta['chunk']} only; run merge-chunks")
    pre = read_prescreen()
    if sorted(map(str, fits["event_id"])) != sorted(map(str, passes(pre)["event_id"])):
        raise SystemExit("fit table does not hold exactly the pre-screen passes; refit")


def run_vet(procs: int) -> Path:
    field = moa.MoaField(FIELD)
    set_field_context(read_prescreen())
    fits = Table.read(out_dir() / f"fits_{FIELD}.ecsv")
    check_complete_fits(fits)
    ok = no_error(fits)
    flags = fits[ok & (np.asarray(fits["dbic_min"], float) < P.flag_dbic)]
    ev = field.events()
    pos = {r["event_id"]: r for r in ev}
    neigh = {str(e): neighbours_of(ev, str(e)) for e in flags["event_id"]}
    ids = set(map(str, flags["event_id"])) | {n for v in neigh.values() for n in v}
    arrays = load_arrays(ids, aux=True)
    fl_ids = [str(e) for e in flags["event_id"]]
    rp, rp_note = gaia_rp_mags(
        [float(pos[e]["ra"]) for e in fl_ids], [float(pos[e]["dec"]) for e in fl_ids]
    )
    jobs = []
    for i, r in enumerate(flags):
        eid = str(r["event_id"])
        t, f, sf, ax = arrays[eid]
        scan = {k: float(r[k]) for k in SCAN_KEYS}
        nb = [(n, *arrays[n][:3]) for n in neigh[eid] if n in arrays]
        d = {"event_id": eid, "ra": float(pos[eid]["ra"]), "dec": float(pos[eid]["dec"]),
             "dophot_mag": float(pos[eid]["dophot_magnitude"]),
             "ref_mag": float(rp[i])}  # fmt: skip
        jobs.append((d, t, f, sf, ax, scan, nb, True, res_from_row(r)))
    print(f"vetting {len(jobs)} flags", flush=True)
    t1 = time.time()
    out = []
    with _pool(procs) as pool:
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
        "fit_errors": sorted(map(str, fits["event_id"][~ok])),  # unscreened: no null limit
        "n_flags": len(flags),
        "wall_time_s": time.time() - t1,
        "variable_xmatch": var,
        "gaia_rp_xmatch": rp_note,
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
    arrays = load_arrays([o["event_id"] for o in vet])
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
            lab = "PSPL (PAR not drawn)" if m == "PAR" else m
            ax.plot(tt - lo, w3.model_flux(mm, fine, r), color=c, lw=0.8, label=lab)
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
INJ_IS = (14.2, 21.4)  # source MOA-Red magnitudes (Koshimoto et al. 2023's range)
LF_SLOPE = 0.319  # dex/mag: Gaia DR3 RP counts in gb22, RP 13–17.5 (ASSUMPTION: extends to 21.4)


def quiet_carriers(pre: Table, n: int, seed: int) -> list[str]:
    """Light curves with no significant notch either way (``|z| < 4``) and ≥ 1,000 epochs: stand-ins
    for the difference light curve of a constant star (ASSUMPTION; constant stars are not Cut-0
    objects, so their own light curves are not in the release)."""
    ids = np.asarray(pre["event_id"][is_quiet(pre)], str)
    rng = np.random.default_rng(seed)
    return list(rng.choice(ids, size=min(n, ids.size), replace=False))


def sample_magnitudes(rng, size: int, sampling: str = "lf") -> np.ndarray:
    """Source magnitudes in ``INJ_IS``: ``"lf"`` draws from the luminosity function
    ∝ 10^(LF_SLOPE·I)
    (D-068: every injection then has the same weight, n_eff = n), ``"uniform"`` as D-062."""
    lo, hi = INJ_IS
    u = rng.uniform(0.0, 1.0, size)
    if sampling == "uniform":
        return lo + (hi - lo) * u
    a = LF_SLOPE * math.log(10.0)
    return np.log(np.exp(a * lo) + u * (np.exp(a * hi) - np.exp(a * lo))) / a


def sampling_density(mag, sampling: str) -> np.ndarray:
    """Probability density of ``sample_magnitudes`` at ``mag``."""
    mag = np.asarray(mag, float)
    lo, hi = INJ_IS
    if sampling == "uniform":
        return np.full(mag.shape, 1.0 / (hi - lo))
    a = LF_SLOPE * math.log(10.0)
    return a * np.exp(a * mag) / (math.exp(a * hi) - math.exp(a * lo))


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
        row["shape"] = bool(prescreen_pass(scan["z_min"], scan["s_min"], scan["z_min2"]))
        sh = (
            epoch_artefact(_POP[0], _POP[1], scan, eid)
            if row["shape"] and _POP is not None
            else {"p": 1.0}
        )
        row["prescreen"] = bool(row["shape"] and sh["p"] >= COINC_P)
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
    procs: int,
    per_cell: int,
    per_ctrl: int,
    seed: int = 60,
    prescreen_only: bool = False,
    sampling: str = "lf",
    n_carriers: int = 300,
) -> Path:
    """Injection-recovery through the whole chain: Cut-0 emulation, pre-screen, fit, vetting
    (binary lens and the VSX/Gaia match excepted). W3: n = 1, ε < 0, u0 ~ U[0, 2); PSPL controls:
    u0 ~ U[0, 1) (they measure how often an ordinary event ends as a W3 survivor)."""
    field = moa.MoaField(FIELD)
    pre = read_prescreen()
    set_field_context(pre)  # the real field's coincidences and baseline calibration
    carriers = quiet_carriers(pre, n_carriers, seed)
    arrays = load_arrays(carriers, aux=True, pre=pre)
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
                        "Is": float(sample_magnitudes(rng, 1, sampling)[0]),
                        "ra": pos[eid][0],
                        "dec": pos[eid][1],
                    }
                    if prescreen_only:
                        prm["prescreen_only"] = True
                    jobs.append((kind, eid, t, f, sf, ax, prm))
    print(f"{len(jobs)} injections on {len(arrays)} carriers", flush=True)
    t1 = time.time()
    rows = []
    with _pool(procs) as pool:
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
        is_sampling=sampling,
        baseline_chi2=baseline_threshold(),
        n_carriers=len(arrays),
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


def field_star_counts(field: str | None = None) -> dict:
    """Monitored stars (10 ≤ I_s ≤ 21.4) of a field: Nunota et al. 2024 Table 1 where published
    (their N_s counts only the subfields they used; the screen covers all 80, so the published value
    is adopted as is — conservative — and N_s × 80 / n_sub is the upper end), else the N_s per
    Cut-0 object model (median, min, max; ASSUMPTION)."""
    f = field_number(field)
    if f in moa.NUNOTA_NS:
        nsub, ns = moa.NUNOTA_NS[f]
        return {"n_s": float(ns), "n_s_low": float(ns), "n_s_high": ns * 80.0 / nsub,
                "n_s_source": f"Nunota et al. 2024 Table 1 ({nsub}/80 subfields)"}  # fmt: skip
    est, lo, hi = moa.star_count_estimate(moa.CUT0_PER_FIELD[f])
    return {
        "n_s": est,
        "n_s_low": lo,
        "n_s_high": hi,
        "n_s_source": "N_s per Cut-0 object model (median of Nunota et al.'s 20 fields)",
    }


def run_limit(efficiency_only: bool = False) -> Path:
    """95 % limit per monitored star per year; needs a complete null vetting (zero survivors).
    ``efficiency_only``: a field with flags still open gets its efficiency table only
    (``efficiency_<field>.ecsv``, no rate columns; never read by ``combine``)."""
    vet = json.loads((out_dir() / f"vetting_{FIELD}.json").read_text())
    open_flags = [o["event_id"] for o in vet["flags"] if o.get("survives")]
    if open_flags and not efficiency_only:
        raise SystemExit(f"no zero-event limit: flags survive: {open_flags}")
    if vet.get("fit_errors", ["vetting record predates fit_errors"]):
        raise SystemExit(f"no zero-event limit: passes without a fit: {vet.get('fit_errors')}")
    inj = Table.read(out_dir() / f"injections_{FIELD}.ecsv")
    if not no_error(inj).all():  # dropping them would bias the efficiency upward
        raise SystemExit(f"{int((~no_error(inj)).sum())} injections failed; rerun inject")
    sampling = inj.meta.get("is_sampling", "uniform")
    ns = field_star_counts()
    n_s, n_lo = ns["n_s"], ns["n_s_low"]
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
            wt = lf_weights(w["Is"]) / sampling_density(w["Is"], sampling)
            eff_inj = float(np.sum(wt * rec) / np.sum(wt)) if len(w) else np.nan
            eff = eff_inj * frac  # stars brighter than 14.2 counted with efficiency 0
            lim = 3.0 / (n_s * years * eff) if eff > 0 else np.inf
            bright = np.asarray(w["Is"]) < 19.0
            rows.append(
                {
                    "field": FIELD,
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
                    "n_recovered": int(rec.sum()),
                    "p_recovered_bright": float(rec[bright].mean()) if bright.any() else np.nan,
                    "eff_per_star": eff,
                    # Kish effective sample size of the LF weights
                    "n_eff_lf": float(wt.sum() ** 2 / np.sum(wt**2)) if len(w) else 0.0,
                    "n_s": n_s,
                    "n_s_low": n_lo,
                    "years": years,
                    "rate95_per_star_yr": lim,
                    "rate95_conservative": 3.0 / (n_lo * years * eff) if eff > 0 else np.inf,
                    "n_ctrl": len(ctrl),
                    "ctrl_false_w3": int(np.sum(np.asarray(ctrl["recovered"], bool))),
                }
            )
    tab = Table(rows)
    tab.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source=f"scripts/w3_moa.py limit: injection-recovery through the full {FIELD} chain",
        n_s_source=ns["n_s_source"],
        baseline_chi2=inj.meta.get("baseline_chi2"),
        is_sampling=sampling,
        assumptions=(
            f"N_s = {n_s:.3g} (range {n_lo:.3g}–{ns['n_s_high']:.3g}) monitored stars, "
            f"10 ≤ I ≤ 21.4 ({ns['n_s_source']}); T = {years:.2f} yr; LF ∝ 10^({LF_SLOPE} I) "
            f"(Gaia DR3 RP in gb22, assumed for every field); injections I_s in 14.2–21.4 "
            f"({sampling} sampling), MOA-Red ≈ I; carrier noise unchanged by the source "
            "(sky/blend dominated); Cut-0 emulated at light-curve level; 95 % Poisson for 0 "
            "events = 3.0; mass: n = 1, D_L = 4 kpc, D_S = 8 kpc, μ_rel = 5 mas/yr "
            "(model_prediction)"
        ),
    )
    if efficiency_only:
        tab.remove_columns(["rate95_per_star_yr", "rate95_conservative"])
        tab.meta["open_flags"] = open_flags
        tab.meta["source"] += " (efficiency only: flags open, no limit)"
        path = out_dir() / f"efficiency_{FIELD}.ecsv"
        tab.write(path, overwrite=True)
        tab.write(results_dir() / f"efficiency_{FIELD}.ecsv", overwrite=True)
        tab.pprint(max_width=250, max_lines=50)
        return path
    path = out_dir() / f"limits_{FIELD}.ecsv"
    tab.write(path, overwrite=True)
    tab.write(results_dir() / f"limits_{FIELD}.ecsv", overwrite=True)  # tracked (small)
    tab.pprint(max_width=250, max_lines=50)
    return path


def run_combine() -> Path:
    """Combined 95 % limit over every field with a tracked zero-survivor limit table:
    Γ₉₅ = 3 / Σ_f N_s,f T ε_f per t_E × ρ cell (``derived``)."""
    files = sorted(results_dir().glob("limits_gb*.ecsv"))
    if not files:
        raise SystemExit("no per-field limit tables in results/w3_moa/")
    tabs = [Table.read(p) for p in files]
    rows = []
    for te in INJ_TE:
        for rho in INJ_RHO:
            exp = exp_lo = 0.0
            n_inj = n_rec = 0
            for t in tabs:
                r = t[(np.asarray(t["tE_days"]) == te) & (np.asarray(t["rho"]) == rho)][0]
                exp += float(r["n_s"]) * float(r["years"]) * float(r["eff_per_star"])
                exp_lo += float(r["n_s_low"]) * float(r["years"]) * float(r["eff_per_star"])
                n_inj += int(r["n_inj"])
                n_rec += int(r["n_recovered"])
            rows.append({
                "tE_days": te, "rho": rho,
                "mass_msun_model": (te / w3.einstein_time_days(1.0)[0]) ** 2,
                "n_fields": len(tabs), "n_inj": n_inj, "n_recovered": n_rec,
                "star_years_eff": exp,
                "rate95_per_star_yr": 3.0 / exp if exp > 0 else np.inf,
                "rate95_conservative": 3.0 / exp_lo if exp_lo > 0 else np.inf,
            })  # fmt: skip
    out = Table(rows)
    out.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source="scripts/w3_moa.py combine: " + ", ".join(str(t["field"][0]) for t in tabs),
        fields=[str(t["field"][0]) for t in tabs],
        assumptions="per-field assumptions in each limits_<field>.ecsv; zero survivors in all",
    )
    path = results_dir() / "limits_combined.ecsv"
    out.write(path, overwrite=True)
    out.pprint(max_width=250)
    return path


def write_vetting_summary() -> Path:
    """Tracked compact vetting record of a field (no fit parameters): funnel, flags, tests."""
    vet = json.loads((out_dir() / f"vetting_{FIELD}.json").read_text())
    pre = read_prescreen()
    rec = {
        "provenance": "derived",
        "field": FIELD,
        "n_light_curves": n_light_curves(pre),
        "n_shape": int(pre.meta.get("n_shape", 0)),
        "n_passes": len(passes(pre)),
        "baseline_chi2": baseline_threshold(),
        "baseline_calibration": calibrate_baseline(pre) if "chi2_const" in pre.colnames else None,
        "n_flags": vet["n_flags"],
        "funnel": vetting_funnel(vet),
        "flags": [
            {
                "event_id": o["event_id"],
                "exotic": o.get("exotic"),
                "dbic_all": o.get("dbic_all"),
                "survives": o.get("survives"),
                "tests": [[n, bool(ok), note] for n, ok, note in o["tests"]],
            }
            for o in vet["flags"]
        ],
    }
    path = results_dir() / f"vetting_{FIELD}.json"
    path.write_text(json.dumps(rec, indent=1, default=float, ensure_ascii=False) + "\n")
    return path


LAST_MODIFIED = {  # HTTP Last-Modified of the pinned files, read 2026-10-08 (the release version)
    "metadata.ipac.tar.gz": "2023-10-13",
    "gb22.tar": "2023-10-10",
}


def field_range_digest(field: str) -> str | None:
    """Content pin of a streamed field tar from its tracked pre-screen chunks: sha256 over the
    ordered sha256s of its 64 MiB ranges (``moa_stream.range_digest``); None until complete."""
    n = n_chunks(field)
    shas = []
    for k in range(n):
        p = prescreen_chunk_path(field, k, n)
        if not p.exists():
            return None
        meta = Table.read(p, format="ascii.ecsv").meta
        if meta.get("segment_bytes") != SEG_BYTES or "segment_sha256" not in meta:
            return None
        shas.extend(meta["segment_sha256"])
    return moa_stream.range_digest(shas)


def write_manifest() -> Path:
    """Pinned MOA-II files as a tracked manifest (URL, sha256, size, retrieval date, release).
    Downloaded files carry a whole-file sha256; streamed field tars carry ``range_sha256`` (sha256
    over the sha256s of consecutive 64 MiB byte ranges, recorded per chunk in results/w3_moa/)."""
    rows = []
    release = "MOA-II 9-year bulge release (2006-2014), NASA Exoplanet Archive"
    for url, (sha, size) in moa.FILES.items():
        name = url.rsplit("/", 1)[-1]
        rows.append({"url": url, "sha256": sha, "range_sha256": "", "size_bytes": size,
                     "retrieved_utc": "2026-10-08", "last_modified": LAST_MODIFIED[name],
                     "release": release})  # fmt: skip
    for f in sorted(moa.TAR_BYTES):
        dig = field_range_digest(f"gb{f}")
        if dig:
            rows.append({"url": moa.tar_url(f), "sha256": "", "range_sha256": dig,
                         "size_bytes": moa.TAR_BYTES[f], "retrieved_utc": "2026-10-08",
                         "last_modified": moa.TAR_LAST_MODIFIED[f],
                         "release": release})  # fmt: skip
    tab = Table(rows)
    tab.meta.update(
        provenance=schema.Provenance.OBSERVED.value,
        source="jwst_anomaly.moa.FILES and streamed field tars (MOA-II 9-year; D-062, D-068)",
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


def _init_worker(ctx: dict | None = None) -> None:
    """Pool initializer: spawned workers (Windows, macOS) re-import this module without running
    `main`, so the MOA fit bounds, the field, its shared-epoch populations and its baseline
    calibration are set here, not inherited."""
    global _POP, _BASELINE, FIELD
    use_moa_fit_params()
    ctx = ctx or {}
    _POP = ctx.get("pop")
    _BASELINE = ctx.get("baseline")
    FIELD = ctx.get("field", FIELD)


def _ctx() -> dict:
    return {"pop": _POP, "baseline": _BASELINE, "field": FIELD}


def _pool(procs: int):
    return Pool(procs, initializer=_init_worker, initargs=(_ctx(),))


FITS_PER_CHUNK = 120  # pre-screen passes per tracked fit chunk (~30 min on 4 cores)


def run_field(procs: int, conns: int, per_cell: int, per_ctrl: int) -> None:
    """Every stage for ``FIELD``, resumable: streamed pre-screen (finished chunks skipped), fits in
    tracked chunks (finished chunks skipped), vetting, contact sheet, injections, limit. Stops
    before the limit if any flag survives (``run_limit`` refuses)."""
    t0 = time.time()
    run_stream_prescreen(FIELD, procs, conns)
    merge_prescreen(FIELD)
    pre = read_prescreen()
    n_pass = len(passes(pre))
    n = max(1, math.ceil(n_pass / FITS_PER_CHUNK))
    ids = list(passes(pre)["event_id"])
    for k in range(n):
        path = results_dir() / f"{chunk_name((k, n))}.gz"
        ok = (
            path.exists()
            and fit_chunk_problem(Table.read(path, format="ascii.ecsv"), ids, k, n) is None
        )
        if not ok:
            run_fit(procs, None, (k, n))
    merge_chunks(n)
    run_vet(procs)
    contact_sheet(out_dir() / f"contact_sheet_{FIELD}.png")
    vet = json.loads((out_dir() / f"vetting_{FIELD}.json").read_text())
    alive = [o["event_id"] for o in vet["flags"] if o.get("survives")]
    set_field_context(pre)
    write_vetting_summary()
    if alive:
        raise SystemExit(f"{FIELD}: flags survive every test: {alive} — stop and report")
    run_inject(procs, per_cell, per_ctrl)
    run_limit()
    print(f"{FIELD}: chain done in {time.time() - t0:.0f} s", flush=True)


def parse_chunk_list(text: str | None) -> list[int] | None:
    """``"3"``, ``"1-5"``, ``"1,4,7-9"`` (1-based) -> 0-based chunk numbers; None = all."""
    if not text:
        return None
    out = []
    for part in text.split(","):
        a, _, b = part.partition("-")
        out.extend(range(int(a) - 1, int(b or a)))
    return out


def set_field(field: str) -> None:
    global FIELD
    if field_number(field) not in moa.CUT0_PER_FIELD:
        raise SystemExit(f"unknown MOA field {field!r}")
    FIELD = field


def main(argv=None) -> int:
    use_moa_fit_params()
    ncpu = os.cpu_count() or 4
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--field", default="gb22", help="MOA field, gb1 … gb22")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prescreen", help="stream the field tar and pre-screen every light curve")
    p.add_argument("--procs", type=int, default=ncpu)
    p.add_argument("--conns", type=int, default=8, help="concurrent HTTP range reads")
    p.add_argument("--source", choices=("auto", "local", "http"), default="auto")
    p.add_argument("--chunks", default=None, help="1-based chunk list, e.g. 1-5,9 (default all)")
    sub.add_parser("merge-prescreen", help="join the tracked pre-screen chunks of the field")
    f = sub.add_parser("fit", help="ordinary and exotic fits of the pre-screen passes")
    f.add_argument("--procs", type=int, default=ncpu)
    f.add_argument("--limit", type=int, default=None)
    f.add_argument("--chunk", default=None, help="K/N: fit every N-th pass from the K-th")
    m = sub.add_parser("merge-chunks", help="join the tracked chunk fit tables 1..N")
    m.add_argument("--n", type=int, required=True)
    v = sub.add_parser("vet", help="vet the flags of the fit stage")
    v.add_argument("--procs", type=int, default=ncpu)
    c = sub.add_parser("sheet", help="contact sheet of the vetted flags")
    c.add_argument("--out", type=Path, default=None)
    c.add_argument("--survivors-only", action="store_true")
    i = sub.add_parser("inject", help="injection-recovery through the whole chain")
    i.add_argument("--procs", type=int, default=ncpu)
    i.add_argument("--per-cell", type=int, default=200)
    i.add_argument("--per-ctrl", type=int, default=40)
    i.add_argument("--carriers", type=int, default=300)
    i.add_argument("--sampling", choices=("lf", "uniform"), default="lf")
    i.add_argument("--seed", type=int, default=60)
    i.add_argument("--prescreen-only", action="store_true", help="Cut-0 and pre-screen only")
    lim = sub.add_parser("limit", help="95 %% rate limit per monitored star per year")
    lim.add_argument("--efficiency-only", action="store_true", help="flags open: ε table only")
    sub.add_parser("combine", help="combined limit over the fields with tracked limit tables")
    sub.add_parser("summary", help="tracked compact vetting record of the field")
    sub.add_parser("manifest", help="write data/manifests/moa_ii.ecsv")
    r = sub.add_parser("run-field", help="every stage for the field (resumable)")
    r.add_argument("--procs", type=int, default=ncpu)
    r.add_argument("--conns", type=int, default=8)
    r.add_argument("--per-cell", type=int, default=200)
    r.add_argument("--per-ctrl", type=int, default=40)
    a = ap.parse_args(argv)
    if a.cmd in ("fit", "vet", "inject", "run-field") and not w3.have_mm():
        # Without MulensModel the parallax refit and the binary-lens test are skipped, which
        # inflates the injection efficiency and leaves real flags incomplete (D-068).
        ap.error(f"{a.cmd} needs MulensModel: pip install -e '.[mulens]'")
    set_field(a.field)
    if a.cmd == "prescreen":
        run_stream_prescreen(
            FIELD, min(a.procs, ncpu), a.conns, parse_chunk_list(a.chunks), a.source
        )
        if all(_chunk_done(FIELD, k, n_chunks()) for k in range(n_chunks())):
            merge_prescreen(FIELD)
    elif a.cmd == "merge-prescreen":
        merge_prescreen(FIELD)
    elif a.cmd == "fit":
        run_fit(min(a.procs, ncpu), a.limit, w3.parse_chunk(a.chunk))
    elif a.cmd == "merge-chunks":
        merge_chunks(a.n)
    elif a.cmd == "vet":
        run_vet(min(a.procs, ncpu))
    elif a.cmd == "inject":
        run_inject(min(a.procs, ncpu), a.per_cell, a.per_ctrl, a.seed, a.prescreen_only,
                   a.sampling, a.carriers)  # fmt: skip
    elif a.cmd == "limit":
        run_limit(a.efficiency_only)
    elif a.cmd == "combine":
        run_combine()
    elif a.cmd == "summary":
        set_field_context(read_prescreen())
        print(write_vetting_summary())
    elif a.cmd == "run-field":
        run_field(min(a.procs, ncpu), a.conns, a.per_cell, a.per_ctrl)
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
