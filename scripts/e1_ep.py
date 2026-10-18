"""E1 + Einstein Probe: an EP/WXT transient node table from the GCN Circulars archive (D-074 E1
node set; frozen plan ``results/e1_ep/PLAN.md``). No public EP event table exists (2026-10-18
search), so each event's first WXT report circular is parsed: onset time, J2000 position and the
90 % error radius.

The archive (~31 MB, rewritten daily) is streamed to ``$JWST_ANOMALY_DATA/e1_ep/``, pinned by
sha256 in ``data/manifests/e1_ep.ecsv`` and may be deleted after use; the parsed table
``results/e1_ep/ep_events.ecsv`` (with circular ids) is the record the counts use.

  python scripts/e1_ep.py build [--refresh]
  python scripts/e1_ep.py count [--jit 10000] [--jitday 10000] [--inject 40] [--cpu 4]

``count`` reuses the GBM and BAT inputs pinned by ``scripts/e1_swift.py``
(``data/manifests/e1_swift.ecsv``).
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import os
import re
import sys
import tarfile
import time
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import requests
from astropy.table import Table

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import e1_events as e1  # noqa: E402
import e1_swift as sw  # noqa: E402

from jwst_anomaly import event_network as en  # noqa: E402
from jwst_anomaly.paths import data_root  # noqa: E402
from jwst_anomaly.schema import Provenance  # noqa: E402

OUT = ROOT / "results" / "e1_ep"
MANIFEST = ROOT / "data" / "manifests" / "e1_ep.ecsv"
ARCHIVE_URL = "https://gcn.nasa.gov/circulars/archive.json.tar.gz"
ARCHIVE = "gcn_circulars_archive.json.tar.gz"
#: 2-D Gaussian: r90 = 2.146 sigma.
R90_TO_SIGMA = 2.146
#: Tested channels (frozen plan), then the positive-control channels (never tested).
TESTED = (("EP", "EP"), ("EP", "BAT"), ("EP", "GBM"))
CHANNELS = TESTED
CONTROL = (("EP", "BAT"), ("EP", "GBM"))
#: ASSUMPTION: EP orbital period (s), ~96 min at ~600 km.
EP_ORBIT = 5760.0
P = en.Params(quantum_s={"GBM": en.FERMI_ORBIT, "BAT": sw.SWIFT_ORBIT, "EP": EP_ORBIT})
P_DAY = en.Params(quantum_s={"GBM": en.FERMI_ORBIT, "BAT": sw.SWIFT_ORBIT, "EP": EP_ORBIT})
#: ASSUMPTIONs of the frozen control rule.
CONTROL_MIN_PAIRS = 5
CONTROL_MAX_P = 1.8e-5
#: ASSUMPTION: family-wise 3 sigma (one-sided), as D-074.
DETECT_P = e1.DETECT_P
#: ASSUMPTION: an onset more than this many days outside the date in the EP name is a parse error.
MAX_NAME_DAYS = 1.0

NAME = re.compile(r"\bEP(\d{6})[a-z]{1,2}\b")
POS = re.compile(
    r"R\.?\s*A\.?\s*(?:\(J2000\))?\s*[=:]\s*([0-9.]+)\s*(?:deg)?[^0-9]{0,30}?"
    r"DEC\.?\s*(?:\(J2000\))?\s*[=:]\s*([-+\u2212]?)\s*([0-9.]+)",
    re.I | re.S,
)
ERR = re.compile(
    r"uncertainty(?:\s+radius)?\s+of\s*(?:about\s*|~\s*)?([0-9.]+)\s*(arcmin|arcsec)", re.I
)
#: a counterpart retraction (e.g. EP250903a, GCN 41696) does not retract the X-ray trigger
RETRACT = re.compile(r"not a real (?:source|astrophysical)", re.I)
TS = r"(20\d\d-\d\d-\d\d)[T ](\d\d:\d\d:\d\d(?:\.\d+)?)"
ONSET = re.compile(r"(?:began|begins|started|starting|onset)[^.]{0,40}?at\s*" + TS, re.I)
ANY_TS = re.compile(TS)


def _mjd(date: str, clock: str) -> float:
    d = datetime.fromisoformat(f"{date}T{clock}").replace(tzinfo=UTC)
    return d.timestamp() / 86400.0 + 40587.0


def parse_circular(subject: str, body: str) -> dict | None:
    """One WXT report -> {name, mjd, ra, dec, r90_deg, time_rule} or None.

    Frozen rules (PLAN.md): the EP name comes from the subject; the body must mention WXT before
    the position; time is the light-curve onset ("began/started ... at <UTC>") if stated, else the
    first UTC timestamp; the 90 % radius is the first "uncertainty [radius] of <x> arcmin" (else the
    first arcsec one, flagged in ``r90_unit``); the onset
    must lie within MAX_NAME_DAYS of the date in the name."""
    n = NAME.search(subject)
    pos = POS.search(body)
    errs = list(ERR.finditer(body))
    # the WXT radius is in arcmin; an arcsec radius (FXT or refined) is used only if none is
    err = next((e for e in errs if e.group(2).lower() == "arcmin"), errs[0] if errs else None)
    if not (n and pos and err) or "WXT" not in body[: pos.start() + 400]:
        return None
    on = ONSET.search(body)
    ts = on.groups() if on else (ANY_TS.search(body).groups() if ANY_TS.search(body) else None)
    if ts is None:
        return None
    mjd = _mjd(*ts)
    yymmdd = n.group(1)
    name_mjd = _mjd(f"20{yymmdd[:2]}-{yymmdd[2:4]}-{yymmdd[4:]}", "00:00:00")
    if not (-MAX_NAME_DAYS <= mjd - name_mjd <= 1.0 + MAX_NAME_DAYS):
        return None
    sign = -1.0 if pos.group(2) in ("-", "\u2212") else 1.0
    ra, dec = float(pos.group(1)), sign * float(pos.group(3))
    r90 = float(err.group(1)) / (60.0 if err.group(2).lower() == "arcmin" else 3600.0)
    if not (0 <= ra < 360 and -90 <= dec <= 90 and r90 > 0):
        return None
    return {
        "name": n.group(0),
        "mjd": mjd,
        "ra": ra,
        "dec": dec,
        "r90_deg": r90,
        "time_rule": "onset" if on else "first",
        "r90_unit": err.group(2).lower(),
    }


def ep_events(circulars) -> Table:
    """Circular dicts -> one row per EP name, from the lowest-numbered circular that parses; names
    with any circular calling the trigger "not a real source" are dropped."""
    circulars = list(circulars)
    retracted = {
        NAME.search(c["subject"]).group(0)
        for c in circulars
        if NAME.search(c.get("subject", ""))
        and RETRACT.search(c.get("subject", "") + c.get("body", ""))
    }
    best: dict[str, tuple[int, dict]] = {}
    for c in circulars:
        if not NAME.search(c.get("subject", "")):
            continue
        r = parse_circular(c["subject"], c.get("body", ""))
        if r is None:
            continue
        cid = int(c["circularId"])
        if r["name"] in retracted:
            continue
        if r["name"] not in best or cid < best[r["name"]][0]:
            best[r["name"]] = (cid, r)
    rows = [dict(r, circular=cid) for cid, r in (best[k] for k in sorted(best))]
    t = Table(
        rows=[
            (r["name"], r["circular"], r["mjd"], r["ra"], r["dec"], r["r90_deg"],
             r["r90_deg"] / R90_TO_SIGMA, r["time_rule"], r["r90_unit"])
            for r in rows
        ],
        names=("name", "circular", "mjd", "ra", "dec", "r90_deg", "sigma", "time_rule", "r90_unit"),
        dtype=(str, int, float, float, float, float, float, str, str),
    )  # fmt: skip
    t.meta = {
        "provenance": str(Provenance.DERIVED),
        "source": "GCN Circulars archive (EP/WXT reports), parsed by scripts/e1_ep.py",
    }
    return t


def iter_archive(path: Path):
    with tarfile.open(path) as tf:
        for m in tf:
            if m.isfile() and m.name.endswith(".json"):
                yield json.load(tf.extractfile(m))


def fetch(refresh: bool) -> Path:
    d = data_root() / "e1_ep"
    d.mkdir(parents=True, exist_ok=True)
    path = d / ARCHIVE
    if not path.exists() or refresh:
        r = requests.get(ARCHIVE_URL, timeout=600)
        r.raise_for_status()
        path.write_bytes(r.content)
    import hashlib

    h = hashlib.sha256(path.read_bytes()).hexdigest()
    if MANIFEST.exists() and not refresh:
        pin = Table.read(MANIFEST, format="ascii.ecsv")["sha256"][0]
        if pin != h:
            raise SystemExit(f"{ARCHIVE}: sha256 {h} != pinned {pin} (rerun with --refresh)")
    else:
        t = Table(rows=[(ARCHIVE_URL, ARCHIVE, h, path.stat().st_size)],
                  names=("uri", "file", "sha256", "size"))  # fmt: skip
        t.meta = {
            "provenance": str(Provenance.OBSERVED),
            "source": "GCN Circulars JSON archive (rewritten daily; the parsed table in "
            "results/e1_ep/ep_events.ecsv is the record)",
            "retrieved": time.strftime("%Y-%m-%d", time.gmtime()),
            "pipeline_version": "0.0.1",
            "reason": "~31 MB; text circulars only; delete after use in cloud sessions",
        }
        t.write(MANIFEST, format="ascii.ecsv", overwrite=True)
    return path


def load_ep(path: Path = OUT / "ep_events.ecsv") -> Table:
    t = Table.read(path, format="ascii.ecsv")
    return en._events(
        t["name"], "EP", t["mjd"], t["ra"], t["dec"], t["sigma"], "results/e1_ep/ep_events.ecsv"
    )


# ------------------------------------------------------------------------------------- counting
_W: dict = {}


def count_all(samples: dict[str, en.Sample], p: en.Params) -> np.ndarray:
    """Counts, shape (n_channels, n_windows, 3)."""
    return np.stack(
        [en.count_channel(samples[a], samples[b], p, same=(a == b)) for a, b in CHANNELS]
    )


def _init(samples):
    os.environ["OMP_NUM_THREADS"] = "1"
    _W["s"] = {k: samples[k] for k in {c for ch in CHANNELS for c in ch}}


def _null_chunk(args):
    kind, seeds = args
    s = _W["s"]
    q = P_DAY if kind == "jitday" else P
    out = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        out.append(count_all({k: en.scramble_jit(v, rng, q) for k, v in s.items()}, P))
    return np.stack(out)


def run_null(kind: str, n: int, samples, cpu: int, base_seed: int) -> np.ndarray:
    seeds = np.arange(n) + base_seed
    chunks = [(kind, seeds[i : i + 250]) for i in range(0, n, 250)]
    t0 = time.time()
    if cpu <= 1:
        _init(samples)
        res = [_null_chunk(c) for c in chunks]
    else:
        with cf.ProcessPoolExecutor(cpu, initializer=_init, initargs=(samples,)) as ex:
            res = list(ex.map(_null_chunk, chunks))
    print(f"null {kind}: {n} scrambles in {time.time() - t0:.0f} s", flush=True)
    return np.concatenate(res)


def test_mask(p: en.Params) -> np.ndarray:
    """Tested family: EP-EP same/wide at all five lag bins; EP-BAT and EP-GBM same/wide at the
    two longest lag bins (1 h-1 d, 1-7 d) only."""
    nwin, nlag = len(en.windows_for(p)), len(p.lag_edges) - 1
    m = np.zeros((len(CHANNELS), nwin, 3), bool)
    for c, ch in enumerate(CHANNELS):
        lo = 0 if ch == ("EP", "EP") else nlag - 2
        m[c, lo:nlag, :2] = True
    return m


def control_sum(counts: np.ndarray) -> np.ndarray:
    """Same-direction EP-BAT + EP-GBM pairs at 0-10 s, 10-100 s and 100 s-1 h, from counts of
    shape (..., n_channels, n_windows, 3): one count array or a null ensemble."""
    idx = [CHANNELS.index(c) for c in CONTROL]
    same = counts[..., :3, 0]  # (..., n_channels, 3); plain indexing keeps the axis order
    return np.take(same, idx, axis=-2).sum(axis=(-1, -2))


def control_passes(obs_pairs: int, null_sums: np.ndarray) -> tuple[bool, float]:
    """Frozen control rule: >= CONTROL_MIN_PAIRS and the summed count's analytic p (jit null
    ensemble of the same sum) <= CONTROL_MAX_P."""
    pa = float(en.analytic_p(np.array([obs_pairs]), np.asarray(null_sums, float)[:, None])[0])
    return (obs_pairs >= CONTROL_MIN_PAIRS) and (pa <= CONTROL_MAX_P), pa


def outcome(cell_p: np.ndarray, n_cells: int) -> str:
    """Frozen outcome: 'excess' if any tested cell's analytic p x n_cells <= DETECT_P."""
    return "excess" if float(np.min(cell_p)) * n_cells <= DETECT_P else "null"


def count(a) -> int:
    p = P
    ev = sw.load(sw.fetch(False))
    ev = {"GBM": ev["GBM"], "BAT": ev["BAT"], "EP": load_ep()}
    samples = {k: en.Sample.from_table(v) for k, v in ev.items()}
    for k, v in ev.items():
        print(k, len(v), f"MJD {v['mjd'].min():.1f}-{v['mjd'].max():.1f}", flush=True)
    obs = count_all(samples, p)
    jit = run_null("jit", a.jit, samples, a.cpu, 4_000_000)
    jday = run_null("jitday", a.jitday, samples, a.cpu, 5_000_000)
    mask = test_mask(p)
    n_cells = int(mask.sum())
    c_obs, c_null = int(control_sum(obs)), control_sum(jit)
    c_ok, c_p = control_passes(c_obs, c_null)

    wins = [w for w, _ in en.windows_for(p)]
    labels = en.lag_labels(p) + wins[len(p.lag_edges) - 1 :]
    z_jit, pa_jit = en.z_score(obs, jit), en.analytic_p(obs, jit)
    p_jit, p_jday, z_jday = (
        en.empirical_p(obs, jit),
        en.empirical_p(obs, jday),
        en.z_score(obs, jday),
    )
    ul = en.upper_limit_95(obs, jit)
    rows = []
    for c, (ca, cb) in enumerate(CHANNELS):
        for k, lab in enumerate(labels):
            for j, cl in enumerate(e1.CLASSES):
                rows.append(
                    (
                        f"{ca}-{cb}", lab, cl, bool(mask[c, k, j]), int(obs[c, k, j]),
                        round(float(jit[:, c, k, j].mean()), 2),
                        round(float(jit[:, c, k, j].std()), 2),
                        round(float(z_jit[c, k, j]), 2), float(p_jit[c, k, j]),
                        float(pa_jit[c, k, j]),
                        round(float(jday[:, c, k, j].mean()), 2),
                        round(float(z_jday[c, k, j]), 2), float(p_jday[c, k, j]),
                        round(float(ul[c, k, j]), 2),
                    )
                )  # fmt: skip
    counts = Table(
        rows=rows,
        names=(
            "channel", "lag", "cls", "tested", "obs", "jit_mean", "jit_sd", "z_jit", "p_jit",
            "p_jit_analytic", "jday_mean", "z_jday", "p_jday", "ul95_pairs",
        ),
    )  # fmt: skip
    counts.meta = {
        "provenance": str(Provenance.DERIVED),
        "source": "scripts/e1_ep.py count on results/e1_ep/ep_events.ecsv and "
        f"data/manifests/e1_swift.ecsv (GBM, BAT); jit null {a.jit}, jitday null {a.jitday}",
    }
    counts.write(OUT / "counts.ecsv", format="ascii.ecsv", overwrite=True)

    pa_t = pa_jit[mask]
    cells = [(f"{x}-{y}", lab, cl) for x, y in CHANNELS for lab in labels for cl in e1.CLASSES]
    tested_cells = [c for c, m in zip(cells, mask.ravel(), strict=True) if m]
    jmin = int(np.argmin(pa_t))
    res = outcome(pa_t, n_cells)

    dup = {}
    for ca, cb in CHANNELS:
        t = e1.close_pairs(ev[ca], ev[cb], p, ca == cb, p.lag_edges[-1])
        t.write(OUT / f"same_dir_pairs_{ca}-{cb}.ecsv", format="ascii.ecsv", overwrite=True)
        dup[f"{ca}-{cb}"] = len(t)

    # injections and limits for EP-EP only (frozen plan), through e1_events' code (in-process)
    # e1.injections sets its per-cell threshold to DETECT_P / mask.sum(); scale DETECT_P so the
    # threshold is DETECT_P / n_cells of the whole family (18), not of the EP-EP subset
    e1.CHANNELS = (("EP", "EP"),)
    e1.DETECT_P = DETECT_P * int(mask[:1].sum()) / n_cells
    inj, p_star = e1.injections(
        samples, p, jit[:, :1], mask[:1], a.inject, np.random.default_rng(13)
    )
    e1.DETECT_P = DETECT_P
    inj.write(OUT / "injections.ecsv", format="ascii.ecsv", overwrite=True)
    lim = e1.limits(counts, inj, ev, p)
    lim.write(OUT / "limits.ecsv", format="ascii.ecsv", overwrite=True)

    summary = {
        "provenance": "derived",
        "plan": "results/e1_ep/PLAN.md",
        "events": {k: len(v) for k, v in ev.items()},
        "ep_mjd_range": [float(ev["EP"]["mjd"].min()), float(ev["EP"]["mjd"].max())],
        "n_jit": a.jit,
        "n_jitday": a.jitday,
        "tested_cells": n_cells,
        "positive_control_ep_grb_same_0_1h": {
            "obs": c_obs,
            "jit_mean": float(c_null.mean()),
            "jit_sd": float(c_null.std()),
            "analytic_p": c_p,
            "passes": bool(c_ok),
            "rule": f">= {CONTROL_MIN_PAIRS} pairs and analytic p <= {CONTROL_MAX_P}",
        },
        "outcome": res if c_ok else "no result (control failed)",
        "min_cell": {
            "cell": " ".join(tested_cells[jmin]),
            "analytic_p": float(pa_t[jmin]),
            "bonferroni_p": float(min(1.0, pa_t[jmin] * n_cells)),
        },
        "global_jit_pooled": dict(
            zip(("min_cell_p", "trials_corrected_p"), en.global_p(obs, jit, mask), strict=True)
        ),
        "global_jitday_pooled": dict(
            zip(("min_cell_p", "trials_corrected_p"), en.global_p(obs, jday, mask), strict=True)
        ),
        "injection_detection_p_star": p_star,
        "same_direction_pairs_within_7d": dup,
        "ep_orbit_s": EP_ORBIT,
    }
    with open(OUT / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=1, default=float)
        fh.write("\n")
    print(json.dumps({k: summary[k] for k in ("events", "positive_control_ep_grb_same_0_1h",
                                              "outcome", "min_cell", "global_jit_pooled")},
                     indent=1))  # fmt: skip
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["build", "count"])
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("--jit", type=int, default=10000)
    ap.add_argument("--jitday", type=int, default=10000)
    ap.add_argument("--inject", type=int, default=40, help="trials per (bin, n)")
    ap.add_argument("--cpu", type=int, default=1)
    a = ap.parse_args(argv)
    if a.cmd == "count":
        return count(a)
    path = fetch(a.refresh)
    circ = list(iter_archive(path))
    ep_names = {NAME.search(c.get("subject", "")).group(0) for c in circ
                if NAME.search(c.get("subject", "")) and "WXT" in c.get("body", "")}  # fmt: skip
    ev = ep_events(circ)
    OUT.mkdir(parents=True, exist_ok=True)
    ev.write(OUT / "ep_events.ecsv", format="ascii.ecsv", overwrite=True)
    summary = {
        "provenance": "derived",
        "circulars": len(circ),
        "ep_names_with_wxt_mention": len(ep_names),
        "parsed_events": len(ev),
        "unparsed_names": sorted(ep_names - set(ev["name"])),
        "retracted_or_unparsed": "see unparsed_names",
        "r90_unit_counts": {k: int(np.sum(ev["r90_unit"] == k)) for k in ("arcmin", "arcsec")},
        "time_rule_counts": {k: int(np.sum(ev["time_rule"] == k)) for k in ("onset", "first")},
        "mjd_range": [float(ev["mjd"].min()), float(ev["mjd"].max())],
        "r90_arcmin_quantiles_5_50_95": [
            round(float(q) * 60, 3) for q in np.quantile(ev["r90_deg"], [0.05, 0.5, 0.95])
        ],
    }
    with open(OUT / "build_summary.json", "w") as fh:
        json.dump(summary, fh, indent=1)
        fh.write("\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "unparsed_names"}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
