"""Gaia DR3 microlensing candidates as a :class:`signatures.LightCurveSurvey` (D-061).

- Events: ``gaiadr3.vari_microlensing`` (Wyrzykowski et al. 2023, A&A 674, A23,
  doi:10.1051/0004-6361/202243756, arXiv:2206.06121): 363 candidates, 2014-07-25 … 2017-05-28.
  The table has no column for the paper's two samples; ``method`` (A, B or A+B) is parsed from
  the paper's Table D.1 in the pinned arXiv v2 source. Sample A is the automated selection
  (membership score and cuts, Appendix C); sample B is OGLE / ASAS-SN / Gaia-alert events matched
  to Gaia and inspected by eye.
- Light curves: Gaia DR3 epoch photometry from the archive's DataLink service, fetched in parallel
  batches once and cached per source under ``$JWST_ANOMALY_DATA/raw/gaia_dr3_mulens/``. Rows
  flagged ``variability_flag_g_reject`` or ``rejected_by_photometry`` are dropped (the DR3
  outlier filter). ``time`` is the archive's ``g_transit_time`` + 2455197.5 (BJD(TCB), full JD,
  so MulensModel parallax works; Gaia is at L2, ~0.01 au from Earth: ASSUMPTION, negligible).
  ``mag_err`` = 2.5 / ln 10 / ``g_transit_flux_over_error``, as published (unscaled).

Everything here is ``observed`` as published; the ``paczynski0_*`` columns are the published
Level 0 PSPL fits (no blending).
"""

from __future__ import annotations

import csv
import http.client
import io
import os
import re
import sys
import tarfile
import time
import urllib.parse
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.table import Table

from jwst_anomaly import paths, photometry, schema
from jwst_anomaly.signatures import standard_light_curve

TAP_SYNC = "https://gea.esac.esa.int/tap-server/tap/sync"
DATALINK = "https://gea.esac.esa.int/data-server/data"
EVENTS_QUERY = "SELECT * FROM gaiadr3.vari_microlensing"
GAIA_T0 = 2455197.5  # Gaia time = BJD(TCB) − 2455197.5
T_FIRST, T_LAST = 2456863.9375, 2457901.86389  # DR3 photometry window (JD; Wyrzykowski+23 §2)
PAPER_SRC = "https://arxiv.org/src/2206.06121v2"
PAPER_SHA = "40c60eee14f1a5e691cf7efb5b1233878fcabd37496eb4a0865fa5e7afa07295"
PAPER_TABLE = "table-events-uniq-sort-mag-j2000.tex"
REFERENCE = "Wyrzykowski et al. 2023, A&A 674, A23"
BATCH = 12  # source ids per DataLink request (~3 s per source on the server)
WORKERS = 8


def cache_dir() -> Path:
    d = paths.data_root() / "raw" / "gaia_dr3_mulens"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _write_atomic(path: Path, data: bytes) -> None:
    """Write via a temporary file so an interrupted run never leaves a truncated cache file."""
    tmp = path.with_name(path.name + f".{os.getpid()}.part")
    tmp.write_bytes(data)
    tmp.replace(path)


def _tap_csv(query: str) -> str:
    data = urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query}
    ).encode()
    with urllib.request.urlopen(TAP_SYNC, data=data, timeout=300) as resp:  # noqa: S310
        return resp.read().decode()


def parse_method_table(tex: str) -> dict[int, str]:
    """``source_id -> method`` (A, B or A+B) from the rows of the paper's Table D.1."""
    out = {}
    for line in tex.splitlines():
        m = re.match(r"\s*\d+\s*&\s*(\d{6,})\s*&.*&\s*(A\+B|A|B)\s*\\\\", line)
        if m:
            out[int(m.group(1))] = m.group(2)
    return out


def parse_epoch_csv(text: str) -> Table:
    """One DataLink epoch-photometry CSV -> G ``standard_light_curve``, rejected rows dropped."""
    rows = list(csv.DictReader(io.StringIO(text)))
    keep = [
        r
        for r in rows
        if r["g_transit_time"]
        and r["g_transit_mag"]
        and r["g_transit_flux_over_error"]
        and float(r["g_transit_flux_over_error"]) > 0
        and r["variability_flag_g_reject"] != "true"
        and r["rejected_by_photometry"] != "true"
    ]
    t = np.array([float(r["g_transit_time"]) for r in keep]) + GAIA_T0
    mag = np.array([float(r["g_transit_mag"]) for r in keep])
    foe = np.array([float(r["g_transit_flux_over_error"]) for r in keep])
    sid = rows[0]["source_id"] if rows else "?"
    lc = standard_light_curve(
        t,
        mag,
        2.5 / np.log(10.0) / foe,
        "G",
        source=f"Gaia DR3 epoch photometry (DataLink) source_id {sid}",
        time_system="BJD(TCB)",
    )
    lc.meta["n_rejected"] = len(rows) - len(keep)
    return lc


class GaiaDR3Microlensing:
    """Adapter for the Gaia DR3 microlensing candidates (``vari_microlensing``)."""

    name = "gaia-dr3-vari-microlensing"

    def __init__(self, cache: Path | None = None):
        self._cache = cache or cache_dir()
        self._events: Table | None = None

    # ------------------------------------------------------------------ files
    def methods(self) -> dict[int, str]:
        tar_path = photometry.fetch_catalog(PAPER_SRC, PAPER_SHA)
        with tarfile.open(tar_path, "r:*") as tar:
            tex = tar.extractfile(PAPER_TABLE).read().decode()
        return parse_method_table(tex)

    def events(self) -> Table:
        if self._events is None:
            path = self._cache / "vari_microlensing.csv"
            if not path.exists():
                _write_atomic(path, _tap_csv(EVENTS_QUERY).encode())
            t = Table.read(path, format="ascii.csv")
            meth = self.methods()
            t["event_id"] = [str(s) for s in t["source_id"]]
            t["method"] = [meth.get(int(s), "") for s in t["source_id"]]
            t["t0_pub"] = np.asarray(t["paczynski0_tmax"], float) + GAIA_T0
            t["tE_pub"] = np.asarray(t["paczynski0_te"], float)
            t["u0_pub"] = np.abs(np.asarray(t["paczynski0_u0"], float))
            t.meta.update(
                provenance=schema.Provenance.OBSERVED.value,
                source=f"Gaia archive TAP `{EVENTS_QUERY}` ({REFERENCE}); method from "
                f"{PAPER_SRC} Table D.1; *_pub = published Level 0 PSPL fit",
            )
            self._events = t
        return self._events

    def radec(self) -> dict[str, tuple[float, float]]:
        """Positions from ``gaiadr3.gaia_source`` (``vari_microlensing`` has none), cached."""
        path = self._cache / "positions.csv"
        if not path.exists():
            ids = ",".join(self.events()["event_id"])
            q = f"SELECT source_id, ra, dec FROM gaiadr3.gaia_source WHERE source_id IN ({ids})"
            _write_atomic(path, _tap_csv(q).encode())
        t = Table.read(path, format="ascii.csv")
        return {str(r["source_id"]): (float(r["ra"]), float(r["dec"])) for r in t}

    def _lc_path(self, event_id: str) -> Path:
        return self._cache / f"{event_id}.csv"

    def _fetch_batch(self, ids: list[str]) -> None:
        body = urllib.parse.urlencode(
            {
                "RETRIEVAL_TYPE": "EPOCH_PHOTOMETRY",
                "DATA_STRUCTURE": "INDIVIDUAL",
                "FORMAT": "CSV",
                "VALID_DATA": "false",
                "ID": ",".join(ids),
            }
        ).encode()
        for attempt in range(3):
            try:
                with urllib.request.urlopen(DATALINK, data=body, timeout=600) as resp:  # noqa: S310
                    raw = resp.read()
                break
            except (OSError, http.client.HTTPException):  # truncated chunked replies happen
                if attempt == 2:
                    raise
                time.sleep(5 * (attempt + 1))
        if not zipfile.is_zipfile(io.BytesIO(raw)):  # one id: the service answers with bare CSV
            if len(ids) == 1 and raw.startswith(b"solution_id,"):
                _write_atomic(self._lc_path(ids[0]), raw)
            else:
                print(f"{self.name}: no epoch photometry for {ids}: {raw[:200]!r}", file=sys.stderr)
            return
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            for name in z.namelist():
                m = re.search(r"(\d{6,})\.csv$", name)
                if m:
                    _write_atomic(self._lc_path(m.group(1)), z.read(name))

    def prefetch(self, ids=None) -> int:
        """Fetch every missing light curve in parallel batches; returns the number now cached."""
        ids = list(ids if ids is not None else self.events()["event_id"])
        todo = [i for i in ids if not self._lc_path(i).exists()]
        batches = [todo[k : k + BATCH] for k in range(0, len(todo), BATCH)]
        with ThreadPoolExecutor(WORKERS) as pool:
            list(pool.map(self._fetch_batch, batches))
        return sum(self._lc_path(i).exists() for i in todo)

    # ------------------------------------------------------------------ protocol
    def light_curve(self, event_id: str) -> Table:
        path = self._lc_path(event_id)
        if not path.exists():
            self._fetch_batch([event_id])
        if not path.exists():
            raise KeyError(f"{self.name}: no epoch photometry for {event_id!r}")
        return parse_epoch_csv(path.read_text())

    def efficiency(self, t_e_days: float) -> float | None:
        """None: only completeness against OGLE-IV is published (Fig. 10), no table."""
        return None
