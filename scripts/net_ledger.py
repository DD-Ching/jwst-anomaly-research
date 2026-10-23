"""Persistent, shared download ledger with reservation (owner rule 2026-10-23, after #474).

Every network read of a unit goes through `fetch`, which reserves the expected size before the
request, counts the bytes actually received (response header + body + a fixed per-request TLS
allowance) while streaming, and cuts the transfer off when the cap would be exceeded. Retries and
repeated requests are separate reservations. Downloads made by other tools (e.g. dependency
installs) enter as `external` records measured by the caller.

The ledger is an append-only JSONL file under a file lock (fcntl; msvcrt on Windows), so
processes share one total and a restart does not reset it: a reservation that was never closed
(crashed or killed process) keeps counting at its reserved size.

Limit: bytes the operating system has buffered but `fetch` has not read when it cuts a transfer off
(up to the TCP receive window) are not counted.

    python scripts/net_ledger.py status LEDGER.jsonl CAP_BYTES
"""

from __future__ import annotations

import json
import os
import ssl
import sys
import time
import urllib.error
import urllib.request
import uuid
from contextlib import contextmanager
from pathlib import Path

try:  # POSIX
    import fcntl

    def _lock(fh) -> None:
        fcntl.flock(fh, fcntl.LOCK_EX)

    def _unlock(fh) -> None:
        fcntl.flock(fh, fcntl.LOCK_UN)

except ImportError:  # Windows (the owner's machine): lock byte 0 of the ledger file
    import msvcrt

    def _lock(fh, wait_s: float = 300.0) -> None:
        fh.seek(0)
        t0 = time.monotonic()
        while True:
            try:
                msvcrt.locking(fh.fileno(), msvcrt.LK_LOCK, 1)
                return
            except OSError:  # LK_LOCK gives up after ~10 s; retry until wait_s, then fail
                if time.monotonic() - t0 > wait_s:
                    raise
                time.sleep(0.05)

    def _unlock(fh) -> None:
        fh.seek(0)
        msvcrt.locking(fh.fileno(), msvcrt.LK_UNLCK, 1)


TLS_ALLOWANCE = (
    8192  # ASSUMPTION: bytes per request for TLS handshake / framing, not visible to urllib
)
HEAD_ALLOWANCE = 4096  # ASSUMPTION: bytes reserved for one response header
CHUNK = 65536


class BudgetExceeded(RuntimeError):
    """A reservation or a streamed transfer would exceed the ledger's cap."""


class Ledger:
    def __init__(self, path: str | Path, cap_bytes: int):
        self.path = Path(path)
        self.cap = int(cap_bytes)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.touch(exist_ok=True)

    @contextmanager
    def _locked(self):
        with open(self.path, "a+", encoding="utf-8") as fh:
            _lock(fh)
            try:
                fh.seek(0)
                yield fh
            finally:
                fh.flush()
                os.fsync(fh.fileno())
                _unlock(fh)

    @staticmethod
    def _used(lines: list[str]) -> int:
        reserved: dict[str, int] = {}
        closed = 0
        for ln in lines:
            if not ln.strip():
                continue
            r = json.loads(ln)
            if r["ev"] == "reserve":
                reserved[r["id"]] = reserved.get(r["id"], 0) + r["bytes"]
            elif r["ev"] == "extend":
                reserved[r["id"]] = reserved.get(r["id"], 0) + r["bytes"]
            elif r["ev"] == "done":
                reserved.pop(r["id"], None)
                closed += r["counted"]
            elif r["ev"] == "external":
                closed += r["bytes"]
        return closed + sum(reserved.values())

    def used(self) -> int:
        with self._locked() as fh:
            return self._used(fh.readlines())

    def _append(self, fh, rec: dict) -> None:
        rec = {"t": round(time.time(), 3), "pid": os.getpid(), **rec}
        fh.write(json.dumps(rec, sort_keys=True) + "\n")

    def reserve(self, nbytes: int, what: str) -> str:
        rid = uuid.uuid4().hex[:12]
        with self._locked() as fh:
            used = self._used(fh.readlines())
            if used + nbytes > self.cap:
                self._append(fh, {"ev": "refused", "bytes": int(nbytes), "what": what})
                raise BudgetExceeded(
                    f"reserve {nbytes} B for {what}: used {used} + {nbytes} > cap {self.cap}"
                )
            self._append(fh, {"ev": "reserve", "id": rid, "bytes": int(nbytes), "what": what})
        return rid

    def extend(self, rid: str, nbytes: int) -> None:
        with self._locked() as fh:
            used = self._used(fh.readlines())
            if used + nbytes > self.cap:
                raise BudgetExceeded(
                    f"extend {rid} by {nbytes} B: used {used} + {nbytes} > cap {self.cap}"
                )
            self._append(fh, {"ev": "extend", "id": rid, "bytes": int(nbytes)})

    def done(self, rid: str, counted: int, status: str) -> None:
        with self._locked() as fh:
            self._append(fh, {"ev": "done", "id": rid, "counted": int(counted), "status": status})

    def external(self, nbytes: int, what: str) -> None:
        """A download measured outside `fetch`; charged even past the cap (caller checks first)."""
        with self._locked() as fh:
            self._append(fh, {"ev": "external", "bytes": int(nbytes), "what": what})


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    """Redirects are not followed (urllib would read their body uncounted); they raise HTTPError."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


def _opener() -> urllib.request.OpenerDirector:
    # NET_LEDGER_CAFILE: a combined CA file (e.g. the proxy bundle + certifi's Mozilla roots) for
    # hosts whose chain ends at a root the system bundle lacks. Verification is never disabled.
    cafile = os.environ.get("NET_LEDGER_CAFILE") or os.environ.get("SSL_CERT_FILE") or None
    ctx = ssl.create_default_context(cafile=cafile)
    return urllib.request.build_opener(urllib.request.HTTPSHandler(context=ctx), _NoRedirect)


def fetch(
    url: str,
    ledger: Ledger,
    est_bytes: int,
    byte_range: tuple[int, int] | None = None,
    retries: int = 1,
    timeout: float = 60.0,
) -> bytes:
    """GET `url` (optionally bytes [a, b] inclusive) through the ledger; returns the body.

    `est_bytes` is the expected body size; it is reserved together with the header and TLS
    allowances. If the body grows beyond the reservation, more is reserved chunk by chunk before it
    is read, and the transfer is cut off with BudgetExceeded when the cap would be passed. Each
    attempt (including retries) is a separate reservation.
    """
    if byte_range is not None:
        est_bytes = byte_range[1] - byte_range[0] + 1
    last: Exception | None = None
    for attempt in range(retries + 1):
        reserved = int(est_bytes) + HEAD_ALLOWANCE + TLS_ALLOWANCE
        rid = ledger.reserve(reserved, f"GET {url} range={byte_range} attempt={attempt}")
        got = TLS_ALLOWANCE
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "jwst-anomaly-research"})
            if byte_range is not None:
                req.add_header("Range", f"bytes={byte_range[0]}-{byte_range[1]}")
            with _opener().open(req, timeout=timeout) as resp:
                got += len(str(resp.headers).encode()) + 32
                parts = []
                while True:
                    if (
                        got + CHUNK > reserved
                    ):  # reserve before reading, so counted bytes never pass the cap
                        ledger.extend(rid, got + CHUNK - reserved)
                        reserved = got + CHUNK
                    buf = resp.read(CHUNK)
                    if not buf:
                        break
                    got += len(buf)
                    parts.append(buf)
            ledger.done(rid, got, "ok")
            return b"".join(parts)
        except BudgetExceeded:
            ledger.done(rid, got, "cut: budget")
            raise
        except urllib.error.HTTPError as e:  # error or redirect: count its header and body
            body = e.read() if e.fp is not None else b""
            got += len(str(e.headers).encode()) + 32 + len(body)
            ledger.done(rid, got, f"error: HTTP {e.code}")
            last = e
        except Exception as e:  # network error: count what arrived, then retry
            ledger.done(rid, max(got, TLS_ALLOWANCE), f"error: {type(e).__name__}")
            last = e
    raise RuntimeError(f"fetch failed after {retries + 1} attempts: {last}")


if __name__ == "__main__":
    if sys.argv[1] == "status":
        lg = Ledger(sys.argv[2], int(sys.argv[3]))
        print(json.dumps({"used_bytes": lg.used(), "cap_bytes": lg.cap}))
