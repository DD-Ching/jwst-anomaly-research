"""Stream MOA-II field tars with parallel HTTP byte-range reads (D-063).

The field tars (3.5–508 GB, uncompressed, ``Accept-Ranges: bytes``) are never written to disk:
a byte range ``[a, b)`` is downloaded, its light-curve members are found by walking the tar
headers (resynchronised on the first valid ustar header at or after ``a``; ``moa.walk_members``),
and each member's gzipped bytes go to the caller. A member belongs to the range its header starts
in, so consecutive ranges cover every member exactly once; the per-field member count is checked
against ``moa.CUT0_PER_FIELD`` by the caller.

``RangeReader`` reads ranges from a URL (one ``requests.Session`` per thread, retries with
exponential back-off on 429/5xx and connection errors) or from a local copy of the tar (tests,
gb22). Everything returned is ``observed`` data as published.
"""

from __future__ import annotations

import hashlib
import random
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from jwst_anomaly import moa

RETRY_STATUS = {429, 500, 502, 503, 504}


class RangeReader:
    """``read(a, b)`` -> bytes ``[a, b)`` of one tar, from ``url`` or a local ``path``."""

    def __init__(self, url: str | None = None, path: Path | None = None, retries: int = 8,
                 timeout: float = 120.0):  # fmt: skip
        if (url is None) == (path is None):
            raise ValueError("give exactly one of url and path")
        self.url, self.path = url, None if path is None else Path(path)
        self.retries, self.timeout = retries, timeout
        self._local = threading.local()
        self.bytes_read = 0
        self.n_retries = 0
        self._lock = threading.Lock()

    def _session(self):
        s = getattr(self._local, "session", None)
        if s is None:
            import requests

            s = self._local.session = requests.Session()
        return s

    def read(self, a: int, b: int) -> bytes:
        if b <= a:
            return b""
        if self.path is not None:
            with open(self.path, "rb") as fh:
                fh.seek(a)
                data = fh.read(b - a)
        else:
            data = self._http(a, b)
        with self._lock:
            self.bytes_read += len(data)
        return data

    def _http(self, a: int, b: int) -> bytes:
        import requests

        delay = 2.0
        for attempt in range(self.retries + 1):
            try:
                r = self._session().get(
                    self.url, headers={"Range": f"bytes={a}-{b - 1}"}, timeout=self.timeout
                )
                if r.status_code == 206 and len(r.content) == b - a:
                    return r.content
                if r.status_code == 200 and a == 0 and len(r.content) >= b:
                    return r.content[:b]  # server ignored the range on a small file
                if r.status_code not in RETRY_STATUS and r.status_code != 206:
                    raise OSError(f"{self.url} bytes {a}-{b - 1}: HTTP {r.status_code}")
            except (
                requests.ConnectionError,
                requests.Timeout,
                requests.exceptions.ChunkedEncodingError,
            ):
                pass
            if attempt == self.retries:
                break
            with self._lock:
                self.n_retries += 1
            time.sleep(delay * (1.0 + random.random()))
            delay = min(delay * 2.0, 120.0)
        raise OSError(f"{self.url} bytes {a}-{b - 1}: failed after {self.retries} retries")


def read_segment(
    reader: RangeReader, start: int, stop: int, total: int, digest=None
) -> list[tuple]:
    """``[(event_id, offset_data, size, gz_view)]`` of the members whose header starts in
    ``[start, stop)``; the last member's data is fetched past ``stop`` as needed. ``gz_view`` is a
    zero-copy ``memoryview`` into the downloaded range (one copy of the bytes in memory).
    ``digest``: a list that receives the sha256 of bytes ``[start, stop)`` (content pin)."""
    buf = reader.read(start, min(stop, total))
    if digest is not None:
        digest.append(hashlib.sha256(buf).hexdigest())
    base, pos, out = start, start, []
    while True:
        need = None
        view = memoryview(buf)
        for eid, off, size in moa.walk_members(view, base, pos, stop):
            if eid == "":
                need, pos = size, off  # resume at this header once buf reaches `need`
                break
            out.append((eid, off, size, view[off - base : off - base + size]))
        if need is None or base + len(buf) >= total:
            return out
        end = min(max(need, base + len(buf) + (1 << 20)), total)
        end += (-end) % 512 if end < total else 0
        buf = buf + reader.read(base + len(buf), min(end, total))  # rare: a member past `stop`


def read_segment_hashed(reader: RangeReader, start: int, stop: int, total: int):
    """``(members, sha256 of [start, stop))``."""
    d: list[str] = []
    return read_segment(reader, start, stop, total, d), d[0]


def range_digest(segment_sha256s) -> str:
    """Content pin of a streamed file: sha256 of its ordered 64 MiB range digests (hex, joined)."""
    return hashlib.sha256("".join(segment_sha256s).encode()).hexdigest()


def segments(start: int, stop: int, seg_bytes: int) -> list[tuple[int, int]]:
    """512-aligned ``[a, b)`` pieces of ``[start, stop)``."""
    seg_bytes -= seg_bytes % 512
    return [(a, min(a + seg_bytes, stop)) for a in range(start, stop, seg_bytes)]


def fetch_members(reader: RangeReader, items, conns: int = 12) -> dict[str, bytes]:
    """``{event_id: gz_bytes}`` for ``items`` = ``[(event_id, offset_data, size)]`` (concurrent
    single-member range reads; for the few hundred to few thousand light curves after the
    pre-screen)."""
    items = list(items)
    if reader.path is not None:  # one local file: sequential reads in offset order
        out = {}
        with open(reader.path, "rb") as fh:
            for eid, off, size in sorted(items, key=lambda x: x[1]):
                fh.seek(int(off))
                out[eid] = fh.read(int(size))
        return out
    with ThreadPoolExecutor(conns) as ex:
        data = ex.map(lambda it: reader.read(int(it[1]), int(it[1]) + int(it[2])), items)
        return {it[0]: d for it, d in zip(items, data, strict=True)}


def fetch_objects(event_ids, conns: int = 8) -> dict[str, bytes]:
    """Per-object uncompressed IPAC files (``moa.object_url``) for objects without a recorded tar
    offset (vetting neighbours)."""
    ids = list(event_ids)
    with ThreadPoolExecutor(conns) as ex:
        return dict(zip(ids, ex.map(_get_whole, ids), strict=True))


def _get_whole(eid: str) -> bytes:
    import requests

    delay = 2.0
    for attempt in range(7):
        try:
            r = requests.get(moa.object_url(eid), timeout=120)
            if r.status_code == 200:
                return r.content
            if r.status_code not in RETRY_STATUS:
                raise OSError(f"{moa.object_url(eid)}: HTTP {r.status_code}")
        except (requests.ConnectionError, requests.Timeout):
            pass
        if attempt < 6:
            time.sleep(delay * (1.0 + random.random()))
            delay *= 2.0
    raise OSError(f"{moa.object_url(eid)}: failed after retries")


def limit_heap_growth() -> None:
    """glibc only: serve allocations ≥ 1 MB by mmap (returned to the OS when freed) and cap the
    malloc arenas, so a process that cycles 64 MB ranges does not keep a fragmented heap."""
    try:
        import ctypes

        libc = ctypes.CDLL("libc.so.6")
        libc.mallopt(-3, 1 << 20)  # M_MMAP_THRESHOLD (also stops the dynamic threshold)
        libc.mallopt(-8, 2)  # M_ARENA_MAX
    except (OSError, AttributeError):
        pass


class CpuMeter:
    """Machine-wide busy fraction between calls (``/proc/stat``; None where unavailable)."""

    def __init__(self):
        self._last = self._read()

    @staticmethod
    def _read():
        try:
            with open("/proc/stat") as fh:
                v = [int(x) for x in fh.readline().split()[1:]]
            return sum(v), v[3] + v[4]  # total, idle + iowait
        except (OSError, ValueError, IndexError):
            return None

    def busy(self) -> float | None:
        now = self._read()
        if now is None or self._last is None:
            return None
        dt, di = now[0] - self._last[0], now[1] - self._last[1]
        self._last = now
        return 100.0 * (1.0 - di / dt) if dt > 0 else None
