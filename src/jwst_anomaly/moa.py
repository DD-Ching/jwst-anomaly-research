"""MOA-II 9-year Galactic-bulge light curves as a :class:`signatures.LightCurveSurvey` (D-062).

The NASA Exoplanet Archive serves the MOA-II 2006–2014 release: 2,409,061 variable objects found
by Cut-0 (Koshimoto et al. 2023, arXiv:2303.08279, Table 2: DAOFIND on difference images, positive
**and negative** PSF profiles, S/N > 2.7, N_continue,8 ≥ 3, σ_x,y ≤ 1 px), each with a
difference-imaging light curve, before any bump or PSPL cut. The archive's own page describes an
older Cut-0 (Sumi et al. 2011: S/N > 5, max N_detect,continue > 2); the 2,409,061 count is
Nunota et al. 2024 (arXiv:2410.23553), Sect. 2.2.

Files (pinned by sha256 in ``FILES``; fetched with ``photometry.fetch_catalog``):

- ``bulk/metadata.ipac.tar.gz``: one row per object (field, chip, subframe, id, x, y, RA/Dec); the
  published selection statistics and PSPL/FSPL fits are filled only for the ~6,000 candidates;
- ``bulk/gb<F>.tar``: an uncompressed tar of ``gb<F>-R-<chip>-<subframe>-<id>.ipac.gz`` IPAC
  tables. Only one field is pinned (gb22, the smallest: 18,599 objects, 3.5 GB). Members are read
  in place through a header index (no extraction to disk).

Light-curve columns (https://exoplanetarchive.ipac.caltech.edu/docs/API_moa_columns.html): ``HJD``
(HJD − 2450000), ``flux`` (difference flux in counts, relative to the reference image; a MOA-Red
magnitude of 20 is 691.8 counts on chip 2 and 1,445 counts on every other chip), ``cor_flux``
(the same, de-trended against seeing and airmass, Koshimoto et al. 2023 Eq. 4; present for a
subset of objects), ``flux_err``, ``included`` (False = bad-quality epoch) and observing-condition
columns (fwhm, sky, airmass, ...). Everything returned here is ``observed`` as published.
"""

from __future__ import annotations

import gzip
import io
import tarfile
import warnings
import zlib
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack

from jwst_anomaly import paths, photometry, schema
from jwst_anomaly.signatures import standard_flux_light_curve

BULK = "https://exoplanetarchive.ipac.caltech.edu/data/Contributed/MOA/bulk/"
COLUMNS_DOC = "https://exoplanetarchive.ipac.caltech.edu/docs/API_moa_columns.html"
REFERENCE = (
    "MOA-II 9-year release, NASA Exoplanet Archive; Koshimoto et al. 2023 (arXiv:2303.08279)"
)

# url -> (sha256, size in bytes), retrieved 2026-10-08 (SOURCES.md "MOA-II 9-year bulge release").
# gb22.tar is 3.5 GB: stated reason (CLAUDE.md) — light curves exist only as per-field tars, and
# gb22 is the smallest field.
FILES: dict[str, tuple[str, int]] = {
    BULK + "metadata.ipac.tar.gz": (
        "b339ec176933f4bfcad09d6e17f08166c7a18ffea2ddc567b6cca83fa25fe55c",
        97332954,
    ),
    BULK + "gb22.tar": (
        "1cb0173e676915dcf602647e5a2dc315f0ab51c8f2a612ffcdcab4e44edad2d0",
        3510138880,
    ),
}

HJD_OFFSET = 2450000.0  # light-curve HJD column is HJD − 2450000
T_START, T_END = 2453824.0, 2456970.0  # HJD span of the release (Nunota et al. 2024, Sect. 2)
N_CUT0_TOTAL = 2_409_061  # variable objects passing Cut-0, all 22 fields (Nunota et al. 2024)
BAND = "MOA-Red"


# Cut-0 objects per field, counted in metadata.ipac (sha256 above; observed). Sum: 2,409,061.
CUT0_PER_FIELD = {
    1: 109715, 2: 95624, 3: 167431, 4: 205886, 5: 259556, 6: 35948, 7: 75328, 8: 120726,
    9: 237123, 10: 181820, 11: 64090, 12: 73386, 13: 100150, 14: 182302, 15: 83145, 16: 64353,
    17: 100448, 18: 83855, 19: 53861, 20: 52005, 21: 43710, 22: 18599,
}  # fmt: skip
# Nunota et al. 2024, Table 1: (subfields used of 80, N_s = source stars with 10 ≤ I_s ≤ 21.4).
# gb6 and gb22 are not in it (no clear red-clump population in the CMD).
NUNOTA_NS = {
    1: (79, 21047010), 2: (79, 17647488), 3: (79, 22711037), 4: (77, 25985143),
    5: (65, 29137851), 7: (78, 16344845), 8: (78, 22263658), 9: (79, 33308780),
    10: (70, 21465124), 11: (76, 10931979), 12: (79, 16090446), 13: (79, 23728248),
    14: (79, 25094851), 15: (62, 10404096), 16: (79, 15028199), 17: (79, 17979175),
    18: (78, 15478756), 19: (78, 12087586), 20: (79, 10730219), 21: (73, 8321495),
}  # fmt: skip


# Field tar sizes in bytes (HTTP HEAD Content-Length, 2026-10-08; Last-Modified 2023-10-10 … 13).
# Only gb22 is downloaded whole; the others are streamed with byte-range reads (``moa_stream``).
TAR_BYTES = {
    1: 70145955840, 2: 60503715840, 3: 214532423680, 4: 280555120640, 5: 508495472640,
    6: 7730944000, 7: 30745886720, 8: 77880094720, 9: 430753607680, 10: 220873328640,
    11: 22442977280, 12: 27467520000, 13: 53574082560, 14: 212026122240, 15: 32830392320,
    16: 22349178880, 17: 52299028480, 18: 41926410240, 19: 16264140800, 20: 15218606080,
    21: 12313610240, 22: 3510138880,
}  # fmt: skip


TAR_LAST_MODIFIED = {  # HTTP Last-Modified of each field tar (HEAD, 2026-10-08)
    1: "2023-10-13", 2: "2023-10-11", 3: "2023-10-12", 4: "2023-10-13", 5: "2023-10-13",
    6: "2023-10-12", 7: "2023-10-12", 8: "2023-10-11", 9: "2023-10-13", 10: "2023-10-12",
    11: "2023-10-11", 12: "2023-10-11", 13: "2023-10-12", 14: "2023-10-12", 15: "2023-10-13",
    16: "2023-10-12", 17: "2023-10-13", 18: "2023-10-12", 19: "2023-10-12", 20: "2023-10-10",
    21: "2023-10-10", 22: "2023-10-10",
}  # fmt: skip


def tar_url(field: int | str) -> str:
    return f"{BULK}gb{int(str(field).removeprefix('gb'))}.tar"


def object_url(event_id: str) -> str:
    """Per-object uncompressed IPAC file (used by the archive viewer; undocumented, SOURCES.md)."""
    f, chip, _s, _i = parse_event_id(event_id)
    return f"https://exoplanetarchive.ipac.caltech.edu/data/Contributed/MOA/gb{f}/R/{chip}/{event_id}.ipac"


def tar_header(block: bytes):
    """``tarfile.TarInfo`` of a 512-byte ustar header block, or None (bad checksum, no magic,
    an end-of-archive zero block). Used to resynchronise on member boundaries inside a byte range
    of an uncompressed tar (gzip member data never passes the checksum and magic by chance in
    practice; the per-field member count is checked against the metadata)."""
    if len(block) < 512 or block[257:262] != b"ustar":
        return None
    try:
        return tarfile.TarInfo.frombuf(block[:512], "utf-8", "surrogateescape")
    except tarfile.HeaderError:
        return None


def member_id(name: str) -> str | None:
    base = name.rsplit("/", 1)[-1]
    return base[: -len(".ipac.gz")] if base.endswith(".ipac.gz") else None


def walk_members(buf, base: int, start: int, stop: int):
    """Light-curve members whose header starts in ``[start, stop)`` (absolute tar offsets) of the
    bytes ``buf`` that begin at absolute offset ``base``: yields ``(event_id, offset_data, size)``
    with absolute offsets, or ``("", next_header_offset, needed_end)`` once a member's data runs
    past ``buf`` (the caller extends ``buf`` and resumes from that header). ``start`` need not be
    on a header: the first valid ustar header at or after it (512-aligned) is used."""
    off = start + (-start) % 512
    while off < stop and tar_header(bytes(buf[off - base : off - base + 512])) is None:
        if off - base + 512 > len(buf):
            return  # no header in the rest of the range (end of archive)
        off += 512
    while off < stop:
        if off - base + 512 > len(buf):
            yield "", off, off + 512
            return
        info = tar_header(bytes(buf[off - base : off - base + 512]))
        if info is None:  # end of archive (zero blocks)
            return
        data = off + 512
        end = data + info.size
        if end > base + len(buf):
            yield "", off, end
            return
        if info.type == tarfile.XGLTYPE:  # pax global header: names nothing
            pass
        elif info.type in (tarfile.GNUTYPE_LONGNAME, tarfile.XHDTYPE):
            # a name header applies to the next header, which may start in the next range,
            # where the resync cannot see it; MOA member names fit in 100 bytes, so refuse
            raise ValueError(f"tar name extension header at offset {off} is not supported")
        else:
            eid = member_id(info.name) if info.isfile() else None
            if eid:
                yield eid, data, info.size
        off = data + info.size + (-info.size) % 512


def stars_per_cut0_object() -> np.ndarray:
    """N_s / (Cut-0 objects × subfields used / 80) for the 20 fields of Nunota et al. (derived)."""
    return np.array([ns / (CUT0_PER_FIELD[f] * nsub / 80.0) for f, (nsub, ns) in NUNOTA_NS.items()])


def star_count_estimate(n_cut0: int) -> tuple[float, float, float]:
    """Monitored source stars (10 ≤ I_s ≤ 21.4) of a field without a published N_s, from its number
    of Cut-0 objects: (median, min, max) over the 20 published fields of N_s per Cut-0 object.
    ASSUMPTION: the ratio holds for the field (gb22 is not in Nunota et al.'s table)."""
    r = stars_per_cut0_object()
    return float(np.median(r) * n_cut0), float(r.min() * n_cut0), float(r.max() * n_cut0)


def counts_at_mag20(chip: int) -> float:
    """Difference-flux counts of a MOA-Red magnitude of 20 (archive column documentation)."""
    return 691.8 if int(chip) == 2 else 1445.0


def mag_to_counts(mag, chip: int):
    return counts_at_mag20(chip) * 10.0 ** (-0.4 * (np.asarray(mag, float) - 20.0))


def counts_to_mag(counts, chip: int):
    return 20.0 - 2.5 * np.log10(np.asarray(counts, float) / counts_at_mag20(chip))


def parse_event_id(event_id: str) -> tuple[int, int, int, int]:
    """``"gb22-R-8-5-98"`` -> (field 22, chip 8, subframe 5, id 98)."""
    parts = event_id.split("-")
    if len(parts) != 5 or not parts[0].startswith("gb") or parts[1] != "R":
        raise ValueError(f"not a MOA object id: {event_id!r}")
    return int(parts[0][2:]), int(parts[2]), int(parts[3]), int(parts[4])


def make_event_id(field: int, chip: int, subframe: int, oid: int) -> str:
    return f"gb{int(field)}-R-{int(chip)}-{int(subframe)}-{int(oid)}"


def _ipac_header(lines: list[str]) -> list[str]:
    return [h.strip() for h in lines[0].strip().strip("|").split("|")]


def parse_lightcurve(raw: bytes | str, columns=None) -> dict[str, np.ndarray]:
    """Columns of one MOA light-curve IPAC table as arrays (``included`` as bool).

    The published tables carry one all-"nan" row after the header; "nan" values stay NaN.
    ``columns``: convert only these (the rest are skipped; ``HJD`` is always read). The fast path
    splits the whole body into tokens once (C level) and converts each wanted column with numpy;
    it gives bit-identical arrays to the line-by-line path, which stays as the fallback for
    ragged tables.
    """
    if isinstance(raw, str):
        raw = raw.encode()
    if bytes(raw[:2]) == b"\x1f\x8b":
        d = zlib.decompressobj(wbits=31)  # one gzip member, CRC and length checked in C
        out = d.decompress(raw)
        raw = out if d.eof and not d.unused_data else gzip.decompress(raw)
    elif isinstance(raw, memoryview):
        raw = raw.tobytes()
    fast = _parse_fixed(raw, columns)
    if fast is not None:
        return fast
    fast = _parse_tokens(raw, columns)
    if fast is not None:
        return fast
    out = _parse_lines(raw.decode())
    if columns is not None:
        want = {"HJD", *columns}
        out = {k: v for k, v in out.items() if k in want}
    return out


def _parse_fixed(raw: bytes, columns=None) -> dict[str, np.ndarray] | None:
    """Fixed-width parse: only the wanted columns' byte ranges are converted. None unless every
    row has the header's width, ends in a newline and is blank under every header bar (then each
    column slice holds exactly the token the whitespace split would give)."""
    pos, bars, names = 0, None, None
    while raw[pos : pos + 1] in (b"|", b"\\"):
        end = raw.find(b"\n", pos)
        if end < 0:
            return None
        if names is None and raw[pos : pos + 1] == b"|":
            line = raw[pos:end].rstrip(b"\r")
            bars = [i for i, c in enumerate(line) if c == 124]
            names = [h.strip().decode() for h in line.strip(b"|").split(b"|")]
        pos = end + 1
    body = raw[pos:]
    width = body.find(b"\n") + 1
    if names is None or len(bars) != len(names) + 1 or width <= bars[-1] or len(body) % width:
        return None
    u8 = np.frombuffer(body, np.uint8).reshape(-1, width)
    if (
        (u8[:, -1] != 10).any()
        or (u8[:, bars[:-1]] != 32).any()
        or (u8[:, bars[-1] : -1] != 32).any()
    ):
        return None
    want = None if columns is None else {"HJD", *columns}
    out: dict[str, np.ndarray] = {}
    with warnings.catch_warnings():  # once per member, not per column
        warnings.simplefilter("error", DeprecationWarning)
        for j, name in enumerate(names):
            if want is not None and name not in want:
                continue
            a, b = bars[j] + 1, bars[j + 1]
            col = np.ascontiguousarray(u8[:, a:b])
            if name == "included":
                out[name] = np.char.strip(col.view(f"S{b - a}").ravel()) == b"True"
            else:  # each slice ends in a blank (the bar column), so tokens stay apart
                try:  # an unparseable token: fall back to the token parser, never raise
                    v = np.fromstring(np.ascontiguousarray(u8[:, a : b + 1]).tobytes(), sep=" ")
                except (DeprecationWarning, ValueError):
                    return None
                if v.size != len(u8):
                    return None
                out[name] = v
    keep = np.isfinite(out["HJD"]) if "HJD" in out else np.ones(len(u8), bool)
    return {k: v[keep] for k, v in out.items()}


def _parse_tokens(raw: bytes, columns=None) -> dict[str, np.ndarray] | None:
    """Whole-member parse (None when the body is not a full rectangle of tokens)."""
    pos, names = 0, None
    while raw[pos : pos + 1] in (b"|", b"\\"):
        end = raw.find(b"\n", pos)
        if end < 0:
            return None
        if names is None and raw[pos : pos + 1] == b"|":
            names = [h.strip().decode() for h in raw[pos:end].strip().strip(b"|").split(b"|")]
        pos = end + 1
    if names is None:
        raise ValueError("no IPAC header")
    toks = raw[pos:].split()
    nc = len(names)
    if not toks or len(toks) % nc:
        return None
    want = None if columns is None else {"HJD", *columns}
    out: dict[str, np.ndarray] = {}
    for j, name in enumerate(names):
        if want is not None and name not in want:
            continue
        col = np.array(toks[j::nc])
        out[name] = col == b"True" if name == "included" else col.astype(float)
    keep = np.isfinite(out["HJD"]) if "HJD" in out else np.ones(len(toks) // nc, bool)
    return {k: v[keep] for k, v in out.items()}


def _parse_lines(raw: str) -> dict[str, np.ndarray]:
    lines = raw.splitlines()
    head = [ln for ln in lines if ln.startswith("|")]
    if not head:
        raise ValueError("no IPAC header")
    names = _ipac_header(head)
    rows = [ln.split() for ln in lines if ln.strip() and not ln.startswith(("|", "\\"))]
    rows = [r for r in rows if len(r) == len(names)]
    arr = np.array(rows, dtype=str).reshape(-1, len(names))
    out: dict[str, np.ndarray] = {}
    for j, name in enumerate(names):
        col = arr[:, j]
        if name == "included":
            out[name] = col == "True"
        else:
            out[name] = col.astype(float)
    keep = np.isfinite(out["HJD"]) if "HJD" in out else np.ones(len(arr), bool)
    return {k: v[keep] for k, v in out.items()}


def select_flux(cols: dict[str, np.ndarray]) -> tuple[np.ndarray, ...]:
    """``(HJD, flux, flux_err, kind)`` of the included epochs.

    ``cor_flux`` (de-trended) when it is finite on every included epoch, else ``flux``. The
    de-trended flux has a different zero point (its baseline is near 0); both are counts.
    """
    inc = cols["included"] & np.isfinite(cols["flux"]) & np.isfinite(cols["flux_err"])
    cor = cols.get("cor_flux")
    if cor is not None and inc.any() and np.all(np.isfinite(cor[inc])):
        f, kind = cor, "detrended difference"
    else:
        f, kind = cols["flux"], "difference"
    t = cols["HJD"][inc] + HJD_OFFSET
    return t, f[inc], cols["flux_err"][inc], kind


def parse_metadata(lines, field: int) -> Table:
    """Rows of ``metadata.ipac`` for one field (``lines``: an iterable of text lines)."""
    return split_metadata(lines, [field])[int(field)]


def split_metadata(lines, fields, block: int = 50_000) -> dict[int, Table]:
    """Rows of ``metadata.ipac`` for several fields in one pass (``observed``); rows are converted
    to columns in blocks so the 2.9 GB table is never held as text."""
    names: list[str] | None = None
    types: list[str] | None = None
    want = {str(int(f)): int(f) for f in fields}
    pending: dict[int, list] = {f: [] for f in want.values()}
    parts: dict[int, list] = {f: [] for f in want.values()}

    def flush(f):
        rows = [r for r in pending[f] if len(r) == len(names)]
        pending[f] = []
        arr = np.array(rows, dtype=str).reshape(-1, len(names))
        out = {}
        for j, (name, typ) in enumerate(zip(names, types or ["double"] * len(names), strict=True)):
            col = arr[:, j]
            out[name] = col if typ == "char" else np.where(col == "null", "nan", col).astype(float)
        parts[f].append(Table(out))

    for ln in lines:
        if ln.startswith("|"):
            if names is None:
                names = _ipac_header([ln])
            elif types is None:
                types = _ipac_header([ln])
            continue
        if ln.startswith("\\") or not ln.strip():
            continue
        tok = ln.split()
        f = want.get(tok[0])
        if f is None:
            continue
        pending[f].append(tok)
        if len(pending[f]) >= block:
            flush(f)
    if names is None:
        raise ValueError("no IPAC header in metadata")
    types = types or ["double"] * len(names)
    out = {}
    for f in want.values():
        flush(f)
        t = vstack(parts[f]) if len(parts[f]) > 1 else parts[f][0]
        for c in ("field", "chip", "subframe", "id"):
            t[c] = t[c].astype(int)
        t.meta.update(
            provenance=schema.Provenance.OBSERVED.value,
            source=f"{BULK}metadata.ipac.tar.gz, field gb{f} rows ({REFERENCE})",
        )
        out[f] = t
    return out


def write_metadata_caches(metadata_path: Path, derived_dir: Path, fields=None) -> list[Path]:
    """One pass over ``metadata.ipac`` writing ``metadata_gb<F>.ecsv`` for ``fields`` (default all
    22; ``observed``)."""
    fields = sorted(CUT0_PER_FIELD) if fields is None else [int(f) for f in fields]
    with tarfile.open(metadata_path, "r:gz") as tar:
        member = next(m for m in tar if m.isfile() and m.name.endswith("metadata.ipac"))
        fh = io.TextIOWrapper(tar.extractfile(member), encoding="ascii")
        tabs = split_metadata(fh, fields)
    derived_dir.mkdir(parents=True, exist_ok=True)
    paths_out = []
    for f, t in tabs.items():
        t["event_id"] = [
            make_event_id(*r)
            for r in zip(t["field"], t["chip"], t["subframe"], t["id"], strict=True)
        ]
        t.rename_columns(["ra_j2000", "dec_j2000"], ["ra", "dec"])
        t = t[
            [
                "event_id",
                "ra",
                "dec",
                *[c for c in t.colnames if c not in ("event_id", "ra", "dec")],
            ]
        ]
        t.meta.update(
            provenance=schema.Provenance.OBSERVED.value,
            source=f"{BULK}metadata.ipac.tar.gz, field gb{f} rows ({REFERENCE})",
        )
        path = derived_dir / f"metadata_gb{f}.ecsv"
        t.write(path, overwrite=True)
        paths_out.append(path)
    return paths_out


class MoaField:
    """Adapter for one MOA-II field: events = every Cut-0 object, light curves in counts.

    ``tar_path`` / ``metadata_path`` override the pinned downloads (tests, other mirrors).
    """

    def __init__(
        self,
        field: str = "gb22",
        raw_dir: Path | None = None,
        derived_dir: Path | None = None,
        tar_path: Path | None = None,
        metadata_path: Path | None = None,
    ):
        if not field.startswith("gb"):
            raise ValueError(f"MOA fields are gb1 … gb22, not {field!r}")
        self.field = field
        self.field_number = int(field[2:])
        self.name = f"moa2-9yr-{field}"
        self._raw_dir = raw_dir or paths.data_root() / "raw" / "moa"
        self._derived_dir = derived_dir or paths.data_root() / "derived" / "moa"
        self._tar_path = tar_path
        self._metadata_path = metadata_path
        self._index: dict[str, tuple[int, int]] | None = None
        self._events: Table | None = None

    # ------------------------------------------------------------------ files
    def _fetch(self, name: str) -> Path:
        url = BULK + name
        if url not in FILES:
            raise KeyError(f"{url} is not pinned in moa.FILES (only one pilot field is)")
        sha, size = FILES[url]
        return photometry.fetch_catalog(url, sha, cache_dir=self._raw_dir, max_bytes=size)

    def tar_path(self) -> Path:
        if self._tar_path is None:
            self._tar_path = self._fetch(f"{self.field}.tar")
        return self._tar_path

    def metadata_path(self) -> Path:
        if self._metadata_path is None:
            self._metadata_path = self._fetch("metadata.ipac.tar.gz")
        return self._metadata_path

    def index(self) -> dict[str, tuple[int, int]]:
        """``event_id -> (byte offset, size)`` of each light-curve member (tar headers only)."""
        if self._index is None:
            idx = {}
            with tarfile.open(self.tar_path(), "r:") as tar:
                for m in tar:
                    base = m.name.rsplit("/", 1)[-1]
                    if m.isfile() and base.endswith(".ipac.gz"):
                        idx[base[: -len(".ipac.gz")]] = (m.offset_data, m.size)
            self._index = idx
        return self._index

    # ------------------------------------------------------------------ protocol
    def events(self) -> Table:
        """Every Cut-0 object of the field (``observed``): ``event_id, ra, dec``, chip, subframe,
        id, pixel x/y and the published statistics and fits (NaN except for the ~6,000 candidates).
        The field's rows are cached as ECSV under ``derived/moa/`` (the metadata is 2.9 GB)."""
        if self._events is not None:
            return self._events
        cache = self._derived_dir / f"metadata_{self.field}.ecsv"
        if not cache.exists():
            write_metadata_caches(self.metadata_path(), self._derived_dir, [self.field_number])
        t = Table.read(cache)
        self._events = t
        return t

    def read_member(self, event_id: str, fh=None) -> dict[str, np.ndarray]:
        """All columns of one light curve (``observed``)."""
        try:
            off, size = self.index()[event_id]
        except KeyError:
            raise KeyError(f"{self.name}: no light curve for {event_id!r}") from None
        close = fh is None
        fh = fh or open(self.tar_path(), "rb")  # noqa: SIM115
        try:
            fh.seek(off)
            raw = fh.read(size)
        finally:
            if close:
                fh.close()
        return parse_lightcurve(raw)

    def light_curve(self, event_id: str) -> Table:
        cols = self.read_member(event_id)
        t, f, sf, kind = select_flux(cols)
        out = standard_flux_light_curve(
            t,
            f,
            sf,
            BAND,
            source=f"{BULK}{self.field}.tar:{event_id}.ipac.gz ({REFERENCE})",
            time_system="HJD",
            flux_unit=f"counts (MOA-Red 20 mag = {counts_at_mag20(parse_event_id(event_id)[1])})",
            flux_kind=kind,
        )
        out.meta["n_excluded"] = int((~cols["included"]).sum())
        return out

    def iter_light_curves(self, ids=None):
        """``(event_id, columns)`` for every light curve in tar order (one sequential pass)."""
        want = None if ids is None else set(ids)
        items = sorted(self.index().items(), key=lambda kv: kv[1][0])
        with open(self.tar_path(), "rb") as fh:
            for eid, _ in items:
                if want is None or eid in want:
                    yield eid, self.read_member(eid, fh)

    def efficiency(self, t_e_days: float) -> float | None:
        """None: the release publishes no machine-readable efficiency (limits need injections)."""
        return None
