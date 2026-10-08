"""MOA-II 9-year Galactic-bulge light curves as a :class:`signatures.LightCurveSurvey` (D-060).

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
from pathlib import Path

import numpy as np
from astropy.table import Table

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


def parse_lightcurve(raw: bytes | str) -> dict[str, np.ndarray]:
    """Columns of one MOA light-curve IPAC table as arrays (``included`` as bool).

    The published tables carry one all-"nan" row after the header; "nan" values stay NaN.
    """
    if isinstance(raw, bytes):
        if raw[:2] == b"\x1f\x8b":
            raw = gzip.decompress(raw)
        raw = raw.decode()
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
    names: list[str] | None = None
    types: list[str] | None = None
    rows = []
    want = str(int(field))
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
        if tok[0] != want:
            continue
        rows.append(tok)
    if names is None:
        raise ValueError("no IPAC header in metadata")
    types = types or ["double"] * len(names)
    rows = [r for r in rows if len(r) == len(names)]
    arr = np.array(rows, dtype=str).reshape(-1, len(names))
    out = {}
    for j, (name, typ) in enumerate(zip(names, types, strict=True)):
        col = arr[:, j]
        if typ == "char":
            out[name] = col
        else:
            out[name] = np.where(col == "null", "nan", col).astype(float)
    t = Table(out)
    for c in ("field", "chip", "subframe", "id"):
        t[c] = t[c].astype(int)
    return t


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
        if cache.exists():
            t = Table.read(cache)
        else:
            with tarfile.open(self.metadata_path(), "r:gz") as tar:
                member = next(m for m in tar if m.isfile() and m.name.endswith("metadata.ipac"))
                fh = io.TextIOWrapper(tar.extractfile(member), encoding="ascii")
                t = parse_metadata(fh, self.field_number)
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
                source=f"{BULK}metadata.ipac.tar.gz, field {self.field} rows ({REFERENCE})",
            )
            self._derived_dir.mkdir(parents=True, exist_ok=True)
            t.write(cache, overwrite=True)
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
