"""Reproducible, cached acquisition with checksummed manifests. Owner: bootstrap unit 1.

Files are downloaded with ``astroquery.mast.Observations.download_file`` (DECISIONS.md D-002)
into a temporary ``*.part`` file next to the target, checked against the size MAST reports,
hashed (sha256) and only then moved into place with ``os.replace``, so a file at its final
path is always complete. When MAST's listed size is stale (a reprocessed product), the
download is accepted only if it matches the Content-Length the download service declares.
The manifest (ECSV, one row per ``dataURI``, sorted, LF line endings) is the tracked
reproducibility record; it is rewritten only when its content changes.
"""

from __future__ import annotations

import hashlib
import io
import logging
import os
import urllib.request
from collections.abc import Iterable, Mapping
from datetime import UTC, datetime
from pathlib import Path, PurePosixPath
from typing import Any
from urllib.parse import quote

import numpy as np
from astropy.io import fits
from astropy.table import Row, Table
from astroquery.mast import conf as mast_conf

from . import paths, query, schema

log = logging.getLogger(__name__)

DEFAULT_MANIFEST_NAME = "mast_products.ecsv"
# paths.cache_dir()/mast, relative to the data root (cache_dir() is data_root()/"cache")
CACHE_PREFIX = PurePosixPath(paths.cache_dir().relative_to(paths.data_root()).as_posix(), "mast")
MAX_DOWNLOAD_BYTES = 200_000_000  # CLAUDE.md: a single download >200 MB needs a stated reason
MANIFEST_SOURCE = (
    "MAST via astroquery.mast.Observations.download_file; size checked against MAST, "
    "sha256 computed locally, pipeline_version read from each file's own header"
)
_INT_COLUMNS = {"size"}  # every other manifest column is a string
_HASH_CHUNK = 1 << 20


class DownloadError(RuntimeError):
    """A product could not be downloaded and verified."""


def fetch_products(
    products: Table,
    data_root: Path | None = None,
    manifest_path: Path | None = None,
    *,
    max_size_bytes: int | None = MAX_DOWNLOAD_BYTES,
    verbose: bool = False,
) -> Table:
    """Download ``products`` (``schema.PRODUCT_COLUMNS``) into the data cache.

    Idempotent: files already present with a matching sha256 are not re-downloaded.
    Writes/updates a manifest (``schema.MANIFEST_COLUMNS``) and returns its rows for
    ``products``. ``data_root`` defaults to ``paths.data_root()``; ``manifest_path``
    defaults to a file under ``paths.manifests_dir()``.

    Files go to ``<data_root>/cache/mast/<obs_id>/<productFilename>``. The returned rows are
    sorted by ``dataURI`` and carry an extra ``status`` column: ``downloaded``, ``cached``
    (sha256 matches the manifest) or ``adopted`` (file present with the MAST size but not in
    the manifest yet; hashed and recorded without re-downloading). Products larger than
    ``max_size_bytes`` are not downloaded (``None`` disables this guard); products for which
    MAST reports no size are refused, since neither completeness nor the guard can be checked.
    A failed product does not stop the others: the manifest is updated with every success,
    then a :class:`DownloadError` lists the failures. Concurrent runs on one manifest are not
    supported (the manifest is re-read just before writing, which narrows but does not close
    the race; a later run re-adopts any file whose row was lost).
    """
    missing = [c for c in schema.PRODUCT_COLUMNS if c not in products.colnames]
    if missing:
        raise ValueError(f"products: missing required columns {missing}")
    root = Path(data_root) if data_root is not None else paths.data_root()
    if manifest_path is None:
        manifest_path = paths.manifests_dir() / DEFAULT_MANIFEST_NAME
    manifest_path = Path(manifest_path)

    # Validate every target path before touching the network.
    planned: dict[str, tuple[Row, PurePosixPath]] = {}
    for prod in products:
        uri = _cell(prod, "dataURI")
        planned.setdefault(uri, (prod, local_relpath(prod)))

    previous = {row["dataURI"]: row for row in _rows(read_manifest(manifest_path))}
    rows: list[dict[str, Any]] = []
    failures: list[str] = []
    for uri in sorted(planned):
        prod, rel = planned[uri]
        try:
            row = _fetch_one(prod, root, rel, previous.get(uri), max_size_bytes, verbose)
        except (DownloadError, OSError) as exc:
            log.error("failed %s: %s", uri, exc)
            failures.append(f"{uri}: {exc}")
            continue
        rows.append(row)

    recorded = {row["dataURI"]: row for row in _rows(read_manifest(manifest_path))}
    recorded.update({row["dataURI"]: {c: row[c] for c in schema.MANIFEST_COLUMNS} for row in rows})
    if write_manifest(recorded.values(), manifest_path):
        log.info("manifest updated: %s (%d rows)", manifest_path, len(recorded))
    else:
        log.info("manifest unchanged: %s", manifest_path)
    if failures:
        raise DownloadError(f"{len(failures)} product(s) failed:\n  " + "\n  ".join(failures))

    out = _manifest_table(rows, extra=("status",))
    out.meta.update(
        provenance=schema.Provenance.OBSERVED.value,
        source=f"{MANIFEST_SOURCE}; manifest {manifest_path.name}",
    )
    return schema.validate(out, schema.MANIFEST_COLUMNS, name="manifest")


def local_relpath(product: Row | Mapping[str, Any]) -> PurePosixPath:
    """Cache location of one product relative to the data root (POSIX separators)."""
    parts = (_cell(product, "obs_id"), _cell(product, "productFilename"))
    for part in parts:
        if not part or part in (".", "..") or any(ch in part for ch in "/\\:"):
            raise ValueError(f"unsafe path component from MAST metadata: {part!r}")
    return CACHE_PREFIX.joinpath(*parts)


def read_pipeline_version(path: str | Path, filename: str | None = None) -> str:
    """Software versions recorded inside a product file, e.g. ``jwst=2.0.1;photutils=2.3.0``.

    ``_cat.ecsv``: the ``version`` dict in the ECSV meta (jwst, photutils, astropy).
    FITS: primary header only (``CAL_VER``, ``CRDS_CTX``). Returns ``""`` when unknown.
    """
    path = Path(path)
    name = (filename or path.name).lower()
    parts: list[str] = []
    try:
        if name.endswith(".ecsv"):
            versions = Table.read(path, format="ascii.ecsv").meta.get("version") or {}
            parts = [
                f"{k}={versions[k]}" for k in ("jwst", "photutils", "astropy") if k in versions
            ]
        elif name.endswith((".fits", ".fits.gz")):
            header = fits.getheader(path, ext=0)
            parts = [
                f"{label}={header[key]}"
                for label, key in (("jwst", "CAL_VER"), ("crds", "CRDS_CTX"))
                if key in header
            ]
    except Exception as exc:  # unreadable metadata must not discard a verified download
        log.warning("could not read pipeline version from %s: %s", path, exc)
        return ""
    return ";".join(parts)


def sha256_file(path: str | Path) -> str:
    """Hex sha256 of a file, read in 1 MiB chunks."""
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(_HASH_CHUNK):
            digest.update(chunk)
    return digest.hexdigest()


def read_manifest(path: str | Path) -> Table:
    """Read a manifest; a missing file gives an empty table with ``schema.MANIFEST_COLUMNS``."""
    path = Path(path)
    if not path.is_file():
        return _manifest_table([])
    table = Table.read(path, format="ascii.ecsv")
    missing = [c for c in schema.MANIFEST_COLUMNS if c not in table.colnames]
    if missing:
        raise ValueError(f"{path}: missing manifest columns {missing}")
    out = _manifest_table(_rows(table))
    out.meta.update(table.meta)
    out.meta.setdefault("provenance", schema.Provenance.OBSERVED.value)
    out.meta.setdefault("source", MANIFEST_SOURCE)
    return out


def write_manifest(rows: Iterable[Mapping[str, Any]], path: str | Path) -> bool:
    """Atomically write manifest ``rows`` sorted by ``dataURI``; return False if unchanged."""
    table = _manifest_table(sorted(rows, key=lambda r: r["dataURI"]))
    table.meta.update(
        provenance=schema.Provenance.OBSERVED.value,
        source=MANIFEST_SOURCE,
        local_path_root="$JWST_ANOMALY_DATA (jwst_anomaly.paths.data_root)",
    )
    return write_ecsv(table, path)


def write_ecsv(table: Table, path: str | Path) -> bool:
    """Write ``table`` as ECSV with LF line endings via temp file + ``os.replace``.

    Returns False (and leaves the file untouched) when the content would not change, so
    re-running an acquisition does not churn tracked files.
    """
    path = Path(path)
    buf = io.StringIO()
    table.write(buf, format="ascii.ecsv")
    data = buf.getvalue().replace("\r\n", "\n").encode("utf-8")  # astropy uses os.linesep
    if path.is_file() and path.read_bytes() == data:
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f"{path.name}.{os.getpid()}.tmp")
    try:
        tmp.write_bytes(data)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)
    return True


def verify_manifest(
    manifest_path: str | Path,
    data_root: Path | None = None,
    *,
    uris: Iterable[str] | None = None,
) -> Table:
    """Check manifest rows (all, or only ``uris``) against the files on disk.

    Returns ``dataURI``, ``local_path`` and ``check`` (``ok``, ``missing``, ``size_mismatch``
    or ``sha256_mismatch``) per row; ``meta["provenance"] == "derived"``.
    """
    root = Path(data_root) if data_root is not None else paths.data_root()
    manifest = read_manifest(manifest_path)
    if uris is not None:
        manifest = manifest[np.isin(manifest["dataURI"], np.asarray(sorted(uris), dtype=str))]
    checks = []
    for row in manifest:
        local = root.joinpath(*PurePosixPath(row["local_path"]).parts)
        if not local.is_file():
            checks.append("missing")
        elif local.stat().st_size != int(row["size"]):
            checks.append("size_mismatch")
        elif sha256_file(local) != row["sha256"]:
            checks.append("sha256_mismatch")
        else:
            checks.append("ok")
    out = Table(
        {
            "dataURI": np.asarray(manifest["dataURI"], dtype=str),
            "local_path": np.asarray(manifest["local_path"], dtype=str),
            "check": np.asarray(checks, dtype=str),
        }
    )
    out.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source=f"sha256/size of files under {root} vs manifest {Path(manifest_path).name}",
    )
    return out


def _fetch_one(
    prod: Row,
    root: Path,
    rel: PurePosixPath,
    prev: Mapping[str, Any] | None,
    max_size_bytes: int | None,
    verbose: bool,
) -> dict[str, Any]:
    """Return the manifest row (plus ``status``) for one product, downloading if needed."""
    uri, fname = _cell(prod, "dataURI"), _cell(prod, "productFilename")
    expected = _expected_size(prod)
    dest = root.joinpath(*rel.parts)

    if dest.is_file():
        on_disk = dest.stat().st_size
        # A size the manifest already verified (download-service Content-Length + sha256) is
        # trusted over a stale MAST listing; the sha256 check below still applies.
        verified = prev is not None and int(prev["size"]) == on_disk
        if expected is not None and on_disk != expected and not verified:
            log.warning(
                "%s: %d B on disk, MAST reports %d B; re-downloading", rel, on_disk, expected
            )
        else:
            digest = sha256_file(dest)
            if prev is not None and prev["sha256"] == digest:
                _warn_if_reprocessed(prod, prev["pipeline_version"])
                return {**prev, "local_path": str(rel), "status": "cached"}
            if prev is None and expected is not None:
                mtime = datetime.fromtimestamp(dest.stat().st_mtime, UTC)
                log.info("adopting %s (already on disk, size matches MAST)", rel)
                return _new_row(prod, rel, dest, on_disk, digest, mtime, "adopted")
            if prev is not None:
                log.warning("%s: sha256 differs from the manifest; re-downloading", rel)

    if expected is None:
        raise DownloadError(
            "MAST reports no size, so neither completeness nor max_size_bytes can be checked"
        )
    if max_size_bytes is not None and expected > max_size_bytes:
        raise DownloadError(
            f"{expected} B exceeds max_size_bytes={max_size_bytes}; raise the limit (and state "
            "why) or read the file remotely via its cloud_uri"
        )
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(f"{dest.name}.{os.getpid()}.part")  # keeps a suffix for astroquery
    try:
        log.info("downloading %s (%d B)", fname, expected)
        status, msg, _url = query.mast_client().download_file(
            uri, local_path=str(tmp), cache=False, verbose=verbose
        )
        if status != "COMPLETE" or not tmp.is_file():
            hint = ""
            if "404" in str(msg):
                hint = " (MAST answers 404 for EXCLUSIVE_ACCESS products without a valid token)"
            raise DownloadError(f"{status}: {msg}{hint}")
        size = tmp.stat().st_size
        if size != expected:  # the only integrity check available: MAST publishes no checksums
            # MAST's product listing can lag a reprocessing (seen for 1176 o241 on 2026-10-01),
            # so a complete transfer of the new file is accepted when the server declares it.
            served = _served_size(uri)
            if served != size:
                raise DownloadError(
                    f"got {size} B, MAST reports {expected} B (download service: {served} B)"
                )
            log.warning(
                "%s: MAST lists %d B but serves %d B (Content-Length matches the download); "
                "the product listing is stale, probably after a reprocessing",
                fname,
                expected,
                size,
            )
        digest = sha256_file(tmp)
        os.replace(tmp, dest)
    finally:
        tmp.unlink(missing_ok=True)
    return _new_row(prod, rel, dest, size, digest, datetime.now(UTC), "downloaded")


def _served_size(uri: str) -> int | None:
    """Content-Length the MAST download service declares for ``uri`` (HTTP HEAD), or None."""
    url = f"{mast_conf.server}/api/v0.1/Download/file?uri={quote(uri, safe=':/')}"
    try:
        with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60) as r:
            length = r.headers.get("Content-Length")
        return int(length) if length is not None else None
    except (OSError, ValueError) as exc:  # URLError is an OSError
        log.warning("could not read the served size of %s: %s", uri, exc)
        return None


def _new_row(
    prod: Row,
    rel: PurePosixPath,
    dest: Path,
    size: int,
    digest: str,
    retrieved: datetime,
    status: str,
) -> dict[str, Any]:
    fname = _cell(prod, "productFilename")
    version = read_pipeline_version(dest, fname)
    _warn_if_reprocessed(prod, version)
    return {
        "dataURI": _cell(prod, "dataURI"),
        "productFilename": fname,
        "local_path": str(rel),
        "size": int(size),
        "sha256": digest,
        "retrieved_utc": retrieved.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "pipeline_version": version,
        "status": status,
    }


def _warn_if_reprocessed(prod: Row, pipeline_version: str) -> None:
    """Warn when the file's jwst version differs from MAST's current ``prvversion``.

    A reprocessed FITS product can keep its byte size, so this is the only cheap signal that a
    cached copy is stale; it is a warning, not a re-download, to avoid looping on archive
    metadata that disagrees with the file itself.
    """
    mast_version = _cell(prod, "prvversion") if "prvversion" in prod.colnames else ""
    parts = dict(p.split("=", 1) for p in pipeline_version.split(";") if "=" in p)
    file_jwst = parts.get("jwst")
    if mast_version and file_jwst and file_jwst != mast_version:
        log.warning(
            "%s: file records jwst=%s but MAST lists prvversion=%s (reprocessed? delete the "
            "local copy to re-fetch)",
            _cell(prod, "productFilename"),
            file_jwst,
            mast_version,
        )


def _manifest_table(rows: Iterable[Mapping[str, Any]], extra: Iterable[str] = ()) -> Table:
    rows = list(rows)
    cols = {}
    for name in (*schema.MANIFEST_COLUMNS, *extra):
        if name in _INT_COLUMNS:
            cols[name] = np.asarray([int(r[name]) for r in rows], dtype=np.int64)
        else:
            cols[name] = np.asarray([str(r[name]) for r in rows], dtype=str)
    return Table(cols)


def _rows(table: Table) -> list[dict[str, Any]]:
    return [{c: _cell(row, c) for c in schema.MANIFEST_COLUMNS} for row in table]


def _cell(row: Row | Mapping[str, Any], name: str) -> Any:
    """One value as a plain Python object; masked string cells become ``""``."""
    value = row[name]
    if np.ma.is_masked(value):
        return ""
    if isinstance(value, np.generic):
        value = value.item()
    return value if isinstance(value, int | float) else str(value)


def _expected_size(prod: Row) -> int | None:
    value = prod["size"]
    if np.ma.is_masked(value):
        return None
    size = int(value)
    return size if size > 0 else None
