"""Durable store for small derived products across ephemeral cloud sessions (owner brief
2026-10-10, S-1).

Reduced products (< 50 MB each, never raw archives) are published as assets of a GitHub Release
tagged ``derived-data-YYYYMMDD``. Each asset is pinned in ``data/manifests/derived_data.ecsv``
(asset name, tag, sha256, size, source URLs, code commit, creation date). Scripts call
:func:`fetch` first and fall back to the original source when the durable copy is missing or does
not match its pin.
"""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

from .paths import repo_root

REPO = "DD-Ching/jwst-anomaly-research"
MANIFEST = repo_root() / "data" / "manifests" / "derived_data.ecsv"
#: ASSUMPTION: assets above this size do not belong in the store (brief: < 50 MB).
MAX_BYTES = 50_000_000


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def pinned(name: str, manifest: Path = MANIFEST) -> dict | None:
    """The newest manifest row for asset ``name`` (as a dict), or None."""
    if not manifest.exists():
        return None
    from astropy.table import Table

    t = Table.read(manifest, format="ascii.ecsv")
    rows = [r for r in t if r["name"] == name]
    if not rows:
        return None
    r = max(rows, key=lambda r: str(r["tag"]))
    return {c: r[c].item() if hasattr(r[c], "item") else r[c] for c in t.colnames}


def asset_url(row: dict) -> str:
    return f"https://github.com/{REPO}/releases/download/{row['tag']}/{row['name']}"


def fetch(name: str, dest: Path, manifest: Path = MANIFEST, timeout: float = 120.0) -> Path | None:
    """Download the pinned durable copy of ``name`` to ``dest`` and verify its sha256.

    Returns ``dest`` when ``dest`` already matches the pin or the download matches; None when
    there is no pin, the host is unreachable, or the bytes do not match (nothing is written)."""
    row = pinned(name, manifest)
    if row is None:
        return None
    dest = Path(dest)
    if dest.exists() and sha256_file(dest) == row["sha256"]:
        return dest
    import requests

    try:
        r = requests.get(asset_url(row), timeout=timeout, allow_redirects=True)
    except requests.RequestException:
        return None
    if r.status_code != 200 or hashlib.sha256(r.content).hexdigest() != row["sha256"]:
        return None
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    tmp.write_bytes(r.content)
    os.replace(tmp, dest)
    return dest


def manifest_row(path: Path, tag: str, sources: str, commit: str, created: str) -> dict:
    """A manifest row for a file about to be published (refuses files above ``MAX_BYTES``)."""
    size = Path(path).stat().st_size
    if size > MAX_BYTES:
        raise ValueError(
            f"{path}: {size} bytes > {MAX_BYTES}; the durable store holds reduced products only"
        )
    return {
        "name": Path(path).name,
        "tag": tag,
        "sha256": sha256_file(path),
        "size": size,
        "sources": sources,
        "code_commit": commit,
        "created": created,
    }
