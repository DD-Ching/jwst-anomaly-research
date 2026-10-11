"""Publish reduced derived products as GitHub Release assets and pin them (owner brief 2026-10-10,
S-1).

Run from a session that may create releases (the owner's machine or a local session). Cloud routine
sessions get HTTP 403 ("Creating, editing, or deleting releases is not permitted for this session
type").

  python scripts/derived_publish.py data/e1_events/gw_skymaps_nside32.npz \
      --sources "Zenodo 6513631, 8177023, 20275769, 20348005 (GWTC PE sky maps, CC BY 4.0)"

It runs ``gh release create derived-data-YYYYMMDD <files> --target <HEAD> --prerelease
--latest=false`` and, only when that succeeds, appends the manifest rows
(``data/manifests/derived_data.ecsv``). The tag is created at the recorded code commit, so it
keeps that commit reachable after a squash merge; push it first (a commit GitHub does not have is
refused and nothing is pinned). Commit the manifest change in a PR afterwards. Only reduced
products below ``derived_store.MAX_BYTES``; never raw archives.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from jwst_anomaly import derived_store as ds  # noqa: E402


def main(argv=None) -> int:
    from astropy.table import Table, vstack

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("files", nargs="+", type=Path)
    ap.add_argument("--sources", required=True, help="where the inputs came from (URLs, DOIs)")
    ap.add_argument("--tag", default=time.strftime("derived-data-%Y%m%d", time.gmtime()))
    ap.add_argument("--dry-run", action="store_true", help="write nothing, run nothing")
    a = ap.parse_args(argv)
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True, cwd=ROOT
    ).stdout.strip()
    created = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    rows = [ds.manifest_row(f, a.tag, a.sources, commit, created) for f in a.files]
    new = Table(rows=rows)
    for r in rows:
        print(f"{r['name']}: {r['size']} B sha256 {r['sha256']} -> {a.tag}")
    if a.dry_run:
        return 0
    notes = (
        "Derived (reduced) data products for ephemeral cloud sessions (S-1). Not a code release. "
        f"Inputs: {a.sources}. Pinned in data/manifests/derived_data.ecsv; code commit {commit}."
    )
    cmd = ["gh", "release", "create", a.tag, *map(str, a.files), "--target", commit]
    cmd += ["--prerelease", "--latest=false"]
    subprocess.run([*cmd, "--title", a.tag, "--notes", notes], check=True, cwd=ROOT)
    # pin only after the release exists (a refused or failed release leaves no dangling row)
    if ds.MANIFEST.exists():
        new = vstack([Table.read(ds.MANIFEST, format="ascii.ecsv"), new])
    new.meta = {
        "provenance": "derived",
        "source": "reduced products published as GitHub Release assets (S-1, derived_store.py)",
    }
    new.write(ds.MANIFEST, format="ascii.ecsv", overwrite=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
