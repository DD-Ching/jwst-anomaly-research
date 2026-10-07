"""Fetch the products of a reference-sample config and record tracked manifests.

Writes under ``data/manifests/`` (both rewritten only when their content changes):

* ``<config stem>.ecsv`` -- one row per downloaded file: dataURI, local path, size, sha256,
  retrieval time and pipeline version (``schema.MANIFEST_COLUMNS``);
* ``<config stem>_products.ecsv`` -- every level-3 product MAST lists for the sample's
  observations, with ``cloud_uri`` for S3 byte-range reads, including products that were
  not downloaded (e.g. NIRCam ``_i2d.fits`` above ``--max-size-mb``).

Examples (``CFG=configs/reference_sample.yaml``)::

    python scripts/fetch_reference_sample.py --config $CFG --catalogs-only
    python scripts/fetch_reference_sample.py --config $CFG --sample smacs0723_miri --dry-run
    python scripts/fetch_reference_sample.py --config $CFG --verify-only
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import yaml
from astropy.table import Table, vstack

from jwst_anomaly import acquire, paths, query, schema

INDEX_COLUMNS = (
    "dataURI",
    "sample_id",
    "obs_id",
    "obsID",
    "productFilename",
    "productSubGroupDescription",
    "calib_level",
    "size",
    "dataRights",
    "prvversion",
    "cloud_uri",
)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    if hasattr(sys.stdout, "reconfigure"):  # keep prints in order with log lines when piped
        sys.stdout.reconfigure(line_buffering=True)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    config_path = Path(args.config)
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    manifest_dir = Path(args.manifest_dir) if args.manifest_dir else paths.manifests_dir()
    manifest_path = manifest_dir / f"{config_path.stem}.ecsv"
    index_path = manifest_dir / f"{config_path.stem}_products.ecsv"
    data_root = Path(args.data_root) if args.data_root else paths.data_root()

    if args.verify_only:
        return _verify(manifest_path, data_root, uris=None)

    archive = config["archive"]
    if archive.get("data_rights", query.PUBLIC) != query.PUBLIC:
        print("error: this script fetches PUBLIC data only", file=sys.stderr)
        return 2
    sample_of = _requested_obs_ids(config["samples"], args.sample)
    if sample_of is None:
        return 2

    obs = query.query_observations(
        obs_collection=archive.get("collection", query.COLLECTION),
        obs_id=sorted(sample_of),
        calib_level=int(archive["calib_level"]),
    )
    obs_ids = query.str_values(obs["obs_id"])
    obs = obs[np.isin(obs_ids, sorted(sample_of))]
    missing = sorted(set(sample_of) - set(obs_ids))
    if missing:
        print(
            f"error: not found on MAST as PUBLIC level-3 observations: {missing}", file=sys.stderr
        )
        return 2
    print(f"{len(obs)} observations in {len(set(sample_of.values()))} sample(s)")

    products = query.list_products(
        obs,
        subgroups=archive["product_subgroups"],
        calib_level=int(archive["calib_level"]),
        cloud_uris=True,
    )
    # list_products guarantees each product's obsID is one of ``obs``
    sample_by_obsid = {
        obsid: sample_of[obs_id]
        for obsid, obs_id in zip(
            query.str_values(obs["obsid"]), query.str_values(obs["obs_id"]), strict=True
        )
    }
    products["sample_id"] = [sample_by_obsid[o] for o in query.str_values(products["obsID"])]
    wanted = ["CAT"] if args.catalogs_only else list(archive["product_subgroups"])
    max_bytes = int(args.max_size_mb * 1e6)
    sizes = np.ma.filled(np.ma.asarray(products["size"]), -1)  # unknown: fetch_products refuses
    in_scope = np.isin(query.str_values(products["productSubGroupDescription"]), wanted)
    too_big = in_scope & (sizes > max_bytes)
    selected = products[in_scope & ~too_big]
    print(
        f"{len(products)} level-3 products listed; {len(selected)} selected "
        f"({np.sum(selected['size']) / 1e6:.1f} MB), {int(np.sum(too_big))} skipped as larger "
        f"than {args.max_size_mb:g} MB (read those remotely via cloud_uri)"
    )
    for row in products[too_big]:
        print(f"  skipped {row['productFilename']} ({row['size'] / 1e6:.1f} MB) {row['cloud_uri']}")

    if args.dry_run:
        for row in selected:
            print(f"  would fetch {row['productFilename']} ({row['size'] / 1e6:.2f} MB)")
        return 0

    source = (
        f"MAST level-3 product lists for the observations in {config_path.name} "
        f"(calib_level={archive['calib_level']}, subgroups={list(archive['product_subgroups'])}, "
        "PUBLIC; level-2 members excluded); cloud_uri from MAST path_lookup"
    )
    if _write_index(products, set(query.str_values(obs["obs_id"])), source, index_path):
        print(f"product index updated: {index_path}")
    if not len(selected):
        print("nothing to fetch")
        return 0
    try:
        rows = acquire.fetch_products(
            selected, data_root=data_root, manifest_path=manifest_path, max_size_bytes=max_bytes
        )
    except acquire.DownloadError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    for row in rows:
        print(
            f"  {row['status']:<10} {row['size']:>12,d} B  {row['sha256'][:12]}  "
            f"{row['pipeline_version'] or '-':<40} {row['productFilename']}"
        )
    counts = {s: int(np.sum(rows["status"] == s)) for s in ("downloaded", "adopted", "cached")}
    print("summary: " + ", ".join(f"{k} {v}" for k, v in counts.items()))
    return _verify(manifest_path, data_root, uris=set(np.asarray(rows["dataURI"], dtype=str)))


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--config", required=True, help="reference-sample YAML config")
    parser.add_argument(
        "--catalogs-only", action="store_true", help="fetch only level-3 source catalogs (CAT)"
    )
    parser.add_argument(
        "--sample", action="append", metavar="ID", help="restrict to this sample id (repeatable)"
    )
    parser.add_argument(
        "--max-size-mb",
        type=float,
        default=acquire.MAX_DOWNLOAD_BYTES / 1e6,
        help="skip products larger than this (default %(default)g; CLAUDE.md asks for a reason "
        "to download more in one file)",
    )
    parser.add_argument("--dry-run", action="store_true", help="list, do not download or write")
    parser.add_argument(
        "--verify-only", action="store_true", help="check the manifest against files on disk"
    )
    parser.add_argument("--data-root", help="override $JWST_ANOMALY_DATA")
    parser.add_argument("--manifest-dir", help="override data/manifests (for testing)")
    return parser.parse_args(argv)


def _requested_obs_ids(samples: list[dict], only: list[str] | None) -> dict[str, str] | None:
    """Map obs_id -> sample id for the selected samples (None after printing an error)."""
    known = [s["id"] for s in samples]
    unknown = sorted(set(only or ()) - set(known))
    if unknown:
        print(f"error: unknown sample id(s) {unknown}; known: {known}", file=sys.stderr)
        return None
    sample_of: dict[str, str] = {}
    for sample in samples:
        if not only or sample["id"] in only:
            for obs_id in sample["obs_ids"]:
                sample_of[str(obs_id)] = sample["id"]
    return sample_of


def _write_index(products: Table, listed_obs_ids: set[str], source: str, path: Path) -> bool:
    """Update the product index at ``path``; False if its content did not change.

    Rows of the observations just listed are replaced (so products MAST no longer lists drop
    out); rows of other observations, e.g. from an earlier ``--sample`` run, are kept.
    """
    index = _index_columns(products)
    if path.is_file():
        old = _index_columns(Table.read(path, format="ascii.ecsv"))
        stale = np.isin(old["obs_id"], sorted(listed_obs_ids)) | np.isin(
            old["dataURI"], index["dataURI"]
        )
        index = vstack([old[~stale], index])
    index.sort("dataURI")
    index.meta.update(provenance=schema.Provenance.OBSERVED.value, source=source)
    schema.validate(index, schema.PRODUCT_COLUMNS, name="product index")
    return acquire.write_ecsv(index, path)


def _index_columns(table: Table) -> Table:
    """``INDEX_COLUMNS`` of ``table`` as plain columns (masked/absent: ``""``, or -1 for ints)."""
    cols = {}
    for name in INDEX_COLUMNS:
        if name in ("calib_level", "size"):
            col = table[name] if name in table.colnames else np.full(len(table), -1)
            cols[name] = np.ma.filled(np.ma.asarray(col), -1).astype(np.int64)
        else:
            col = table[name] if name in table.colnames else np.full(len(table), "")
            cols[name] = query.str_values(col)
    return Table(cols)


def _verify(manifest_path: Path, data_root: Path, uris: set[str] | None) -> int:
    if not manifest_path.is_file():
        print(f"error: no manifest at {manifest_path}", file=sys.stderr)
        return 1
    checks = acquire.verify_manifest(manifest_path, data_root, uris=uris)
    if not len(checks):
        print(f"error: nothing to verify in {manifest_path.name}", file=sys.stderr)
        return 1
    bad = checks[checks["check"] != "ok"]
    print(f"verified {len(checks) - len(bad)}/{len(checks)} files against {manifest_path.name}")
    for row in bad:
        print(f"  {row['check']}: {row['local_path']}", file=sys.stderr)
    return 1 if len(bad) else 0


if __name__ == "__main__":
    sys.exit(main())
