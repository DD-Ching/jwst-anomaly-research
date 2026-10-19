"""Fetch the committed pipeline catalogues from S3 (cloud_uri in the products manifest) to the
data root, checking sha256 against the fetch manifest. Run from the repo root: fetch_cats.py FIELD...
"""

import hashlib
import sys
from concurrent.futures import ThreadPoolExecutor

import fsspec
from astropy.table import Table

from jwst_anomaly import paths


def one(row, sha):
    dest = paths.data_root() / str(row["local_path"])
    if dest.exists() and hashlib.sha256(dest.read_bytes()).hexdigest() == sha:
        return dest.name, "cached"
    dest.parent.mkdir(parents=True, exist_ok=True)
    with fsspec.open(str(row["cloud_uri"]), "rb", anon=True) as f:
        data = f.read()
    got = hashlib.sha256(data).hexdigest()
    if got != sha:
        raise ValueError(f"{dest.name}: sha256 {got} != manifest {sha}")
    dest.write_bytes(data)
    return dest.name, len(data)


jobs = []
for field in sys.argv[1:]:
    man = Table.read(f"data/manifests/dimming_{field}.ecsv")
    prods = Table.read(f"data/manifests/dimming_{field}_products.ecsv")
    uri = {str(r["productFilename"]): str(r["cloud_uri"]) for r in prods}
    for r in man:
        jobs.append(
            (
                {"local_path": r["local_path"], "cloud_uri": uri[str(r["productFilename"])]},
                str(r["sha256"]),
            )
        )
with ThreadPoolExecutor(8) as ex:
    for name, n in ex.map(lambda j: one(*j), jobs):
        print(name, n, flush=True)
