"""List MAST HLSP provenance names around lensing-cluster centres (2026-10-25 search).

Usage: python hlsp_scan.py out.json  (cloud sessions: SSL_CERT_FILE = proxy bundle + certifi).
The output is the `provenance_names_by_field` block of search.json; the product checks there were
made by hand with `Observations.get_product_list`. Exits 1 if any field query failed.
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor

from astroquery.mast import Observations

FIELDS = {
    "SMACS0723": "110.8375 -73.4544",
    "A2744": "3.5875 -30.3972",
    "ElGordo": "15.7208 -49.2494",
    "WHL0137": "24.3554 -8.4572",
    "MACS0416": "64.0381 -24.0675",
    "MACS1149": "177.3988 22.3985",
    "A370": "39.9713 -1.5823",
    "MACS0647": "101.9623 70.2475",
    "MACS0417": "64.3946 -11.9089",
    "MACS1423": "215.9495 24.0784",
    "Bullet": "104.6292 -55.9469",
    "RXJ2129": "322.4165 0.0892",
    "MACS1931": "292.9568 -26.5758",
    "AS1063": "342.1832 -44.5309",
    "MACS0717": "109.3982 37.7452",
    "RXJ0152": "28.1824 -13.9625",
    "MACS0553": "88.3304 -33.7076",
    "SPT0615": "93.9656 -57.7801",
}


def query(item):
    name, coords = item
    try:
        t = Observations.query_criteria(
            coordinates=coords, radius="0.05 deg", obs_collection="HLSP"
        )
        return name, sorted({str(p) for p in t["provenance_name"]})
    except Exception as exc:  # recorded, not raised: one field failing should not stop the scan
        return name, ["ERR " + repr(exc)[:120]]


if __name__ == "__main__":
    with ThreadPoolExecutor(6) as ex:
        out = dict(ex.map(query, FIELDS.items()))
    with open(sys.argv[1], "w") as fh:
        json.dump(out, fh, indent=1)
    failed = [k for k, v in out.items() if any(p.startswith("ERR ") for p in v)]
    print(sorted({p for v in out.values() for p in v if not p.startswith("ERR ")}))
    if failed:
        sys.exit(f"query failed for: {failed}")
