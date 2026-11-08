"""Compare two epoch-inventory snapshots: epochs added, removed or newly public per field.

Works on the ``inventory.json`` of ``epoch_inventory.py`` and ``epoch_inventory_ext.py``.
Field-centre shifts above 0.5'' are reported when both snapshots store ``centre_deg``.

    python scripts/epoch_inventory_diff.py OLD.json NEW.json --out diff.json
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def diff_inventories(old: dict, new: dict, centre_tol_arcsec: float = 0.5) -> dict:
    out: dict[str, dict] = {}
    for field, v in new["fields"].items():
        o = old["fields"].get(field, {})
        oe, ne = o.get("epochs", {}), v.get("epochs", {})
        rec = {
            "added": sorted(set(ne) - set(oe)),
            "removed": sorted(set(oe) - set(ne)),
            "newly_public": sorted(
                k for k, e in ne.items() if e["public"] and not oe.get(k, {}).get("public", False)
            ),
        }
        if "centre_deg" in v and "centre_deg" in o:
            (ra0, de0), (ra1, de1) = o["centre_deg"], v["centre_deg"]
            dra = (ra0 - ra1 + 180.0) % 360.0 - 180.0
            d = 3600 * math.hypot(dra * math.cos(math.radians(de0)), de0 - de1)
            if d > centre_tol_arcsec:
                rec["centre_shift_arcsec"] = round(d, 2)
        out[field] = rec
    for field in sorted(set(old["fields"]) - set(new["fields"])):
        out[field] = {"missing_in_new": True}
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("old", type=Path)
    ap.add_argument("new", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    old, new = (json.loads(p.read_text()) for p in (a.old, a.new))
    d = diff_inventories(old, new)
    meta = {"old": str(a.old), "new": str(a.new), "old_meta": old["meta"], "new_meta": new["meta"]}
    a.out.write_text(json.dumps({"meta": meta, "fields": d}, indent=1) + "\n")
    changed = {
        f: r for f, r in d.items() if any(r.get(k) for k in ("added", "removed", "newly_public"))
    }
    print(f"{len(d)} fields, {len(changed)} with epoch changes: {sorted(changed)}")


if __name__ == "__main__":
    main()
