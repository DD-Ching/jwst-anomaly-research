"""F200W / F444W forced-photometry cutouts of one field's targets, all epochs (blank: no
coverage). Usage (repo root): sheet.py FIELD OUT PNG
"""

import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, "scripts")
import dimming_screen as ds  # noqa: E402

name, out, png = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
field = ds.load_field(ds.DEFAULT_CONFIG, name)
t = Table.read(out / "forced.ecsv")
rows = np.flatnonzero(~np.asarray(t["control"], bool))
files = {}
for b in field["bands"]:
    for k in range(len(field["epochs"])):
        d = next((out / "cutouts" / f"{b}_e{k}").glob("*"), None)
        files[(k, b)] = [
            str(d / f"{u}.fits") if d and (d / f"{u}.fits").exists() else "" for u in t["uid"]
        ]
ds._epoch_sheet(t, rows, files, field, png)
