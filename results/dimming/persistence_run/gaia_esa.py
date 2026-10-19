import sys
import urllib.parse
import urllib.request

import numpy as np
from astropy.table import Table

sys.path.insert(0, "scripts")
import dimming_screen as ds

root = sys.argv[1]  # output root shared with whatif.py
for name in sys.argv[2:]:
    ra, dec = ds.load_field(ds.DEFAULT_CONFIG, name)["centre"]
    q = (
        f"SELECT source_id, ra, dec, phot_g_mean_mag FROM gaiadr3.gaia_source WHERE "
        f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},0.1))"
    )
    url = "https://gea.esac.esa.int/tap-server/tap/sync?" + urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q}
    )
    t = Table.read(urllib.request.urlopen(url, timeout=120).read().decode(), format="ascii.csv")
    out = Table(
        {
            "source_id": np.asarray(t["source_id"], np.int64),
            "ra": np.asarray(t["ra"], float),
            "dec": np.asarray(t["dec"], float),
            "gmag": np.ma.filled(np.ma.asarray(t["phot_g_mean_mag"], float), np.nan),
        }
    )
    out.meta = {
        "provenance": "observed",
        "source": "ESA Gaia archive TAP gaiadr3.gaia_source (J2016.0); VizieR timed out",
    }
    p = f"{root}/{name}/gaia.ecsv"
    import os

    os.makedirs(os.path.dirname(p), exist_ok=True)
    out.write(p, overwrite=True)
    print(name, len(out))
