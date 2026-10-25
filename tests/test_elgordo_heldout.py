import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import elgordo_heldout as eh  # noqa: E402

TEX = r"""
%    ID   &     RA    &    DEC     &  Zs    &   Zg   & Zphot  \\
    1.1a  & 15.722308 & -49.254551 & 2.5635 &   2.64\ddag &  3.46& 2.35 &  201.2 & A  \\
    1.1c$^{\dag}$  & 15.730896 & -49.250114 &  \no   &  \no   &  3.30& 2.35 &  292.3 & A  \\
    43a & 15.728325 & -49.266068 &  \no   &   3.75 &  -1  & 3.33   &  0.9  & B  \\
    43b & 15.725151 & -49.267990 &  \no   &  \no   &  -1  & -1     &  -1   & B  \\
\hline
"""


def test_parse_table_family_z():
    t = eh.parse_table(TEX)
    assert list(t["image_id"]) == ["1.1a", "1.1c", "43a", "43b"]
    assert np.allclose(t["z"], [2.5635, 2.5635, 3.75, 3.75])
    assert list(t["z_spec"]) == [True, True, False, False]


def test_primary_gate_uses_spec_rows_only():
    images = Table({"image_id": [f"{k}a" for k in range(1, 8)], "z_spec": [True] * 5 + [False] * 2})
    r = Table(
        {
            "image_id": images["image_id"],
            "dtheta_arcsec": [0.2, 0.3, 0.4, 0.5, 0.6, 9.0, 9.0],
            "inside": [True] * 7,
        }
    )
    p = eh.primary_gate(r, images)
    assert p["n_rows"] == 5
    assert p["outcome"] == "reproduces"
