"""Offline checks of the MACS0647 image-list parsers and the committed tables (input check)."""

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import macs0647_inputs as mi  # noqa: E402

from jwst_anomaly import lensmodel  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / mi.OUT

MEENA = r"""
 1.1 & 101.9660445 & 70.2558166 & $2.13_{-0.18}^{+0.33}$ & \citet{x} \\
 1.3 & 101.9666816 & 70.2483074 & -- & -- \vspace{0.02in} \\
% 2.2a & 102.0014192 & 70.2496564 & $4.81_{-0.19}^{+0.15}$ & New pair \\
 6.12 & 101.9713212 & 70.2397022 & $10.58_{-0.37}^{+0.32}$ & -- \\
 8.3? & 101.9433492 & 70.2385144 & $2.37_{-0.30}^{+0.23}$ & -- \vspace{0.02in} \\
"""

NISHIDA = r"""
 1.1 & 101.966807 &  70.248308 & $\checkmark$ &   $ 2.134$ & $\checkmark$ & \nodata & 1,2\\
 1.2 & 101.955537 &  70.249205 & $\checkmark$ &            & $\checkmark$ & & \\
 3.1 & 101.974507 &  70.243384 & $\checkmark$ & $3.0\pm0.5$ & $\checkmark$ & \nodata & 1,2\\
 3.2 & 101.972537 &  70.242639 &              &            & $\checkmark$ & & \\
 4.1 & 101.9 &  70.2 & $\checkmark$ &  \nodata   & & & 2 \\
%20.3 & 101.974230 &  70.237479 & $\checkmark$ &            & & & %\\
"""


def test_parse_meena_rows():
    t = mi.parse_lensed_systems(MEENA)
    assert list(t["paper_id"]) == ["1.1", "1.3", "6.12", "8.3?"]  # commented row left out
    assert list(t["image_id"]) == ["1.1", "1.3", "6-1.2", "8.3"]
    assert list(t["system"]) == ["1", "1", "6-1", "8"]
    assert list(t["uncertain"]) == [False, False, False, True]
    assert t["z_phot"][0] == 2.13 and np.isnan(t["z_phot"][1]) and t["z_phot"][2] == 10.58
    assert not any(t["z_spec"]) and not any(t["new_in_jwst"])


def test_parse_nishida_rows():
    t = mi.parse_nishida25(NISHIDA)
    assert list(t["paper_id"]) == ["1.1", "1.2", "3.1", "4.1"]  # no JWST mark: 3.2 left out
    assert list(t["z_jwst"][:3]) == [2.134, 2.134, 3.0]  # carried from the family's first row
    assert np.isnan(t["z_jwst_pm"][0]) and t["z_jwst_pm"][2] == 0.5
    assert np.isnan(t["z_jwst"][3]) and np.isnan(t["z_jwst_pm"][3])


def test_committed_tables_match_the_papers():
    """Counts the papers state: 23 systems (Meena+2023, JD's three components as families) and
    86 images from 28 sources (Nishida+2025); positions pinned by the summary's hashes."""
    m = Table.read(OUT / mi.IMAGES)
    n = Table.read(OUT / mi.IMAGES_N25)
    assert (len(m), len(set(m["system"]))) == (72, 25)
    assert (len(n), len(set(n["system"]))) == (86, 28)
    assert int(np.isnan(m["z_phot"]).sum()) == 13  # the rows printed "--"
    assert list(dict.fromkeys(m["system"][m["new_in_jwst"]])) == [str(i) for i in range(13, 24)]
    for t in (m, n):
        assert [lensmodel.image_family(i) for i in t["image_id"]] == list(t["system"])
        assert t.meta["provenance"] == "observed"
        assert np.all((t["dec"] > 70.23) & (t["dec"] < 70.27))
    import json

    s = json.loads((OUT / "summary.json").read_text())
    assert mi.table_sha256(m) == s["images"]["table_sha256"]
    assert mi.table_sha256(n) == s["images_nishida25"]["table_sha256"]


def test_parsers_refuse_unparsed_rows():
    with pytest.raises(ValueError):
        mi.parse_lensed_systems(r" 2.2a & 102.0014192 & 70.2496564 & -- & -- \\")
    with pytest.raises(ValueError):
        mi.parse_lensed_systems(r" 2.2 & 102.0014192 & 70.2496564 & $>6$ & -- \\")
    with pytest.raises(ValueError):
        mi.parse_nishida25(r" 5.1 & 101.9 & 70.2 & $\checkmark$ & $3.0^{+0.4}$ & & & \\")


def test_nishida_z_from_a_first_row_without_jwst_mark():
    tex = (
        " 7.1 & 101.9 & 70.21 & & $2.5$ & $\\checkmark$ & & 1\\\\\n"
        " 7.2 & 101.9 & 70.22 & $\\checkmark$ & & $\\checkmark$ & & \\\\\n"
    )
    t = mi.parse_nishida25(tex)
    assert list(t["paper_id"]) == ["7.2"] and t["z_jwst"][0] == 2.5
