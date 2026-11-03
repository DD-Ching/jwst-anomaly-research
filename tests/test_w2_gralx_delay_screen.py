"""Offline tests for scripts/w2_gralx_delay_screen.py (parsing, matching, stored summary)."""

import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "w2gralx", ROOT / "scripts" / "w2_gralx_delay_screen.py"
)
w = importlib.util.module_from_spec(spec)
sys.modules["w2gralx"] = w
spec.loader.exec_module(w)
OUT = ROOT / "results" / "w2_gralx_delay_screen"

TEX = r"""
\caption{\label{timedelay_table}Time delays of known gravitational lens systems. } \\
Name & $\Delta$$t_{AB}$ & $\Delta$$t_{AC}$ & $\Delta$$t_{AD}$ & x & x & x & Author \\
\multicolumn{8}{c}{[days]}\\
  HE0047-1756 & $-10.8 \pm 1.0$ &  &  &  &  &  & 2020Millon\_b\\
%  B0218+357 & $11.3 \pm 0.2$ &  &  &  &  &  & 2018Biggs\\
%  APM08279+5255 & $900.0$ &  &  &  &  &  & 2016Saturni\\
%  VDESJ2325-5229 & $43.8 \pm 4.5$ &  &  &  &  &  & 2020Millon\_b\\\hline
\end{longtable}
"""


def test_parse_rows_reads_printed_and_commented_rows():
    rows = w.parse_rows(TEX)
    assert [r["name"] for r in rows] == [
        "HE0047-1756",
        "B0218+357",
        "APM08279+5255",
        "VDESJ2325-5229",
    ]
    assert rows[0]["delays"]["AB"] == {"dt": -10.8, "err": 1.0}
    assert rows[0]["ref"] == "2020Millon_b"
    assert rows[2]["delays"]["AB"] == {"dt": 900.0, "err": None}
    assert rows[3]["ref"] == "2020Millon_b"


def test_name_coords_and_designation_match():
    ra, dec = w.name_coords("DESJ0029-3814")
    assert ra == pytest.approx(15 * (0 + 29 / 60)) and dec == pytest.approx(-(38 + 14 / 60))
    assert w.name_coords("PG1115+080") is None  # B1950: resolved through SIMBAD instead
    # same HHMM±DDMM designation (±1 in the last digit) counts as one system
    assert w._desig(169.5706, 7.7662) == w._desig(169.5704, 7.7660)
    a, b = w._desig(141.0, 2.3167), w._desig(141.2325, 2.3236)
    assert a[1] == b[1] and abs(a[0] - b[0]) == 0


@pytest.mark.skipif(not (OUT / "summary.json").exists(), reason="no stored summary")
def test_stored_summary_is_consistent():
    s = json.loads((OUT / "summary.json").read_text())
    t2 = json.loads((OUT / "table2.json").read_text())
    d = json.loads((OUT / "deflectors.json").read_text())["systems"]
    assert s["n_systems"] == t2["n_systems"] == len(t2["rows"]) == 73
    assert sum(s["counts"].values()) == 73 and s["not_curated"] == []
    assert set(d) == {r["name"] for r in t2["rows"]}
    assert all(v["category"] in w.CATEGORIES for v in d.values())
    # a system qualifies only with a published limit, a reliable delay and high-resolution imaging
    want = sorted(
        n
        for n, v in d.items()
        if v["category"] == "quantitative_limit"
        and v.get("delay_reliable") is True
        and v.get("imaging_high_res") is True
    )
    assert sorted(s["qualifying"]) == want
