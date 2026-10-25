"""Offline checks of the SLICE x RELICS held-out inputs (PLAN.md)."""

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import slice_relics_heldout as sr  # noqa: E402

_TEX = r"""
\sidehead {{X\,J0000}}\hline\\[-5pt]
1.1.1 & 10.0 & -5.0 & 1.645 & \cite{a} & \\
1.1.2 & 10.1 & -5.1 &  &  &  \\
1.2.1 & 10.2 & -5.2 & \purple{ $3.04_{-0.04}^{+0.04} $} & & \\
1.2.2 & 10.3 & -5.3 &  & & \\
\hline
2.1.1 & 11.0 & -6.0 & \purple{ $5.76_{-0.86}^{+0.15}$ } &  &  \\
2.1.2 &11.1 & -6.1 &  &  &  \\
  \hline\hline\\
\sidehead {{Y}}\hline\\[-5pt]
9.1.1 & 1.0 & 1.0 & 2.0 & x & \\
"""


def test_parse_block_keeps_clump1_and_labels_z():
    t = sr.parse_block(_TEX, r"X\,J0000")
    assert list(t["image_id"]) == ["1.1", "1.2", "2.1", "2.2"]
    assert list(t["z"]) == [1.645, 1.645, 5.76, 5.76]
    assert list(t["z_spec"]) == [True, True, False, False]
    assert list(t["paper_id"]) == ["1.1.1", "1.1.2", "2.1.1", "2.1.2"]


def test_parse_block_raises_on_an_unparsed_row():
    bad = _TEX.replace("2.1.2 &11.1 & -6.1", "2.1.2 & 11.1 & $-$6.1")
    with pytest.raises(ValueError, match="unparsed Table 8 row"):
        sr.parse_block(bad, r"X\,J0000")


def test_stored_image_list_matches_the_frozen_families():
    t = Table.read(sr.IMAGES, format="ascii.ecsv")
    # held-out and spectroscopic families as read from the SLICE text (PLAN.md table)
    expect = {
        "rxc0232": (set(), {"5", "7"}),
        "rxc2211": ({str(i) for i in range(4, 11)}, set()),
        "macs0553": ({str(i) for i in range(11, 22)}, {"1", "5", "10", "13"}),
    }
    for cluster, (held, spec) in expect.items():
        c = sr.load_cluster(t, cluster)
        fam = np.asarray(c["system"], str)
        assert set(fam[np.asarray(c["held_out"], bool)]) == set(held)
        assert set(fam[np.asarray(c["z_spec"], bool)]) == set(spec)
    n = {c: int(np.sum(np.asarray(t["cluster"], str) == c)) for c in sr.CLUSTERS}
    assert n == {
        "rxc0232": 26,
        "rxc2211": 33,
        "macs0553": 57,
    }  # text: 26 / 32 / 57 images (RXC2211: Table 8 has 33)
