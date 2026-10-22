"""Unit tests of the MACS1149 containing-source rule (plan rule 4) and the gold-list copy."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import m0417_arcref as ar  # noqa: E402
import m1149_canucsref as m  # noqa: E402

DEC0 = 22.4


def _cat(rows):
    """rows: (dx East, dy North, a, b, pa) in arcsec / deg, around (180, DEC0)."""
    r = np.asarray(rows, float)
    return {
        "ra": 180.0 + r[:, 0] / 3600.0 / np.cos(np.radians(DEC0)),
        "dec": DEC0 + r[:, 1] / 3600.0,
        "a": r[:, 2],
        "b": r[:, 3],
        "pa": r[:, 4],
    }


def test_containing_finds_long_arc_with_far_centroid():
    # 0: arc centred 1.5" north, major axis N-S, 3a = 3"; 1: compact source 0.5" east;
    # 2: arc 1.5" east with major axis N-S (does not contain the origin across its minor axis)
    cat = _cat([(0, 1.5, 1.0, 0.1, 0.0), (0.5, 0, 0.05, 0.05, 0.0), (1.5, 0, 1.0, 0.1, 0.0)])
    assert ar.containing(180.0, DEC0, cat).tolist() == [0]


def test_compared_sources_union_and_reference_excluded():
    cat = _cat([(0, 1.5, 1.0, 0.1, 0.0), (0.5, 0, 0.05, 0.05, 0.0), (1.5, 0, 1.0, 0.1, 0.0)])
    assert ar.compared_sources(180.0, DEC0, 1.0, cat, -1) == [1]
    assert ar.compared_sources(180.0, DEC0, 1.0, cat, -1, contain=True) == [0, 1]
    assert ar.compared_sources(180.0, DEC0, 1.0, cat, 0, contain=True) == [1]
    assert ar.compared_sources(180.0, DEC0, 0.2, cat, -1, contain=True) == [0]


def test_written_images_are_spectroscopic(tmp_path):
    t = m.write_images(tmp_path)
    from astropy.table import Table

    tab = Table.read(t, format="ascii.ecsv")
    assert len(tab) == 23 and bool(np.all(tab["z_spec"]))
