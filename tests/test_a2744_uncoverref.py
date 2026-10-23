"""Unit tests of the Abell 2744 UNCOVER image-list copy (PLAN.md Inputs)."""

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

from jwst_anomaly import lensmodel

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import a2744_uncover as au  # noqa: E402
from a2744_uncoverref import image_list, lenstool_id  # noqa: E402


def test_lenstool_id_keeps_the_knot_level_source():
    for uid in ("1.11", "1.13", "78.21", "83.1", "83.3", "12.2"):
        assert lensmodel.image_family(lenstool_id(uid)) == au.uncover_system(uid)
    assert lenstool_id("78.22") == "78.2b"
    assert lenstool_id("83.3") == "83c"


def test_image_list_flags_spec_and_drops_none():
    raw = Table(
        {
            "ID": [1.11, 1.12, 83.1, 83.2, 90.1],
            "RA": [3.5, 3.6, 3.7, 3.8, 3.9],
            "Dec": [-30.4] * 5,
            "z": [2.6, 2.6, 4.0, 4.0, 0.0],
            "Redshift_type": ["spec", "spec", "phot", "phot", "None"],
        }
    )
    t = image_list(raw)
    assert list(t["image_id"]) == ["1.1a", "1.1b", "83a", "83b"]
    assert list(t["z_spec"]) == [True, True, False, False]
    assert np.all(np.isfinite(t["ra"]))


def test_ambiguous_ids_fail_loudly():
    with pytest.raises(ValueError):
        lenstool_id("34.1" + "0")
    raw = Table(
        {
            "ID": [1.11, 1.11],
            "RA": [3.5, 3.6],
            "Dec": [-30.4] * 2,
            "z": [2.6, 2.6],
            "Redshift_type": ["spec", "spec"],
        }
    )
    with pytest.raises(ValueError):
        image_list(raw)


def test_gate_k_wls_decides_with_the_unchanged_gate(monkeypatch):
    import a2744_uncoverref as ar

    sp = Table({"image_id": ["1a", "1b", "2a", "2b", "2c"], "system": ["1", "1", "2", "2", "2"]})
    monkeypatch.setattr(ar.lensmodel, "backtrace_images", lambda model, t, zml: t)

    def fake_imageplane(model, grid, bt, d=(0.2, 0.3, 0.1, 0.4, 0.5)):
        return Table(
            {
                "image_id": bt["image_id"],
                "system": bt["system"],
                "dtheta_arcsec": list(d),
                "n_predicted": [3, 3, 2, 2, 2],
                "shared_match": [False, False, False, True, True],
            }
        )

    monkeypatch.setattr(ar.au, "imageplane", fake_imageplane)
    _, k = ar.gate_k_wls(None, None, sp)
    assert k["gate_K"]["passes"] and k["gate_K"]["n_images"] == 5
    assert k["not_reproduced"] == {} and k["shared_match"] == 2
    assert k["extra_predicted"] == {"1": (2, 3)} and k["fewer_predicted"] == {"2": (3, 2)}

    # one image beyond 3.0'' (80 % matched still passes on fraction) but rms > 1.0'' fails
    monkeypatch.setattr(
        ar.au, "imageplane", lambda m, g, bt: fake_imageplane(m, g, bt, (0.2, 0.3, 0.1, 0.4, 5.0))
    )
    _, k = ar.gate_k_wls(None, None, sp)
    assert not k["gate_K"]["passes"] and k["gate_K"]["matched_fraction"] == 0.8
    assert k["not_reproduced"] == {"2c": 5.0}

    # an unmatched (NaN) image counts against the fraction
    monkeypatch.setattr(
        ar.au,
        "imageplane",
        lambda m, g, bt: fake_imageplane(m, g, bt, (0.2, np.nan, 0.1, np.nan, 0.5)),
    )
    _, k = ar.gate_k_wls(None, None, sp)
    assert not k["gate_K"]["passes"] and k["gate_K"]["matched_fraction"] == 0.6
