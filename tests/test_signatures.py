from __future__ import annotations

import numpy as np
import pytest
from astropy.table import Table

from jwst_anomaly import signatures as sg


def test_registry_holds_the_searchable_d047_signatures():
    assert sorted(sg.REGISTRY) == ["W1", "W2", "W3", "W5"]
    assert [s.code for s in sg.for_kind("light_curve")] == ["W3"]
    assert {s.code for s in sg.for_kind("catalogue")} == {"W1", "W2"}
    assert [s.code for s in sg.for_kind("count_map")] == ["W5"]
    w3 = sg.get("W3")
    lc = w3.predict(np.array([0.0, 1.0]), 0.0, 10.0, 0.5)
    assert lc[0] == 0.0  # negative mass bound: u0 = 0.5 is inside the umbra (u < 2)
    assert sg.get("W2").lens == {"n": 2.0, "sign": 1} and sg.get("W1").lens["sign"] == -1
    w1 = sg.get("W1").inject([3.0], [0.0], 1.0)
    assert len(w1) == 2 and (w1["dx"] > 0).all()  # both images on the source's side
    assert sg.get("W5").limits_doc == "docs/exotic_limits.md"
    assert "D-063" in sg.get("W5").decisions
    with pytest.raises(KeyError, match="registered"):
        sg.get("W9")
    with pytest.raises(ValueError, match="already registered"):
        sg.register(w3)
    with pytest.raises(ValueError, match="unknown data kinds"):
        sg.Signature("WX", "x", ("radio",), None, None, (), (), "", ())


def test_standard_light_curve_sorts_drops_bad_rows_and_labels_provenance():
    lc = sg.standard_light_curve(
        [3.0, 1.0, 2.0, np.nan],
        [18.0, 17.5, np.nan, 19.0],
        [0.01, 0.02, 0.01, 0.01],
        "I",
        source="test",
        time_system="HJD - 2450000",
    )
    assert list(lc["time"]) == [1.0, 3.0]
    assert list(lc["band"]) == ["I", "I"]
    assert lc.meta["provenance"] == "observed" and lc.meta["time_system"] == "HJD - 2450000"
    assert lc.meta["n_dropped"] == 2
    with pytest.raises(ValueError, match="one shape"):
        sg.standard_light_curve([1.0], [1.0, 2.0], [0.1], "I", source="t", time_system="t")


def test_adapters_are_structural_protocols():
    class Toy:
        name = "toy"

        def events(self):
            return Table({"event_id": ["a"], "ra": [0.0], "dec": [0.0]})

        def light_curve(self, event_id):
            return sg.standard_light_curve([0.0], [18.0], [0.01], "I", "toy", "MJD")

        def efficiency(self, t_e_days):
            return None

    assert isinstance(Toy(), sg.LightCurveSurvey)
    assert not isinstance(Toy(), sg.CatalogueSurvey)
