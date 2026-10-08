from __future__ import annotations

import numpy as np
import pytest
from astropy.table import Table

from jwst_anomaly import signatures as sg


def test_registry_holds_the_searchable_d047_signatures():
    assert sorted(sg.REGISTRY) == ["W1", "W2", "W3", "W5"]
    assert [s.code for s in sg.for_kind("light_curve")] == ["W3"]
    assert {s.code for s in sg.for_kind("catalogue")} == {"W1", "W2", "W5"}
    w3 = sg.get("W3")
    lc = w3.predict(np.array([0.0, 1.0]), 0.0, 10.0, 0.5, n=1.0, sign=-1, rho=0.0)
    assert lc.shape == (2,)  # the closed form is callable through the registry
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
