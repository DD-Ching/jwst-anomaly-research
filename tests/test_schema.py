import importlib

import pytest
from astropy.table import Table

from jwst_anomaly import schema

STAGE_MODULES = [
    "paths",
    "query",
    "acquire",
    "catalog",
    "features",
    "rank",
    "cutouts",
    "crossmatch",
    "candidates",
    "provenance",
    "pipeline",
    "cli",
]


@pytest.mark.parametrize("name", STAGE_MODULES)
def test_stage_modules_import(name):
    importlib.import_module(f"jwst_anomaly.{name}")


def test_band_column():
    assert schema.band_column("F200W", "aper50_abmag") == "f200w_aper50_abmag"


def test_validate_accepts_contract():
    t = Table({"source_uid": ["a"], "score": [1.0], "rank": [1]})
    t.meta.update(provenance=schema.Provenance.MODEL_PREDICTION.value, source="unit test")
    assert schema.validate(t, schema.SCORE_COLUMNS) is t


@pytest.mark.parametrize(
    "cols, meta, match",
    [
        ({"source_uid": ["a"]}, {"provenance": "derived", "source": "x"}, "missing"),
        ({"source_uid": ["a"]}, {"provenance": "made-up", "source": "x"}, "provenance"),
        ({"source_uid": ["a"]}, {"provenance": "derived"}, "source"),
    ],
)
def test_validate_rejects(cols, meta, match):
    t = Table(cols)
    t.meta.update(meta)
    required = schema.SCORE_COLUMNS if match == "missing" else schema.FEATURE_ID_COLUMNS
    with pytest.raises(ValueError, match=match):
        schema.validate(t, required)
