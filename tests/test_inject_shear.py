"""Tests for scripts/inject_shear.py (synthetic inputs only)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(_DIR))
sys.path.insert(0, str(Path(__file__).resolve().parent))
_spec = importlib.util.spec_from_file_location("inject_shear", _DIR / "inject_shear.py")
ish = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ish)

from test_inject_radial import SCALE, _model, _shapes  # noqa: E402


def _field(n=900, seed=0):
    rng = np.random.default_rng(seed)
    model = _model()
    x, y = rng.uniform(-40, 40, (2, n))
    shapes = _shapes(
        model,
        x,
        y,
        pa_obs=rng.uniform(0, 180, n),
        ellipticity=rng.uniform(0.1, 0.5, n),
        semimajor_px=rng.uniform(3.0, 6.0, n),
    )
    args = SimpleNamespace(
        **(
            ish.SHEAR_DEFAULTS
            | {"max_radius": 40.0, "recover_tol": 2.0, "n_random": 40, "responsivity": 1.0}
        )
    )
    return model, shapes, args


def test_shear_defaults_come_from_the_screen():
    assert ish.SHEAR_DEFAULTS == ish.es.shear_defaults()
    assert (
        ish.SHEAR_DEFAULTS["filter"] == "schirmer" and ish.SHEAR_DEFAULTS["aperture_arcsec"] == 10
    )


def test_a_massive_w1_lens_is_recovered_and_no_lens_is_not():
    model, shapes, args = _field()
    inj = ish.ShearInjector(model, shapes, args, psf_sigma=1.0)
    assert inj.r == 1.0  # synthetic shapes carry no cluster shear to measure R from
    inj.draw_null(np.random.default_rng(1))
    big = np.where(inj.background, 6.0, np.nan)  # theta_E 6": dozens of radial images
    hit = inj.trial(15.0, -10.0, big, SCALE)
    assert hit["n_lensed"] > 5 and hit["n_image_used"] > 5
    assert hit["recovered"] and hit["s_best"] > np.max(inj.null)
    none = np.full(len(shapes), np.nan)  # no lensable rows: nothing painted
    miss = inj.trial(15.0, -10.0, none, SCALE)
    assert miss["n_lensed"] == 0 and miss["s_best"] < hit["s_best"]
