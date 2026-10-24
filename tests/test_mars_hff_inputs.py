import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import mars_hff_inputs as mh  # noqa: E402

CAT = (
    "# Strong lensing multiple images ...\n"
    "# The redshifts with ** indicate the redshift from the MARS algorithm.\n\n"
    "# Multiple images of Abell 370\t\t\t\n"
    "# ID\t\tR.A. (J2000)\t\tDecl. (J2000)\t\tz\t\tClass\t\tDelta_rms (arcsec)\n"
    "1.1\t\t39.967047\t\t-1.5769172\t\t0.8041\t\tgold\t\t0.01\n"
    "1.2\t\t39.976273\t\t-1.5760558\t\t0.8041\t\tgold\t\t0.0059\n"
    "8.1\t\t39.964485\t\t-1.5698065\t\t3.3379^{**} \tsilver\t\t0.0021\n"
    "# Multiple image of Abell S1063\n"
    "2.1\t\t342.18\t\t-44.53\t\t1.2\t\tgold\t\t0.02\n"
    "2.2\t\t342.19\t\t-44.53\t\t1.2\t\tgold\t\t0.02\n"
    "3.1\t\t342.18\t\t-44.531\t\t2.0^{**}\t\tgold\t\t0.02\n"
    "3.2\t\t342.19\t\t-44.531\t\t2.0^{**}\t\tgold\t\t0.02\n"
    "4.1\t\t342.18\t\t-44.532\t\t3.0\t\tgold\t\t0.02\n"
)


def test_parse_catalogue_fields_and_fitted_z():
    t = mh.parse_catalogue(CAT)
    assert list(t["field"][:3]) == ["Abell 370"] * 3
    assert list(t["field"][3:]) == ["Abell S1063"] * 5
    assert list(t["fitted_z"]) == [False, False, True, False, False, True, True, False]
    assert np.isclose(t["z"][2], 3.3379)
    assert list(t["system"][:3]) == ["1", "1", "8"]


def test_spec_images_gold_specz_min_two():
    t = mh.parse_catalogue(CAT)
    s = mh.spec_images(t, "Abell S1063")
    # system 3 has fitted z, system 4 has one image
    assert list(s["image_id"]) == ["2.1", "2.2"]
    assert list(mh.spec_images(t, "Abell 370")["image_id"]) == ["1.1", "1.2"]


def _header(n=41, pix=1.0):
    h = fits.Header()
    h.update(
        CTYPE1="RA---TAN",
        CTYPE2="DEC--TAN",
        CRVAL1=39.97,
        CRVAL2=-1.58,
        CRPIX1=(n + 1) / 2,
        CRPIX2=(n + 1) / 2,
        CDELT1=-pix / 3600,
        CDELT2=pix / 3600,
    )
    return h


def _sis(n=41, pix=1.0, theta_e=8.0):
    """Deflection of an SIS (arcsec at D_ds / D_s = 1), centred on the map."""
    c = (n - 1) / 2
    yy, xx = np.mgrid[0:n, 0:n]
    dx, dy = (xx - c) * pix, (yy - c) * pix
    r = np.hypot(dx, dy)
    r[r == 0] = 1e-9
    return theta_e * dx / r, theta_e * dy / r


def _images_for(model, z, r_pix):
    """Two images on the x axis at +-r_pix pixels from the centre, as RA / Dec."""
    from astropy.wcs.utils import pixel_to_skycoord

    sc = pixel_to_skycoord(np.array([20 + r_pix, 20 - r_pix]), np.array([20.0, 20.0]), model.wcs)
    t = Table()
    t["image_id"] = ["1.1", "1.2"]
    t["system"] = ["1", "1"]
    t["ra"], t["dec"] = sc.ra.deg, sc.dec.deg
    t["z"] = [z, z]
    t["a"] = [1.0, 1.0]
    return t


def test_coverage_field_fraction():
    ax, ay = _sis()
    m = mh.build_model(ax, ay, _header(), False, 0.375, "x", "t")
    t = _images_for(m, 2.0, 5.0)
    far = t[:1].copy()
    far["ra"], far["image_id"] = [far["ra"][0] + 0.05], ["9.9"]
    from astropy.table import vstack

    c = mh.coverage_field(m, vstack([t, far]), min_inside=0.6)
    assert c["inside"] == 2 and c["outside"] == ["9.9"] and c["passes"]
    assert not mh.coverage_field(m, vstack([t, far]), min_inside=0.9)["passes"]


def test_unit_readings_pick_arcsec_for_sis():
    """A symmetric pair at r = theta_E * D_ds/D_s back-traces to one source in the arcsec reading
    only (pixel 2" makes the pixel reading twice too strong)."""
    ax, ay = _sis(pix=2.0)
    a = mh.build_model(ax, ay, _header(pix=2.0), False, 0.375, "x", "A")
    p = mh.build_model(ax, ay, _header(pix=2.0), True, 0.375, "x", "P")
    assert np.isclose(p.pixel_arcsec, 2.0)
    zs = 2.0
    w = (
        1.0
        - (
            a.cosmology.comoving_transverse_distance(0.375)
            / a.cosmology.comoving_transverse_distance(zs)
        ).value
    )
    t = _images_for(a, zs, 8.0 * w / 2.0)  # pixel units: radius in pixels = arcsec / 2
    ua, up = mh.per_system_u(a, t)["1"], mh.per_system_u(p, t)["1"]
    assert ua < 0.05 and up > 1.0


def test_record_outcome_counts_all_three():
    f = {
        "a": {"check_C": {"passes": True}, "check_U": {"passes": True}, "gate_K": {"passes": True}},
        "b": {
            "check_C": {"passes": True},
            "check_U": {"passes": True},
            "gate_K": {"passes": False},
        },
        "c": {"check_C": {"passes": False}},
    }
    o = mh.record_outcome(f)
    assert o["usable_fields"] == ["a"] and o["n_usable"] == 1 and o["n_fields"] == 3
