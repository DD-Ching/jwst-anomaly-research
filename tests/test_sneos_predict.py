import sys
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import sneos_predict as sp  # noqa: E402

COSMO = FlatLambdaCDM(H0=70.0, Om0=0.3)


def test_sis_delay_matches_analytic():
    # SIS: ψ = θ_E |θ|; images at β ± θ_E on the axis; Δt = D_dt/c · 2 β θ_E
    theta_e, beta = 10.0, 1.5
    x = np.array([beta + theta_e, beta - theta_e])
    psi = theta_e * np.abs(x)
    t = sp.delay_days(sp.fermat(x, np.zeros(2), beta, 0.0, psi), 0.35, 5.0, COSMO)
    analytic = sp.delay_days(np.array([2.0 * beta * theta_e]), 0.35, 5.0, COSMO)[0]
    assert np.isclose(t[1] - t[0], analytic, rtol=1e-12)
    assert t[1] > t[0]  # the image nearer the centre (saddle) arrives later


def test_delay_scale_is_physical():
    # 1 arcsec² of Fermat potential at z_l 0.35, z_s 5: tens of days (D_dt ~ 2.6 Gpc)
    d = sp.delay_days(np.array([1.0]), 0.35, 5.0, COSMO)[0]
    assert 20.0 < d < 60.0


def test_potential_scale_recovers_factor():
    j, i = np.indices((60, 60), dtype=float)
    psi = 0.01 * ((i - 30) ** 2 + 0.5 * (j - 25) ** 2)
    dj, di = np.gradient(psi)
    mask = np.ones_like(psi, bool)
    mask[0, :] = mask[-1, :] = mask[:, 0] = mask[:, -1] = False
    mask &= np.hypot(di, dj) > 0.05
    g = sp.potential_scale(3.0 * di, 3.0 * dj, psi, mask)
    assert np.isclose(g["s"], 3.0) and g["passes"]
    bad = sp.potential_scale(3.0 * dj, 3.0 * di, psi, mask)  # swapped axes
    assert not bad["passes"]


def test_match_one_to_one_within_radius():
    m = sp.match(np.array([0.0, 0.8, 5.0]), np.zeros(3), np.array([0.1, 0.7]), np.zeros(2))
    assert m == {0: 0, 1: 1}
    assert sp.match(np.array([0.0]), np.zeros(1), np.array([0.6]), np.zeros(1)) == {}
