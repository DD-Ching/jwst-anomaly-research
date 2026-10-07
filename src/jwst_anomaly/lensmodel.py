"""Published Lenstool lens models evaluated at sky positions (M3 lens-model stage, D-024).

A Lenstool ``.par`` file (e.g. ``best.par`` of a published model) is parsed with
:func:`parse_lenstool_par` and evaluated by :class:`LensModel` at sky positions for a source
redshift: deflection, convergence, shear, magnification and the expected arc orientation.
Every table returned carries ``meta["provenance"] = "model_prediction"``.

Supported: one lens plane whose potentials are all dPIE/PIEMD (Lenstool ``profile 81``), in
absolute coordinates (``runmode`` ``reference 3 ra dec``), with a flat LambdaCDM cosmology.
Anything else raises :class:`UnsupportedModelError`. Scaling-relation galaxy catalogs
(``potfile``) are not expanded: use the ``best.par`` of a run, which lists every potential.

Lenstool conventions used here (checked against the Mahler et al. 2022 ICLv2 convergence map,
D-024):

* frame: ``x = -(ra - ra0) cos(dec0) * 3600`` (arcsec, positive towards West) and
  ``y = (dec - dec0) * 3600`` (positive towards North), relative to the reference point;
* ``angle_pos``: degrees counter-clockwise from the +x axis (West towards North);
* ``ellipticity`` (``ellipticite``): Lenstool's mass ellipticity ``(a^2 - b^2) / (a^2 + b^2)``,
  converted to the potential's ``epot = (1 - q) / (1 + q)`` with ``q = sqrt((1 - e) / (1 + e))``;
* ``v_disp``: Lenstool's fiducial velocity dispersion; lens strength
  ``b0 = 6 pi (v_disp / c)^2 D_LS / D_S`` (radians);
* ``core_radius`` / ``cut_radius`` in arcsec (the ``_kpc`` variants are converted with the
  model's cosmology at ``z_lens`` when the arcsec values are absent).
"""

# The dPIE deflection (``_ci05f``), Hessian (``_mdci05``) and the Lenstool parameter conversion
# (ellipticity -> epot, v_disp -> b0) are ported to plain numpy from PyAutoGalaxy
# (autogalaxy==2026.10.7.1, autogalaxy/profiles/mass/total/dual_pseudo_isothermal_mass.py),
# which ports them from Lenstool's C code (Kassiola & Kovner 1993; Eliasdottir et al. 2007).
# That file is distributed under the MIT License:
#
#   MIT License
#
#   Copyright (c) 2019 James Nightingale
#
#   Permission is hereby granted, free of charge, to any person obtaining a copy
#   of this software and associated documentation files (the "Software"), to deal
#   in the Software without restriction, including without limitation the rights
#   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#   copies of the Software, and to permit persons to whom the Software is
#   furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in all
#   copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#   SOFTWARE.

from __future__ import annotations

import hashlib
import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table

from jwst_anomaly import schema

C_KM_S = 299792.458
ARCSEC_PER_RAD = 648000.0 / np.pi
# Lenstool clamps the potential ellipticity away from 0 (the KK93 integral is singular there).
MIN_EPS, MAX_EPS = 1e-5, 0.99999

#: Mahler et al. (2022; ApJ 945, 49, 2023) SMACS J0723 Lenstool model, version ICLv2 (the
#: repository README's final model), CC0-1.0, pinned commit (tag v2). name -> (url, sha256).
MAHLER22_COMMIT = "f36a41c365865df252a6da08cd79955f9cd68f11"
_MAHLER22_RAW = (
    "https://raw.githubusercontent.com/guillaumemahler/SMACS0723-mahler2022/" + MAHLER22_COMMIT
)
SMACS0723_MAHLER22_ICLV2: dict[str, tuple[str, str]] = {
    "best.par": (
        f"{_MAHLER22_RAW}/ICLv2/best.par",
        "73cc7881b224e4baa492eae56269d90c154521ffc7e2dd66ed174773098a7715",
    ),
    "arcs.dat": (
        f"{_MAHLER22_RAW}/ICLv2/arcs.dat",
        "af7e5368e4bdaa8903fd86239fc0d5e2eb9ad7da23e0cf3f2c0d86e0e49040b9",
    ),
    "input.par": (
        f"{_MAHLER22_RAW}/ICLv2/input.par",
        "db30f227ea69714d5b069863fd849af6fd58337967a297b90cbc8fc1010966dd",
    ),
    # Best-model convergence (MCMC sample 0000), 3000 x 3000 px of 0.01334", D_LS/D_S = 1.
    "kappa_0000": (
        "https://github.com/guillaumemahler/SMACS0723-mahler2022/raw/"
        f"{MAHLER22_COMMIT}/ICLv2/tmp_k/0000_k.fits.tar.xz",
        "4de5d733587dd4dea12233c25069d8e1307936fcd50de316cc4eb967a53edc2c",
    ),
}

# French Lenstool keywords (input files) -> the English ones written in best.par.
_ALIASES = {
    "potentiel": "potential",
    "profil": "profile",
    "ellipticite": "ellipticity",
    "cosmologie": "cosmology",
    "fini": "finish",
    "grille": "grid",
    "champ": "field",
    "nlentille": "nlens",
    "nombre": "number",
    "polaire": "polar",
}
_POTENTIAL_KEYS = (
    "x_centre",
    "y_centre",
    "ellipticity",
    "angle_pos",
    "core_radius",
    "core_radius_kpc",
    "cut_radius",
    "cut_radius_kpc",
    "v_disp",
    "z_lens",
)
# Potential keywords that carry no mass information (``mag`` feeds Lenstool's scaling relations,
# whose result best.par already lists). Any other unknown keyword raises UnsupportedModelError.
_IGNORED_POTENTIAL_KEYS = ("profile", "mag")
SUPPORTED_PROFILES = (81,)


class UnsupportedModelError(ValueError):
    """A Lenstool model feature this module does not implement (it is never approximated)."""


_INTEGER_ID = re.compile(r"^[+-]?\d+(?:\.0*)?$")


def system_key(system_id: str | float) -> str:
    """Canonical multiple-image system id.

    Integer-valued ids (``"4"``, ``"4.0"``, ``4.0``; ``best.par`` writes ``z_m_limit`` systems as
    floats) give ``"4"``; any other id is kept as written (``"4.10"``, ``"c2"``, ``"1a"``).
    """
    text = str(system_id).strip()
    return str(int(float(text))) if _INTEGER_ID.match(text) else text


def image_redshifts(images: Table, z_m_limit: dict[str, float]) -> np.ndarray:
    """Redshift used for each image of an ``arcs.dat`` table.

    The model's fixed ``z_m_limit`` value for the image's system comes first, because Lenstool
    applies it to every image of that system. Otherwise the catalogued ``z`` is used when > 0,
    else NaN.
    """
    return np.array(
        [
            z_m_limit.get(str(s), zc if zc > 0 else np.nan)
            for zc, s in zip(images["z"], images["system"], strict=True)
        ],
        float,
    )


def parse_lenstool_par(path: str | Path) -> dict[str, Any]:
    """Parse a Lenstool ``.par`` file (French or English keywords).

    Returns a dict with:

    * ``reference``: ``{"mode", "ra", "dec"}`` (only mode 3, absolute degrees, is accepted);
    * ``cosmology``: ``{"H0", "Om0"}`` (flat LambdaCDM only);
    * ``potentials``: one dict per ``potential`` section with ``name``, ``profile`` and the
      keys in ``_POTENTIAL_KEYS`` (NaN when absent);
    * ``potfiles``: names of ``potfile`` sections (not expanded);
    * ``z_m_limit``: fixed (best-fit) redshifts of multiple-image systems, by :func:`system_key`;
    * ``sigpos_arcsec``: the image-plane position error ``sigposArcsec``, or None;
    * ``source`` and ``sha256`` of the file.

    Raises :class:`UnsupportedModelError` for a potential profile other than 81 (dPIE), a
    reference mode other than 3, or a non-flat / non-LambdaCDM cosmology, and ``ValueError``
    for a malformed file.
    """
    path = Path(path)
    raw = path.read_bytes()
    sections: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    for n, line in enumerate(raw.decode("latin-1").splitlines(), 1):
        stripped = line.split("#", 1)[0].strip()
        if not stripped:
            continue
        tokens = stripped.split()
        key = _ALIASES.get(tokens[0].lower(), tokens[0].lower())
        if current is None:
            if key == "finish":
                break
            current = {"kind": key, "name": " ".join(tokens[1:]), "line": n, "entries": []}
        elif key == "end":
            sections.append(current)
            current = None
        else:
            current["entries"].append((key, tokens[1:], n))
    if current is not None:
        raise ValueError(f"{path}:{current['line']}: section {current['kind']!r} has no 'end'")

    out: dict[str, Any] = {
        "source": str(path),
        "sha256": hashlib.sha256(raw).hexdigest(),
        "reference": None,
        "cosmology": None,
        "potentials": [],
        "potfiles": [],
        "z_m_limit": {},
        "sigpos_arcsec": None,
    }
    for sec in sections:
        kind = sec["kind"]
        entries = sec["entries"]
        if kind == "runmode":
            for key, vals, n in entries:
                if key == "reference":
                    if len(vals) < 3 or vals[0] != "3":
                        raise UnsupportedModelError(
                            f"{path}:{n}: 'reference {' '.join(vals)}' - only mode 3 "
                            "(absolute RA/Dec in degrees) is supported"
                        )
                    out["reference"] = {
                        "mode": 3,
                        "ra": float(vals[1]),
                        "dec": float(vals[2]),
                    }
        elif kind == "cosmology":
            out["cosmology"] = _parse_cosmology(path, entries)
        elif kind == "image":
            for key, vals, n in entries:
                if key == "z_m_limit":
                    # z_m_limit <n> <system> <flag> <z or zmin> [<zmax> <step>]; flag 0 = fixed.
                    if len(vals) < 4:
                        raise ValueError(f"{path}:{n}: malformed z_m_limit {vals}")
                    if vals[2] == "0":
                        out["z_m_limit"][system_key(vals[1])] = float(vals[3])
                elif key.lower() == "sigposarcsec":
                    out["sigpos_arcsec"] = float(vals[0])
        elif kind == "potential":
            out["potentials"].append(_parse_potential(path, sec))
        elif kind.startswith("potfile"):
            out["potfiles"].append(kind)
    if out["reference"] is None:
        raise UnsupportedModelError(f"{path}: no 'runmode' reference point")
    return out


def _parse_cosmology(path: Path, entries: list[tuple[str, list[str], int]]) -> dict[str, float]:
    vals = {key.lower(): float(v[0]) for key, v, _ in entries if v}
    h0, om = vals.get("h0"), vals.get("omegam")
    ox = vals.get("omegax", 1.0 - (om or 0.0))
    w0, wa = vals.get("wx", -1.0), vals.get("wa", 0.0)
    if h0 is None or om is None:
        raise ValueError(f"{path}: cosmology section needs H0 and omegaM")
    if abs(om + ox - 1.0) > 1e-6 or w0 != -1.0 or wa != 0.0:
        raise UnsupportedModelError(
            f"{path}: only flat LambdaCDM is supported (omegaM={om}, omegaX={ox}, wX={w0}, wa={wa})"
        )
    return {"H0": h0, "Om0": om}


def _parse_potential(path: Path, sec: dict[str, Any]) -> dict[str, Any]:
    name = sec["name"] or f"line{sec['line']}"
    values = {key: vals for key, vals, _ in sec["entries"]}
    if "profile" not in values:
        raise UnsupportedModelError(f"{path}:{sec['line']}: potential {name} has no profile")
    profile = int(float(values["profile"][0]))
    if profile not in SUPPORTED_PROFILES:
        raise UnsupportedModelError(
            f"{path}:{sec['line']}: potential {name} has Lenstool profile {profile}; only "
            f"{SUPPORTED_PROFILES} (dPIE/PIEMD) are supported"
        )
    unknown = sorted(set(values) - set(_POTENTIAL_KEYS) - set(_IGNORED_POTENTIAL_KEYS))
    if unknown:
        raise UnsupportedModelError(
            f"{path}:{sec['line']}: potential {name} has keywords this module does not interpret: "
            f"{unknown}"
        )
    pot: dict[str, Any] = {"name": name, "profile": profile}
    for key in _POTENTIAL_KEYS:
        pot[key] = float(values[key][0]) if key in values else float("nan")
    return pot


@dataclass(frozen=True)
class DPIE:
    """One Lenstool dPIE (PIEMD, profile 81) in Lenstool's own parameters.

    ``x``/``y`` are arcsec in the model frame (see the module docstring), ``angle_pos`` degrees
    counter-clockwise from +x, ``r_core``/``r_cut`` arcsec, ``v_disp`` km/s.
    """

    name: str
    x: float
    y: float
    ellipticity: float
    angle_pos: float
    r_core: float
    r_cut: float
    v_disp: float
    z_lens: float

    @property
    def eps(self) -> float:
        """Potential ellipticity ``(1 - q) / (1 + q)`` from Lenstool's mass ellipticity."""
        e = self.ellipticity
        q = np.sqrt((1.0 - e) / (1.0 + e))
        return float(np.clip((1.0 - q) / (1.0 + q), MIN_EPS, MAX_EPS))

    @property
    def b0(self) -> float:
        """Lens strength in arcsec at D_LS/D_S = 1: ``6 pi (v_disp / c)^2`` radians."""
        return 6.0 * np.pi * (self.v_disp / C_KM_S) ** 2 * ARCSEC_PER_RAD

    def _local(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, float]:
        t = np.deg2rad(self.angle_pos)
        c, s = np.cos(t), np.sin(t)
        dx, dy = x - self.x, y - self.y
        return c * dx + s * dy, -s * dx + c * dy, c, s

    def kappa(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Convergence at D_LS/D_S = 1 (analytic, Eliasdottir et al. 2007 eq. A9 form)."""
        xr, yr, _, _ = self._local(np.asarray(x, float), np.asarray(y, float))
        eps, a, s = self.eps, self.r_core, self.r_cut
        rem2 = xr**2 / (1.0 + eps) ** 2 + yr**2 / (1.0 - eps) ** 2
        return (
            0.5
            * self.b0
            * s
            / (s - a)
            * (1.0 / np.sqrt(a * a + rem2) - 1.0 / np.sqrt(s * s + rem2))
        )

    def deflection(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Deflection (arcsec) at D_LS/D_S = 1 in the model frame."""
        xr, yr, c, s = self._local(np.asarray(x, float), np.asarray(y, float))
        z = _ci05f(xr, yr, self.eps, self.r_core, self.r_cut)
        f = self.b0 * self.r_cut / (self.r_cut - self.r_core)
        axr, ayr = f * z.real, f * z.imag
        return c * axr - s * ayr, s * axr + c * ayr

    def hessian(self, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Second derivatives of the potential (psi_xx, psi_xy, psi_yy) at D_LS/D_S = 1."""
        xr, yr, c, s = self._local(np.asarray(x, float), np.asarray(y, float))
        f = self.b0 * self.r_cut / (self.r_cut - self.r_core)
        a_xx, a_xy, a_yy = _mdci05(xr, yr, self.eps, self.r_core)
        b_xx, b_xy, b_yy = _mdci05(xr, yr, self.eps, self.r_cut)
        hxx, hxy, hyy = f * (a_xx - b_xx), f * (a_xy - b_xy), f * (a_yy - b_yy)
        # Rotate the tensor back to the model frame: H = R H' R^T.
        rxx = c * c * hxx - 2 * c * s * hxy + s * s * hyy
        ryy = s * s * hxx + 2 * c * s * hxy + c * c * hyy
        rxy = c * s * (hxx - hyy) + (c * c - s * s) * hxy
        return rxx, rxy, ryy


def _ci05f(x: np.ndarray, y: np.ndarray, eps: float, rcore: float, rcut: float) -> np.ndarray:
    """Deflection of a dPIE per unit ``b0 rcut / (rcut - rcore)``, as ``alpha_x + i alpha_y``.

    Kassiola & Kovner (1993) eq. 4.1.2 for core radius ``rcore`` minus the same for ``rcut``,
    in the profile's frame (major axis along x). Lenstool ``ci05f``, via PyAutoGalaxy.
    """
    sqe = np.sqrt(eps)
    q = (1.0 - eps) / (1.0 + eps)
    rem2 = x * x / (1.0 + eps) ** 2 + y * y / (1.0 - eps) ** 2
    zci = -0.5j * (1.0 - eps * eps) / sqe
    znum_rc = q * x + 1j * (2.0 * sqe * np.sqrt(rcore * rcore + rem2) - y / q)
    zden_rc = x + 1j * (2.0 * rcore * sqe - y)
    znum_rcut = q * x + 1j * (2.0 * sqe * np.sqrt(rcut * rcut + rem2) - y / q)
    zden_rcut = x + 1j * (2.0 * rcut * sqe - y)
    return zci * np.log((znum_rc * zden_rcut) / (zden_rc * znum_rcut))


def _mdci05(
    x: np.ndarray, y: np.ndarray, eps: float, rcore: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Hessian (psi_xx, psi_xy, psi_yy) per unit ``b0`` of a PIEMD with core ``rcore``.

    Kassiola & Kovner (1993) eq. 4.1.4 in the profile's frame. Lenstool ``mdci05``, via
    PyAutoGalaxy ``_mdci05``.
    """
    sqe = np.sqrt(eps)
    q = (1.0 - eps) / (1.0 + eps)
    qi = 1.0 / q
    cxro = (1.0 + eps) ** 2
    cyro = (1.0 - eps) ** 2
    ci = 0.5 * (1.0 - eps * eps) / sqe
    wrem = np.sqrt(rcore * rcore + x * x / cxro + y * y / cyro)
    num1 = 2.0 * sqe * wrem - y * qi
    den1 = q * q * x * x + num1 * num1
    num2 = 2.0 * rcore * sqe - y
    den2 = x * x + num2 * num2
    dxx = ci * (
        q * (2.0 * sqe * x * x / cxro / wrem - 2.0 * sqe * wrem + y * qi) / den1 + num2 / den2
    )
    dxy = ci * ((2.0 * sqe * x * y * q / cyro / wrem - x) / den1 + x / den2)
    dyy = ci * (
        (
            2.0 * sqe * wrem * qi
            - y * qi * qi
            - 4.0 * eps * y / cyro
            + 2.0 * sqe * y * y / cyro / wrem * qi
        )
        / den1
        - num2 / den2
    )
    return dxx, dxy, dyy


def axis_offset_deg(pa_a: Any, pa_b: Any) -> np.ndarray:
    """Smallest angle between two axes (orientations mod 180 deg), in [0, 90]."""
    d = np.mod(np.asarray(pa_a, float) - np.asarray(pa_b, float), 180.0)
    return np.minimum(d, 180.0 - d)


class LensModel:
    """A single-plane sum of dPIE potentials with a reference point and cosmology.

    Sky-position methods take ``ra``, ``dec`` (degrees; scalars or arrays) and a source
    redshift ``z_s`` (scalar or one per position) and return an ``astropy`` Table with
    ``meta["provenance"] = "model_prediction"``.
    """

    def __init__(
        self,
        components: Sequence[DPIE],
        ra0: float,
        dec0: float,
        cosmology: FlatLambdaCDM,
        source: str = "",
        sha256: str = "",
    ) -> None:
        if not components:
            raise ValueError("a lens model needs at least one potential")
        z_lens = {c.z_lens for c in components}
        if len(z_lens) != 1:
            raise UnsupportedModelError(f"multi-plane models are not supported (z_lens {z_lens})")
        self.components = tuple(components)
        self.z_lens = z_lens.pop()
        self.ra0, self.dec0 = float(ra0), float(dec0)
        self.cosmology = cosmology
        self.source = source
        self.sha256 = sha256
        self._cos0 = np.cos(np.deg2rad(self.dec0))

    @classmethod
    def from_par(cls, par: str | Path | dict[str, Any]) -> LensModel:
        """Build from a Lenstool ``.par`` path or the dict of :func:`parse_lenstool_par`."""
        parsed = par if isinstance(par, dict) else parse_lenstool_par(par)
        if parsed["potfiles"]:
            raise UnsupportedModelError(
                f"{parsed['source']}: potfile sections {parsed['potfiles']} are not expanded; "
                "use the run's best.par, which lists every potential"
            )
        if parsed["cosmology"] is None:
            raise UnsupportedModelError(f"{parsed['source']}: no cosmology section")
        cosmo = FlatLambdaCDM(H0=parsed["cosmology"]["H0"], Om0=parsed["cosmology"]["Om0"])
        comps = [_dpie_from_dict(p, cosmo, parsed["source"]) for p in parsed["potentials"]]
        ref = parsed["reference"]
        return cls(comps, ref["ra"], ref["dec"], cosmo, parsed["source"], parsed["sha256"])

    # --- geometry ---------------------------------------------------------------------------
    def to_frame(self, ra: Any, dec: Any) -> tuple[np.ndarray, np.ndarray]:
        """Sky (deg) -> model frame (arcsec): Lenstool's ``x = -dRA cos(dec0)``, ``y = dDec``."""
        dra = (np.asarray(ra, float) - self.ra0 + 180.0) % 360.0 - 180.0
        return -dra * self._cos0 * 3600.0, (np.asarray(dec, float) - self.dec0) * 3600.0

    def to_sky(self, x: Any, y: Any) -> tuple[np.ndarray, np.ndarray]:
        """Model frame (arcsec) -> sky (deg); the inverse of :meth:`to_frame`."""
        ra = self.ra0 - np.asarray(x, float) / 3600.0 / self._cos0
        return ra % 360.0, self.dec0 + np.asarray(y, float) / 3600.0

    @staticmethod
    def frame_angle_to_pa(phi_deg: Any) -> np.ndarray:
        """Direction angle in the model frame (deg CCW from +x = West) -> sky PA mod 180.

        PA is measured East of North. +x points West, so East = -cos(phi), North = sin(phi).
        """
        p = np.deg2rad(np.asarray(phi_deg, float))
        pa = np.mod(np.rad2deg(np.arctan2(-np.cos(p), np.sin(p))), 180.0)
        return np.where(np.isclose(pa, 180.0, rtol=0.0, atol=1e-9), 0.0, pa)  # keep [0, 180)

    def dls_ds(self, z_s: Any) -> np.ndarray:
        """Distance ratio D_LS / D_S for source redshift(s) ``z_s`` (0 in front of the lens).

        Flat cosmology: D_LS / D_S = 1 - D_M(z_lens) / D_M(z_s), with D_M the transverse comoving
        distance.
        """
        z = np.atleast_1d(np.asarray(z_s, float))
        out = np.where(np.isfinite(z), 0.0, np.nan)  # an unknown redshift stays unknown
        behind = np.isfinite(z) & (z > self.z_lens)
        if behind.any():
            d_l = self.cosmology.comoving_transverse_distance(self.z_lens)
            d_s = self.cosmology.comoving_transverse_distance(z[behind])
            out[behind] = 1.0 - (d_l / d_s).to_value(u.dimensionless_unscaled)
        return out.reshape(np.shape(z_s))

    # --- frame-level fields (D_LS/D_S = 1) --------------------------------------------------
    def kappa_xy(self, x: Any, y: Any) -> np.ndarray:
        """Analytic convergence at D_LS/D_S = 1 at model-frame positions (arcsec)."""
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        total = np.zeros(x.shape)
        for comp in self.components:
            total += comp.kappa(x, y)
        return total

    def fields_xy(self, x: Any, y: Any) -> dict[str, np.ndarray]:
        """Deflection and Hessian at D_LS/D_S = 1: ``alpha_x, alpha_y, psi_xx, psi_xy, psi_yy``."""
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        out = {k: np.zeros(x.shape) for k in ("alpha_x", "alpha_y", "psi_xx", "psi_xy", "psi_yy")}
        for comp in self.components:
            ax, ay = comp.deflection(x, y)
            hxx, hxy, hyy = comp.hessian(x, y)
            out["alpha_x"] += ax
            out["alpha_y"] += ay
            out["psi_xx"] += hxx
            out["psi_xy"] += hxy
            out["psi_yy"] += hyy
        return out

    # --- sky-level predictions --------------------------------------------------------------
    def evaluate(
        self, ra: Any, dec: Any, z_s: Any, fields: dict[str, np.ndarray] | None = None
    ) -> Table:
        """All predictions at sky positions for source redshift(s) ``z_s``.

        Columns: ``ra, dec, z_s, dls_ds, x, y`` (model frame, arcsec), ``alpha_x, alpha_y``
        (arcsec, model frame), ``beta_ra, beta_dec`` (source-plane position, deg), ``kappa``,
        ``gamma1, gamma2`` (model frame), ``gamma``, ``reduced_shear`` (|gamma| / |1 - kappa|),
        ``magnification`` (signed; negative = odd parity) and ``tangential_pa`` (deg E of N,
        mod 180: the major axis expected for the image of a small round source).

        ``fields`` may pass :meth:`fields_xy` already computed at these positions; the fields do not
        depend on ``z_s``, so callers evaluating several redshifts compute them once.
        """
        ra_a, dec_a = np.broadcast_arrays(
            np.atleast_1d(np.asarray(ra, float)), np.atleast_1d(np.asarray(dec, float))
        )
        zs = np.broadcast_to(np.asarray(z_s, float), ra_a.shape)
        scale = self.dls_ds(zs)
        x, y = self.to_frame(ra_a, dec_a)
        f = self.fields_xy(x, y) if fields is None else fields
        ax, ay = scale * f["alpha_x"], scale * f["alpha_y"]
        hxx, hxy, hyy = scale * f["psi_xx"], scale * f["psi_xy"], scale * f["psi_yy"]
        kappa = 0.5 * (hxx + hyy)
        g1, g2 = 0.5 * (hxx - hyy), hxy
        gamma = np.hypot(g1, g2)
        det = (1.0 - hxx) * (1.0 - hyy) - hxy * hxy
        with np.errstate(divide="ignore", invalid="ignore"):
            mu = 1.0 / det
            reduced = gamma / np.abs(1.0 - kappa)
        # Stretching axis: eigenvector of A = I - H with the smaller |eigenvalue|. That is the
        # shear axis 0.5 atan2(g2, g1) where kappa < 1, and the perpendicular where kappa > 1.
        phi = 0.5 * np.rad2deg(np.arctan2(g2, g1)) + np.where(kappa > 1.0, 90.0, 0.0)
        pa = np.where(gamma > 0, self.frame_angle_to_pa(phi), np.nan)
        beta_ra, beta_dec = self.to_sky(x - ax, y - ay)
        t = Table(
            {
                "ra": ra_a,
                "dec": dec_a,
                "z_s": np.asarray(zs, float),
                "dls_ds": scale,
                "x": x,
                "y": y,
                "alpha_x": ax,
                "alpha_y": ay,
                "beta_ra": beta_ra,
                "beta_dec": beta_dec,
                "kappa": kappa,
                "gamma1": g1,
                "gamma2": g2,
                "gamma": gamma,
                "reduced_shear": reduced,
                "magnification": mu,
                "tangential_pa": pa,
            }
        )
        for col in ("x", "y", "alpha_x", "alpha_y"):
            t[col].unit = u.arcsec
        for col in ("ra", "dec", "beta_ra", "beta_dec", "tangential_pa"):
            t[col].unit = u.deg
        t.meta.update(self._meta())
        return t

    def _subset(self, ra: Any, dec: Any, z_s: Any, cols: Sequence[str]) -> Table:
        full = self.evaluate(ra, dec, z_s)
        out = full[["ra", "dec", "z_s", "dls_ds", *cols]]
        out.meta.update(full.meta)
        return out

    def kappa(self, ra: Any, dec: Any, z_s: Any) -> Table:
        """Convergence at sky positions (column ``kappa``)."""
        return self._subset(ra, dec, z_s, ["kappa"])

    def gamma(self, ra: Any, dec: Any, z_s: Any) -> Table:
        """Shear (``gamma1``, ``gamma2`` in the model frame, ``gamma`` = modulus)."""
        return self._subset(ra, dec, z_s, ["gamma1", "gamma2", "gamma"])

    def magnification(self, ra: Any, dec: Any, z_s: Any) -> Table:
        """Signed magnification (column ``magnification``; negative = odd parity)."""
        return self._subset(ra, dec, z_s, ["magnification"])

    def deflection(self, ra: Any, dec: Any, z_s: Any) -> Table:
        """Deflection (arcsec, model frame) and the source-plane position ``beta_ra/dec``."""
        return self._subset(ra, dec, z_s, ["alpha_x", "alpha_y", "beta_ra", "beta_dec"])

    def tangential_pa(self, ra: Any, dec: Any, z_s: Any) -> Table:
        """Expected arc orientation (deg E of N, mod 180) with ``kappa`` and ``gamma``."""
        return self._subset(ra, dec, z_s, ["tangential_pa", "kappa", "gamma"])

    def _meta(self) -> dict[str, Any]:
        return {
            "provenance": schema.Provenance.MODEL_PREDICTION.value,
            "source": f"Lenstool model {self.source}",
            "model_sha256": self.sha256,
            "n_potentials": len(self.components),
            "z_lens": self.z_lens,
            "reference_radec": [self.ra0, self.dec0],
            "cosmology": (
                f"FlatLambdaCDM(H0={self.cosmology.H0.value:g}, Om0={self.cosmology.Om0:g})"
            ),
            "frame": "x = -dRA cos(dec0) * 3600 (West), y = dDec * 3600 (North), arcsec",
        }


def _dpie_from_dict(pot: dict[str, Any], cosmo: FlatLambdaCDM, source: str) -> DPIE:
    def radius(key: str) -> float:
        val, kpc = pot[key], pot[f"{key}_kpc"]
        per_arcsec = cosmo.kpc_proper_per_arcmin(pot["z_lens"]).to_value(u.kpc / u.arcmin) / 60.0
        if np.isfinite(val) and np.isfinite(kpc) and abs(kpc / per_arcsec - val) > 0.02 * val:
            # Both given and inconsistent: which one Lenstool used is ambiguous, so refuse.
            raise UnsupportedModelError(
                f"{source}: potential {pot['name']} {key} {val} arcsec disagrees with "
                f"{key}_kpc {kpc} ({kpc / per_arcsec:.4g} arcsec in the model cosmology)"
            )
        if np.isfinite(val):
            return float(val)
        if not np.isfinite(kpc):
            raise ValueError(f"{source}: potential {pot['name']} has neither {key} nor {key}_kpc")
        return float(kpc / per_arcsec)

    required = ("x_centre", "y_centre", "v_disp", "z_lens")
    missing = [k for k in required if not np.isfinite(pot[k])]
    if missing:
        raise ValueError(f"{source}: potential {pot['name']} lacks {missing}")
    r_core, r_cut = radius("core_radius"), radius("cut_radius")
    if not 0 <= r_core < r_cut:
        raise ValueError(f"{source}: potential {pot['name']} needs 0 <= core < cut radius")
    ell = pot["ellipticity"] if np.isfinite(pot["ellipticity"]) else 0.0
    angle = pot["angle_pos"] if np.isfinite(pot["angle_pos"]) else 0.0
    return DPIE(
        name=str(pot["name"]),
        x=float(pot["x_centre"]),
        y=float(pot["y_centre"]),
        ellipticity=float(ell),
        angle_pos=float(angle),
        r_core=r_core,
        r_cut=r_cut,
        v_disp=float(pot["v_disp"]),
        z_lens=float(pot["z_lens"]),
    )


def load_lenstool_images(path: str | Path) -> Table:
    """Lenstool ``arcs.dat`` (absolute coordinates): ``image_id, system, ra, dec, z``.

    ``z`` is the catalogued redshift (0 means "free in the model": take it from the model's
    ``z_m_limit``). Provenance ``observed`` (published constraint positions).
    """
    path = Path(path)
    raw = path.read_bytes()
    ids, ras, decs, zs = [], [], [], []
    for n, line in enumerate(raw.decode("latin-1").splitlines(), 1):
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            parts = stripped.lstrip("#").split()
            if parts and parts[0].upper() == "REFERENCE" and len(parts) > 1 and parts[1] != "0":
                raise UnsupportedModelError(f"{path}: relative coordinates (REFERENCE {parts[1]})")
            continue
        parts = stripped.split()
        if len(parts) < 7:
            raise ValueError(f"{path}:{n}: expected 'id ra dec a b theta z ...', got {stripped!r}")
        ids.append(parts[0])
        ras.append(float(parts[1]))
        decs.append(float(parts[2]))
        zs.append(float(parts[6]))
    t = Table(
        {
            "image_id": ids,
            "system": [system_key(i.rsplit(".", 1)[0]) if "." in i else system_key(i) for i in ids],
            "ra": ras,
            "dec": decs,
            "z": zs,
        }
    )
    t.meta.update(
        provenance=schema.Provenance.OBSERVED.value,
        source=str(path),
        sha256=hashlib.sha256(raw).hexdigest(),
    )
    return t


def backtrace_images(model: LensModel, images: Table, z_m_limit: dict[str, float]) -> Table:
    """Trace each image to the source plane and measure each system's scatter.

    Per image: the source-plane position, its offset from the system's mean source position
    (``dbeta_arcsec``), and that offset mapped to the image plane with the local inverse
    magnification matrix (``dtheta_arcsec``, the usual source-plane approximation of the
    image-plane residual). Images whose redshift is neither catalogued nor in ``z_m_limit``
    get NaN (:func:`image_redshifts`). Provenance ``model_prediction``.
    """
    z = image_redshifts(images, z_m_limit)
    ok = np.isfinite(z)
    pred = model.evaluate(images["ra"][ok], images["dec"][ok], z[ok])
    bx = np.full(len(images), np.nan)
    by = np.full(len(images), np.nan)
    bx[ok] = pred["x"] - pred["alpha_x"]
    by[ok] = pred["y"] - pred["alpha_y"]
    # Inverse magnification matrix A = I - H at each image (H from the scaled Hessian).
    a11 = np.full(len(images), np.nan)
    a12 = np.full(len(images), np.nan)
    a22 = np.full(len(images), np.nan)
    kap, g1, g2 = (np.asarray(pred[c], float) for c in ("kappa", "gamma1", "gamma2"))
    a11[ok] = 1.0 - kap - g1  # psi_xx = kappa + gamma1, psi_yy = kappa - gamma1, psi_xy = gamma2
    a12[ok] = -g2
    a22[ok] = 1.0 - kap + g1
    systems = np.asarray(images["system"]).astype(str)
    dbx = np.full(len(images), np.nan)
    dby = np.full(len(images), np.nan)
    for sys_id in np.unique(systems):
        m = (systems == sys_id) & ok
        if m.sum() < 2:
            continue
        dbx[m] = bx[m] - bx[m].mean()
        dby[m] = by[m] - by[m].mean()
    det = a11 * a22 - a12 * a12
    # theta offset = A^-1 (beta offset), with A^-1 = [[a22, -a12], [-a12, a11]] / det.
    dtx = (a22 * dbx - a12 * dby) / det
    dty = (-a12 * dbx + a11 * dby) / det
    out = Table(
        {
            "image_id": images["image_id"],
            "system": systems,
            "ra": images["ra"],
            "dec": images["dec"],
            "z_used": z,
            "beta_x": bx,
            "beta_y": by,
            "dbeta_arcsec": np.hypot(dbx, dby),
            "magnification": 1.0 / det,
            "dtheta_arcsec": np.hypot(dtx, dty),
        }
    )
    out.meta.update(model._meta())
    out.meta["source"] = f"{images.meta.get('source', 'images')} traced with {model.source}"
    return out


@dataclass(frozen=True)
class DeflectionGrid:
    """A model's deflection at D_LS/D_S = 1 on a square model-frame grid (arcsec).

    Deflection scales linearly with D_LS/D_S, so one grid serves every source redshift.
    """

    x: np.ndarray  # 1-D grid coordinates (arcsec), shared by both axes
    alpha_x: np.ndarray  # shape (len(x), len(x)), indexed [iy, ix]
    alpha_y: np.ndarray
    model_sha256: str

    @classmethod
    def compute(
        cls, model: LensModel, half_width: float = 60.0, step: float = 0.1
    ) -> DeflectionGrid:
        g = np.arange(-half_width, half_width + step / 2, step)
        xx, yy = np.meshgrid(g, g)
        ax = np.zeros_like(xx)
        ay = np.zeros_like(xx)
        for comp in model.components:
            dx, dy = comp.deflection(xx, yy)
            ax += dx
            ay += dy
        return cls(g, ax, ay, model.sha256)

    @classmethod
    def cached(
        cls, model: LensModel, path: str | Path, half_width: float = 60.0, step: float = 0.1
    ) -> DeflectionGrid:
        """Load ``path`` (``.npz``) if it was computed for this model and grid, else compute it."""
        path = Path(path)
        if path.exists():
            d = np.load(path)
            g = d["x"]
            if (
                str(d["model_sha256"]) == model.sha256
                and np.isclose(g[-1], half_width)
                and np.isclose(g[1] - g[0], step)
            ):
                return cls(g, d["alpha_x"], d["alpha_y"], model.sha256)
        grid = cls.compute(model, half_width, step)
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            path, x=grid.x, alpha_x=grid.alpha_x, alpha_y=grid.alpha_y, model_sha256=model.sha256
        )
        return grid


def find_images(
    model: LensModel,
    grid: DeflectionGrid,
    beta_x: float,
    beta_y: float,
    z_s: float,
    newton_steps: int = 12,
    tol_arcsec: float = 1e-5,
) -> Table:
    """All image positions of a point source at ``(beta_x, beta_y)`` (model frame, arcsec).

    Every grid cell is split into two triangles and mapped to the source plane; a triangle that
    contains the source seeds a Newton iteration on the lens equation with the analytic model.
    Converged solutions closer than 0.05" are merged. Columns: ``x, y`` (arcsec), ``ra, dec``,
    ``magnification`` (signed; negative = odd parity), ``residual_arcsec`` (source-plane misfit).
    Images outside the grid, or in cells where the map folds below the grid scale, can be missed.
    Provenance ``model_prediction``.
    """
    s = float(model.dls_ds(z_s))
    g = grid.x
    bx = g[None, :] - s * grid.alpha_x
    by = g[:, None] - s * grid.alpha_y
    seeds = []
    n = len(g) - 1
    for (a0, a1), (b0, b1), (c0, c1) in (
        ((0, 0), (0, 1), (1, 0)),
        ((1, 1), (1, 0), (0, 1)),
    ):
        ax_, ay_ = bx[a0 : a0 + n, a1 : a1 + n], by[a0 : a0 + n, a1 : a1 + n]
        bx_, by_ = bx[b0 : b0 + n, b1 : b1 + n], by[b0 : b0 + n, b1 : b1 + n]
        cx_, cy_ = bx[c0 : c0 + n, c1 : c1 + n], by[c0 : c0 + n, c1 : c1 + n]
        det = (by_ - cy_) * (ax_ - cx_) + (cx_ - bx_) * (ay_ - cy_)
        with np.errstate(divide="ignore", invalid="ignore"):
            l1 = ((by_ - cy_) * (beta_x - cx_) + (cx_ - bx_) * (beta_y - cy_)) / det
            l2 = ((cy_ - ay_) * (beta_x - cx_) + (ax_ - cx_) * (beta_y - cy_)) / det
        hit = (l1 >= 0) & (l2 >= 0) & (1 - l1 - l2 >= 0)
        for iy, ix in zip(*np.nonzero(hit), strict=True):
            w1, w2 = l1[iy, ix], l2[iy, ix]
            w3 = 1 - w1 - w2
            seeds.append(
                (
                    w1 * g[ix + a1] + w2 * g[ix + b1] + w3 * g[ix + c1],
                    w1 * g[iy + a0] + w2 * g[iy + b0] + w3 * g[iy + c0],
                )
            )
    found: list[tuple[float, float, float, float]] = []
    for x0, y0 in seeds:
        x, y = np.array([x0]), np.array([y0])
        for _ in range(newton_steps):
            f = model.fields_xy(x, y)
            rx = beta_x - (x - s * f["alpha_x"])
            ry = beta_y - (y - s * f["alpha_y"])
            a11, a12, a22 = 1 - s * f["psi_xx"], -s * f["psi_xy"], 1 - s * f["psi_yy"]
            det = a11 * a22 - a12 * a12
            x = x + (a22 * rx - a12 * ry) / det
            y = y + (-a12 * rx + a11 * ry) / det
        f = model.fields_xy(x, y)
        res = float(np.hypot(beta_x - (x - s * f["alpha_x"]), beta_y - (y - s * f["alpha_y"]))[0])
        if not res < tol_arcsec:
            continue
        a11, a12, a22 = 1 - s * f["psi_xx"], -s * f["psi_xy"], 1 - s * f["psi_yy"]
        mu = float(1.0 / (a11 * a22 - a12 * a12)[0])
        if all(np.hypot(x[0] - u_, y[0] - v_) > 0.05 for u_, v_, _, _ in found):
            found.append((float(x[0]), float(y[0]), mu, res))
    found.sort(key=lambda r: -abs(r[2]))
    xs = np.array([r[0] for r in found])
    ys = np.array([r[1] for r in found])
    ra, dec = model.to_sky(xs, ys)
    out = Table(
        {
            "x": xs,
            "y": ys,
            "ra": ra,
            "dec": dec,
            "magnification": np.array([r[2] for r in found]),
            "residual_arcsec": np.array([r[3] for r in found]),
        }
    )
    out.meta.update(model._meta())
    out.meta.update(z_s=float(z_s), beta_xy=[float(beta_x), float(beta_y)])
    return out
