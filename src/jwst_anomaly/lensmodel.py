"""Published Lenstool lens models evaluated at sky positions (M3 lens-model stage, D-024).

A Lenstool ``.par`` file (e.g. ``best.par`` of a published model) is parsed with
:func:`parse_lenstool_par` and evaluated by :class:`LensModel` at sky positions for a source
redshift: deflection, convergence, shear, magnification and the expected arc orientation.
Every table returned carries ``meta["provenance"] = "model_prediction"``.

Supported: one lens plane whose potentials are all dPIE/PIEMD (Lenstool ``profile 81``), in
absolute coordinates (``runmode`` ``reference 3 ra dec``), with a flat LambdaCDM cosmology.
Anything else raises :class:`UnsupportedModelError`. :meth:`LensModel.split_planes` builds a
:class:`MultiPlaneLensModel` from a parsed model by moving named potentials to their own
redshifts (e.g. a foreground galaxy fitted as a cluster member, D-046).
Scaling-relation galaxy catalogs (``potfile``) are not expanded: use the ``best.par`` of a run,
which lists every potential.

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

import copy
import hashlib
import io
import re
from collections.abc import Sequence
from dataclasses import dataclass, replace
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
        f"{_MAHLER22_RAW}/ICLv2/tmp_k/0000_k.fits.tar.xz",
        "4de5d733587dd4dea12233c25069d8e1307936fcd50de316cc4eb967a53edc2c",
    ),
}

#: Caminha et al. (2023, A&A 678, A3) El Gordo Lenstool model, CDS J/A+A/678/A3. Image-plane
#: optimised (Chi2pos 80.22). No convergence map; the best-fit magnification map at z_s = 2 is
#: the map check.
_CAMINHA23_CDS = "https://cdsarc.cds.unistra.fr/ftp/J/A+A/678/A3"
ELGORDO_CAMINHA23: dict[str, tuple[str, str]] = {
    "best.par": (
        f"{_CAMINHA23_CDS}/files/best_fit.par",
        "7b0153ae0ee02f057f6aaa6f46b1b698502e6fc427266ac9a09d241ddc63a472",
    ),
    "arcs.dat": (
        f"{_CAMINHA23_CDS}/files/obs_arcs_v1_new_IDs.dat",
        "d631743921266c34689a1d509f08e53dc3c90bc88064393d7b8fd524a3d5c700",
    ),
}

#: Bergamini et al. (2023b, ApJ 952, 84) Abell 2744 Lenstool model from the authors' page.
#: Image-plane optimised (Chi2pos 146.60). No published FITS maps.
_BERGAMINI23_WEB = "https://www.fe.infn.it/astro/lensing/A2744_Bergamini23"
ABELL2744_BERGAMINI23: dict[str, tuple[str, str]] = {
    "best.par": (
        f"{_BERGAMINI23_WEB}/best.par",
        "7245368f96ad9c7159eb9c8d0045030eda0804554312ee86d84f2e052025b9fb",
    ),
    "arcs.dat": (
        f"{_BERGAMINI23_WEB}/obs_arcs.cat",
        "d02c231f4ee8c81f47335a99182a9f64a4c553e14818b1c9bd07314c2f4f5e1c",
    ),
}

#: MCMC samples of the two models above (Lenstool ``bayes.dat``; D-045). Kept apart from the
#: model file sets because the Abell 2744 chain is 70 MB and only ``posterior`` needs it.
ELGORDO_CAMINHA23_BAYES = (
    f"{_CAMINHA23_CDS}/files/bayes.dat",
    "2d3f736218c1dbba56fdd39a7f8151d307c46bd37abf39bac8831a9e6a81de44",
)
ABELL2744_BERGAMINI23_BAYES = (
    f"{_BERGAMINI23_WEB}/bayes.dat",
    "bf6ae6702a021d185b70305a83d45f2fca7dacf2a60138f815f5d69ff7524c1b",
)

# CANUCS DR1 Lenstool best fits (doi:10.17909/18nv-np70): MACS0416 by Rihtarsic et al. 2025
# (A&A, doi:10.1051/0004-6361/202451117; image-plane chi2pos 344.30) and Abell 370 by Gledhill
# et al. 2025 (ApJ, doi:10.3847/1538-4357/ad684a; source-plane fit). ``input.par`` is the
# Lenstool input file (its sigposArcsec); the multiple-image files are their ``multfile``s (D-044).
_CANUCS_MODEL = "https://archive.stsci.edu/hlsps/canucs/dr1/{0}/model/hlsp_canucs_jwst-hst_multi_{0}-{1}_multi_v1_model.txt"
MACS0416_CANUCS: dict[str, tuple[str, str]] = {
    "best.par": (
        _CANUCS_MODEL.format("macs0416", "lenstool-bestparam"),
        "f3d9a8415044ff8d5ab4774573fcf3670c9dee9169abfca46086bdbf3dd45bed",
    ),
    "arcs.dat": (
        _CANUCS_MODEL.format("macs0416", "lenstool-multim"),
        "ce00444dcc9239f0fd72d5fb37e35cbbea281647803ef4c30eb7b187262f5507",
    ),
    "input.par": (
        _CANUCS_MODEL.format("macs0416", "lenstool-param"),
        "0f1fb7d6947d467b28d8b74485799321b09a1cc63fd7e7df82f01620ddb9b337",
    ),
}
ABELL370_CANUCS: dict[str, tuple[str, str]] = {
    "best.par": (
        _CANUCS_MODEL.format("a370", "lenstool-bestparam"),
        "3c1eea91755ea532e424b1b143e39b35a9d89beae2958d3b042876af50443580",
    ),
    "arcs.dat": (
        _CANUCS_MODEL.format("a370", "lenstool-multim"),
        "d72c3d98e675e5bc00cfbdd9b84d1b8528b22e36311924fea072293af23ef9e2",
    ),
    "input.par": (
        _CANUCS_MODEL.format("a370", "lenstool-param"),
        "aacc2dadd442645d2222c23ee2c3f9f6a76fddaa73996a26a6e88c0691a6bf3d",
    ),
}

#: RELICS (Cerny et al. 2018, ApJ 859, 159; HLSP DOI 10.17909/T9SP45) Lenstool v1 maps of
#: WHL0137-08 (Sunrise), deflection in arcsec at D_LS/D_S = 1 (accessed 2026-10-08). The lens
#: redshift 0.566 and H0 = 70, Om0 = 0.3 reproduce the published z = 6.2 magnification map.
_RELICS_WHL0137 = (
    "https://archive.stsci.edu/hlsps/relics/whl0137m08/models/lenstool/v1/"
    "hlsp_relics_model_model_whl0137m08_lenstool_v1"
)
WHL0137_RELICS_LENSTOOL: dict[str, tuple[str, str]] = {
    "alpha_x": (
        f"{_RELICS_WHL0137}_x-arcsec-deflect.fits",
        "9666a25c06f6f5f078d60ab7093f79948d6fb0adbee548c36fb761dac35f9694",
    ),
    "alpha_y": (
        f"{_RELICS_WHL0137}_y-arcsec-deflect.fits",
        "85aa9f058850b983bdac7566dcbc090d2a88beeae2b251005c97912a64f78c0b",
    ),
    "kappa_map": (
        f"{_RELICS_WHL0137}_kappa.fits",
        "cce65a0723cd497406f99e25365e606d11d9147f370810ab95487a02dc03ca86",
    ),
    "mag_z6.2": (
        f"{_RELICS_WHL0137}_z06p2-magnif.fits",
        "3598be3c0caacb259dfe45fb49b80e5c05e1414f5c212d2a43f01eea028c1be0",
    ),
}

#: HFF CATS Lenstool models (Jauzac, Richard, Mahler et al.; HLSP
#: https://archive.stsci.edu/prepds/frontier/lensmodels/, accessed 2026-10-08): deflection maps
#: (arcsec, D_LS/D_S = 1), κ, z = 2 magnification and the multiple-image list. cluster ->
#: (version, z_lens, {name: (url, sha256)}). The lens redshifts reproduce the published z = 2 maps.
_HFF = "https://archive.stsci.edu/pub/hlsp/frontier"
HFF_CATS: dict[str, tuple[str, float, dict[str, tuple[str, str]]]] = {
    "macs0416": (
        "v4.1",
        0.396,
        {
            "alpha_x": (
                f"{_HFF}/macs0416/models/cats/v4.1/hlsp_frontier_model_macs0416_cats_v4.1_x-arcsec-deflect.fits",
                "538789ef75b2660f880b573be0d45091599d85307bb24300c90453fa81333e91",
            ),
            "alpha_y": (
                f"{_HFF}/macs0416/models/cats/v4.1/hlsp_frontier_model_macs0416_cats_v4.1_y-arcsec-deflect.fits",
                "4fdef438de6384ea87c73a7a59ff5ade737055877aa7089b74b7a60982713d13",
            ),
            "kappa_map": (
                f"{_HFF}/macs0416/models/cats/v4.1/hlsp_frontier_model_macs0416_cats_v4.1_kappa.fits",
                "890ebc82c47e38e7154db5221813bd64d48c9ba7baa55f5e2441e63ba1368296",
            ),
            "mag_z2": (
                f"{_HFF}/macs0416/models/cats/v4.1/hlsp_frontier_model_macs0416_cats_v4.1_z02-magnif.fits",
                "cd01824f7317fdfd019921538380fa434e7d62646a7751f96f799a1524299a82",
            ),
            "arcs.dat": (
                f"{_HFF}/macs0416/models/cats/v4.1/hlsp_frontier_model_macs0416_cats_v4.1_arcs.txt",
                "59b8b56527f80b2ca5c1b72ab36f608b78143f1d24b86c081d9b50f471f80c15",
            ),
            "params.txt": (
                f"{_HFF}/macs0416/models/cats/v4.1/hlsp_frontier_model_macs0416_cats_v4.1_params.txt",
                "7149d8a55792c1b42aae9badb25ff9b0934cfbc780ca37b7ef1f0cb21d71f996",
            ),
        },
    ),
    "macs1149": (
        "v4.1",
        0.543,
        {
            "alpha_x": (
                f"{_HFF}/macs1149/models/cats/v4.1/hlsp_frontier_model_macs1149_cats_v4.1_x-arcsec-deflect.fits",
                "a39dda08daefce24000756fb30b0483a466dd9c1f2dda23fc5d6acea32603c7f",
            ),
            "alpha_y": (
                f"{_HFF}/macs1149/models/cats/v4.1/hlsp_frontier_model_macs1149_cats_v4.1_y-arcsec-deflect.fits",
                "5b3257f22d6b4f373793f7a8924f9768286e710d6f596425109532d5335aa57a",
            ),
            "kappa_map": (
                f"{_HFF}/macs1149/models/cats/v4.1/hlsp_frontier_model_macs1149_cats_v4.1_kappa.fits",
                "37005861ec7a7af7fd137fd74c1c5fe382bbf28998f721b0b779c028145a1257",
            ),
            "mag_z2": (
                f"{_HFF}/macs1149/models/cats/v4.1/hlsp_frontier_model_macs1149_cats_v4.1_z02-magnif.fits",
                "6b750d2d18aaa6faaa5b3a25e7a11f5b62adaf0e0fa8a5940d7c6224b3251397",
            ),
            "arcs.dat": (
                f"{_HFF}/macs1149/models/cats/v4.1/hlsp_frontier_model_macs1149_cats_v4.1_arcs.txt",
                "1a12c542fada53e0ba56c29ed73012c41331ce8f149357d6f5c2a10b8d808cf5",
            ),
            "params.txt": (
                f"{_HFF}/macs1149/models/cats/v4.1/hlsp_frontier_model_macs1149_cats_v4.1_params.txt",
                "37ccf2eccbd47ce3985a1907d1764fab775e0f20ca041fdc1249432cd81e7541",
            ),
        },
    ),
    "abell370": (
        "v4",
        0.375,
        {
            "alpha_x": (
                f"{_HFF}/abell370/models/cats/v4/hlsp_frontier_model_abell370_cats_v4_x-arcsec-deflect.fits",
                "1a4c7ee4f1c0e1330cfd4011587add2dd619680c0f30104668f5b8ef209a2c3d",
            ),
            "alpha_y": (
                f"{_HFF}/abell370/models/cats/v4/hlsp_frontier_model_abell370_cats_v4_y-arcsec-deflect.fits",
                "94153ce7850b2d02505df066ede18641b6e95143416a99a3469de2f41309020c",
            ),
            "kappa_map": (
                f"{_HFF}/abell370/models/cats/v4/hlsp_frontier_model_abell370_cats_v4_kappa.fits",
                "5eee21d023e80a04741cd54116272800211c720ea450dd7cc9ea121905ad915d",
            ),
            "mag_z2": (
                f"{_HFF}/abell370/models/cats/v4/hlsp_frontier_model_abell370_cats_v4_z02-magnif.fits",
                "8773f6044147f46c06dd274faf9e3f0f452b15e450a6260557929dd095276741",
            ),
            "arcs.dat": (
                f"{_HFF}/abell370/models/cats/v4/hlsp_frontier_model_abell370_cats_v4_arcs.txt",
                "27cd08c5b0087b841e5b1775fe5f64e7db70f0f3cbd8e4c98d9e5af416f13052",
            ),
        },
    ),
    "macs0717": (
        "v4.1",
        0.545,
        {
            "alpha_x": (
                f"{_HFF}/macs0717/models/cats/v4.1/hlsp_frontier_model_macs0717_cats_v4.1_x-arcsec-deflect.fits",
                "4916865e8513eed33cea613fb96007e2d87535ba44359ced6babaee4959ac5ad",
            ),
            "alpha_y": (
                f"{_HFF}/macs0717/models/cats/v4.1/hlsp_frontier_model_macs0717_cats_v4.1_y-arcsec-deflect.fits",
                "f3e69604862a51211064d1fb29777dbf2bfac1eea48f7f7430eefd90a195b458",
            ),
            "kappa_map": (
                f"{_HFF}/macs0717/models/cats/v4.1/hlsp_frontier_model_macs0717_cats_v4.1_kappa.fits",
                "6919ef95fd876a4a3d6a43a7cd1b7c35511533d3e7433fa986b7d94fca7af579",
            ),
            "mag_z2": (
                f"{_HFF}/macs0717/models/cats/v4.1/hlsp_frontier_model_macs0717_cats_v4.1_z02-magnif.fits",
                "f0ada68b3db2cbe933af63af697b945b1da37d3d929be87f9ceed0adb8b8d6d8",
            ),
            "arcs.dat": (
                f"{_HFF}/macs0717/models/cats/v4.1/hlsp_frontier_model_macs0717_cats_v4.1_arcs.txt",
                "abc9f95f89467160a59d87e7252eea2667c97ff292d25d67880ce8b1565ad73e",
            ),
            "params.txt": (
                f"{_HFF}/macs0717/models/cats/v4.1/hlsp_frontier_model_macs0717_cats_v4.1_params.txt",
                "7ea3af8dbb3e6dd9917f9155a17a24903d0337ddd6fb63ff17522ca262ed22e4",
            ),
        },
    ),
    "abells1063": (
        "v4.1",
        0.348,
        {
            "alpha_x": (
                f"{_HFF}/abells1063/models/cats/v4.1/hlsp_frontier_model_abells1063_cats_v4.1_x-arcsec-deflect.fits",
                "7aadd23ef613099a522319962b1c9f98606d54c8d895adf9fe715dd9a7b460a7",
            ),
            "alpha_y": (
                f"{_HFF}/abells1063/models/cats/v4.1/hlsp_frontier_model_abells1063_cats_v4.1_y-arcsec-deflect.fits",
                "5b9a81cd1741b13fc09c5a164dc110a01b45ae50915fe41db4b9a013d769dbe7",
            ),
            "kappa_map": (
                f"{_HFF}/abells1063/models/cats/v4.1/hlsp_frontier_model_abells1063_cats_v4.1_kappa.fits",
                "8bae5b69c2a6667b10a3c2fab065dccc8a87192b715be36ba4e869a39c537fe0",
            ),
            "mag_z2": (
                f"{_HFF}/abells1063/models/cats/v4.1/hlsp_frontier_model_abells1063_cats_v4.1_z02-magnif.fits",
                "2ae0c461add316de449d95078725bd37e64af5a913105141c96e3a77331e1ddb",
            ),
            "arcs.dat": (
                f"{_HFF}/abells1063/models/cats/v4.1/hlsp_frontier_model_abells1063_cats_v4.1_arcs.txt",
                "83b207e1b0085678f1b5b9768347f94f352a5390800aa5d11de80038cda4e410",
            ),
        },
    ),
    "abell2744": (
        "v4.1",
        0.308,
        {
            "alpha_x": (
                f"{_HFF}/abell2744/models/cats/v4.1/hlsp_frontier_model_abell2744_cats_v4.1_x-arcsec-deflect.fits",
                "931b6603f31fc32368cc1f41098ec4cbceffdf910ca99231aa323cfcb4b08e2c",
            ),
            "alpha_y": (
                f"{_HFF}/abell2744/models/cats/v4.1/hlsp_frontier_model_abell2744_cats_v4.1_y-arcsec-deflect.fits",
                "4a131c57085384b00fe56edbb2c16027734bfd8cd5b8ba57e2f6439c5c2843a3",
            ),
            "kappa_map": (
                f"{_HFF}/abell2744/models/cats/v4.1/hlsp_frontier_model_abell2744_cats_v4.1_kappa.fits",
                "bd1cf25162bd71ced4f8b75d063713b12d4054e83e409e8fd48e98205597d47b",
            ),
            "mag_z2": (
                f"{_HFF}/abell2744/models/cats/v4.1/hlsp_frontier_model_abell2744_cats_v4.1_z02-magnif.fits",
                "d74639912a9b1abc488d6ae5cbe1f8336f53b3e34d138260c65de20e39008132",
            ),
            "arcs.dat": (
                f"{_HFF}/abell2744/models/cats/v4.1/hlsp_frontier_model_abell2744_cats_v4.1_arcs.txt",
                "2e16e2884f4d1a386ae305737be4d97880e6cb3954117f4efbd62fc61eb28b19",
            ),
        },
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


_LETTER_SUFFIX = re.compile(r"^(.*\d)([A-Za-z]+)$")


def image_family(image_id: str | float) -> str:
    """Multiple-image family (one point source) of a Lenstool image id.

    Trailing letters mark the images of one family (``"23a"`` -> ``"23"``, ``"1.1a"`` ->
    ``"1.1"``, ``"A200.1a"`` -> ``"A200.1"``); otherwise the last ``.N`` does (``"4.1"`` ->
    ``"4"``). Ids without either are a family name already (``"4"``, ``"4.0"`` -> ``"4"``, as in
    ``z_m_limit``).
    """
    text = str(image_id).strip()
    match = _LETTER_SUFFIX.match(text)
    if match:
        return system_key(match.group(1))
    if "." in text and not _INTEGER_ID.match(text):
        return system_key(text.rsplit(".", 1)[0])
    return system_key(text)


def _z_m_limit_entry(vals: list[str], where: str) -> dict[str, float]:
    """Fixed redshifts of one ``z_m_limit`` line (values after the keyword), by image family.

    ``z_m_limit <n> <image id>... <flag> <z or zmin> <zmax> <step>``; flag 0 = fixed. Several
    ids may share one redshift (Bergamini+2023b "A200.1a B200.2a"). Flags 1-4 and -n
    (parabolic) are free redshifts and are skipped."""
    if len(vals) < 4:
        raise ValueError(f"{where}: malformed z_m_limit {vals}")
    names, flag, z = (
        (vals[1:-4], vals[-4], vals[-3]) if len(vals) >= 6 else (vals[1:2], vals[2], vals[3])
    )
    if not re.fullmatch(r"-?\d+", flag):
        raise ValueError(f"{where}: malformed z_m_limit {vals}")
    return {image_family(name): float(z) for name in names} if int(flag) == 0 else {}


def read_z_m_limit(path: str | Path) -> dict[str, float]:
    """Only the fixed ``z_m_limit`` redshifts of a Lenstool parameter file, without parsing (or
    supporting) the rest of it: for map models whose published ``params.txt`` uses features
    this module does not implement (e.g. ``potfile``). Same rules as :func:`parse_lenstool_par`."""
    out: dict[str, float] = {}
    for n, line in enumerate(Path(path).read_text(encoding="latin-1").splitlines(), 1):
        vals = line.split("#", 1)[0].split()
        if vals and vals[0] == "z_m_limit":
            out.update(_z_m_limit_entry(vals[1:], f"{path}:{n}"))
    return out


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
    * ``z_m_limit``: fixed (best-fit) redshifts of multiple-image families, by :func:`image_family`
      of each listed image id (one line may list several ids);
    * ``z_m_limit_groups``: the families of each such line (they share one redshift);
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
        "z_m_limit_groups": [],
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
                    entry = _z_m_limit_entry(vals, f"{path}:{n}")
                    out["z_m_limit"].update(entry)
                    if entry:  # families of one line share one (sampled) redshift
                        out["z_m_limit_groups"].append(sorted(entry))
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
    # ``mag`` marks a potfile (scaling-relation) member; :func:`posterior_par` rescales those.
    pot["mag"] = float(values["mag"][0]) if "mag" in values else float("nan")
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

    def shift_frame(self, dra_arcsec: float, ddec_arcsec: float) -> None:
        """Move the whole model on the sky by (dRA cos dec, dDec) arcsec, in place (D-034).

        Model-frame quantities are unchanged; only the frame's sky anchor moves.
        """
        self.ra0 += dra_arcsec / 3600.0 / self._cos0
        self.dec0 += ddec_arcsec / 3600.0
        self._cos0 = np.cos(np.deg2rad(self.dec0))

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

    def deflection_xy(self, x: Any, y: Any) -> tuple[np.ndarray, np.ndarray]:
        """Deflection at D_LS/D_S = 1 at model-frame positions (arcsec)."""
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        ax, ay = np.zeros(x.shape), np.zeros(x.shape)
        for comp in self.components:
            dx, dy = comp.deflection(x, y)
            ax += dx
            ay += dy
        return ax, ay

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

    # --- lens equation (shared with MultiPlaneLensModel) -------------------------------------
    def lens_map(self, x: Any, y: Any, z_s: Any) -> dict[str, np.ndarray]:
        """Source-plane position and Jacobian ``A = d beta / d theta`` at model-frame positions.

        Keys ``beta_x, beta_y`` (arcsec) and ``a11, a12, a21, a22``; ``z_s`` is a scalar or one
        redshift per position.
        """
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        s = self.dls_ds(z_s)
        f = self.fields_xy(x, y)
        return {
            "beta_x": x - s * f["alpha_x"],
            "beta_y": y - s * f["alpha_y"],
            "a11": 1.0 - s * f["psi_xx"],
            "a12": -s * f["psi_xy"],
            "a21": -s * f["psi_xy"],
            "a22": 1.0 - s * f["psi_yy"],
        }

    def source_points(self, x: Any, y: Any, z_s: Any) -> tuple[np.ndarray, np.ndarray]:
        """Source-plane position only (no Hessian) at model-frame positions."""
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        s = self.dls_ds(z_s)
        ax, ay = self.deflection_xy(x, y)
        return x - s * np.asarray(ax), y - s * np.asarray(ay)

    def source_grid(self, grid: DeflectionGrid, z_s: float) -> tuple[np.ndarray, np.ndarray]:
        """Source-plane position of every node of ``grid`` (indexed ``[iy, ix]``)."""
        if grid.alpha_x.shape != (len(grid.x), len(grid.x)):
            raise ValueError("the deflection grid was not computed for this single-plane model")
        s = float(self.dls_ds(z_s))
        g = grid.x
        return g[None, :] - s * grid.alpha_x, g[:, None] - s * grid.alpha_y

    def split_planes(
        self, z_lens: dict[str, float], v_disp: dict[str, float] | None = None
    ) -> MultiPlaneLensModel:
        """A multi-plane copy with the named potentials moved to their own redshifts (D-046).

        ``z_lens`` maps potential names to their redshift (e.g. a spectroscopic redshift for a
        galaxy that the model treats as a cluster member); ``v_disp`` optionally replaces their
        velocity dispersion (km/s). Potentials at the same redshift share a plane; the rest stay
        on this model's plane.

        A single-plane fit places every potential at its observed (image-plane) position. A
        moved potential with another plane in front of it is put where the ray through its
        observed centre crosses its own plane (delensed through the foreground planes, nearest
        first), so it is still seen where it was fitted. Its radii, ellipticity and angle are kept
        (ASSUMPTION: the foreground's distortion of its shape is neglected), and potentials left
        on this model's plane keep their fitted positions.
        """
        names = [c.name for c in self.components]
        unknown = sorted(set(z_lens) - set(names)) + sorted(set(v_disp or {}) - set(z_lens))
        if unknown:
            raise ValueError(f"split_planes: unknown or unmoved potentials {unknown}")
        by_z: dict[float, list[DPIE]] = {}
        for comp in self.components:
            z = float(z_lens.get(comp.name, comp.z_lens))
            sig = (v_disp or {}).get(comp.name, comp.v_disp)
            by_z.setdefault(z, []).append(replace(comp, z_lens=z, v_disp=float(sig)))

        def ident(comps: list[DPIE]) -> str:  # a plane holds a subset: not the file's hash
            if not self.sha256:  # in memory: no identity, as for the unsplit model
                return ""
            text = f"{self.sha256}{comps!r}{_cosmology_key(self.cosmology)}"
            return hashlib.sha256(text.encode()).hexdigest()

        def build(groups: dict[float, list[DPIE]]) -> list[LensModel]:
            return [
                LensModel(comps, self.ra0, self.dec0, self.cosmology, self.source, ident(comps))
                for comps in groups.values()
            ]

        # delens moved potentials plane by plane, in redshift order, so each one is traced
        # through foreground planes that are already final
        fitted_z = {c.name: c.z_lens for c in self.components}
        for z in sorted(by_z):
            comps = by_z[z]
            # only potentials that changed redshift (a v_disp-only change stays where fitted)
            moved_here = [
                k for k, c in enumerate(comps) if c.name in z_lens and z != fitted_z[c.name]
            ]
            front = [zz for zz in by_z if zz < z]
            if not moved_here or not front:
                continue
            fore = MultiPlaneLensModel(build({zz: by_z[zz] for zz in sorted(front)}))
            # the ray through the fitted centre, traced to this plane as if it were a source
            px, py = fore.source_points(
                np.array([comps[k].x for k in moved_here]),
                np.array([comps[k].y for k in moved_here]),
                z,
            )
            for k, x_new, y_new in zip(moved_here, px, py, strict=True):
                comps[k] = replace(comps[k], x=float(x_new), y=float(y_new))
        planes = build(by_z)
        moved = ", ".join(f"{k}@z{v:g}" for k, v in sorted(z_lens.items()))
        sigs = ", ".join(f"{k}:{v:g} km/s" for k, v in sorted((v_disp or {}).items()))
        return MultiPlaneLensModel(
            planes, source=f"{self.source} [planes: {moved}{'; ' + sigs if sigs else ''}]"
        )

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
        # 2 % plus 1e-6": best.par writes radii with 6 decimals, so a 2e-5" core is coarse.
        tol = 0.02 * val + 1e-6
        if np.isfinite(val) and np.isfinite(kpc) and abs(kpc / per_arcsec - val) > tol:
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
    """Lenstool ``arcs.dat`` (absolute coordinates): ``image_id, system, ra, dec, a, z``.

    ``system`` is the image's family (:func:`image_family`). ``a`` is the file's first shape
    column, which image lists for image-plane models use as the position error (arcsec). ``z`` is
    the catalogued redshift (0 means "free in the model": take it from the model's
    ``z_m_limit``). Provenance ``observed`` (published constraint positions).
    """
    path = Path(path)
    raw = path.read_bytes()
    ids, ras, decs, errs, zs = [], [], [], [], []
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
        errs.append(float(parts[3]))
        zs.append(float(parts[6].strip("()")))  # "(2.16)": a redshift the model fitted
    t = Table(
        {
            "image_id": ids,
            "system": [image_family(i) for i in ids],
            "ra": ras,
            "dec": decs,
            "a": errs,
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
    x, y = model.to_frame(np.asarray(images["ra"], float)[ok], np.asarray(images["dec"], float)[ok])
    lm = model.lens_map(x, y, z[ok])
    bx = np.full(len(images), np.nan)
    by = np.full(len(images), np.nan)
    bx[ok], by[ok] = lm["beta_x"], lm["beta_y"]
    # Inverse magnification matrix A = d beta / d theta at each image (I - H on one plane).
    a11, a12, a21, a22 = (np.full(len(images), np.nan) for _ in range(4))
    a11[ok], a12[ok], a21[ok], a22[ok] = lm["a11"], lm["a12"], lm["a21"], lm["a22"]
    systems = np.asarray(images["system"]).astype(str)
    dbx = np.full(len(images), np.nan)
    dby = np.full(len(images), np.nan)
    for sys_id in np.unique(systems):
        m = (systems == sys_id) & ok
        if m.sum() < 2:
            continue
        dbx[m] = bx[m] - bx[m].mean()
        dby[m] = by[m] - by[m].mean()
    det = a11 * a22 - a12 * a21
    # theta offset = A^-1 (beta offset), with A^-1 = [[a22, -a12], [-a21, a11]] / det.
    dtx = (a22 * dbx - a12 * dby) / det
    dty = (-a21 * dbx + a11 * dby) / det
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

    Deflection scales linearly with D_LS/D_S, so one grid serves every source redshift. For a
    :class:`MultiPlaneLensModel` the arrays hold each plane's ``alpha_i(theta_i)`` along the ray
    through every node, shape ``(n_planes, len(x), len(x))`` (the rays do not depend on z_s).
    Use a grid only with the model it was computed for (``model.source_grid`` checks the shape).
    """

    x: np.ndarray  # 1-D grid coordinates (arcsec), shared by both axes
    alpha_x: np.ndarray  # shape (len(x), len(x)) [iy, ix]; (n_planes, len(x), len(x)) multi-plane
    alpha_y: np.ndarray
    model_sha256: str

    @classmethod
    def compute(
        cls, model: LensModel, half_width: float = 60.0, step: float = 0.1
    ) -> DeflectionGrid:
        g = np.arange(-half_width, half_width + step / 2, step)
        xx, yy = np.meshgrid(g, g)
        if isinstance(model, MultiPlaneLensModel):  # one deflection per plane along each ray
            ax, ay = model.plane_deflections(xx, yy)
        else:
            ax, ay = model.deflection_xy(xx, yy)
        return cls(g, ax, ay, model.sha256)

    @classmethod
    def cached(
        cls, model: LensModel, path: str | Path, half_width: float = 60.0, step: float = 0.1
    ) -> DeflectionGrid:
        """Load ``path`` (``.npz``) if it was computed for this model and grid, else compute it."""
        path = Path(path)
        if not model.sha256:  # an in-memory model has no identity to key the cache on
            return cls.compute(model, half_width, step)
        if path.exists():
            with np.load(path) as d:
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


def _triangle_seeds(
    bx: np.ndarray, by: np.ndarray, gx: np.ndarray, gy: np.ndarray, beta_x: float, beta_y: float
) -> tuple[np.ndarray, np.ndarray]:
    """Seeds from the triangles of quadrilateral meshes that contain ``(beta_x, beta_y)``.

    ``bx, by`` (mapped source-plane coordinates) and ``gx, gy`` (image-plane coordinates) have
    shape ``(..., n, m)``: one or more meshes of ``n × m`` nodes. Every mesh cell is split into two
    triangles; a triangle that contains beta gives the barycentric image-plane point as a seed.
    """

    # Pre-filter: a triangle can contain the source only if its cell's mapped corners bracket
    # beta in both coordinates (inclusive, so boundary hits of the barycentric test are kept).
    def brackets(m: np.ndarray, value: float) -> np.ndarray:
        corners = (m[..., :-1, :-1], m[..., :-1, 1:], m[..., 1:, :-1], m[..., 1:, 1:])
        return (np.minimum.reduce(corners) <= value) & (value <= np.maximum.reduce(corners))

    idx = np.nonzero(brackets(bx, beta_x) & brackets(by, beta_y))
    lead, rows, cols = idx[:-2], idx[-2], idx[-1]
    seeds_x: list[np.ndarray] = []
    seeds_y: list[np.ndarray] = []
    for (a0, a1), (b0, b1), (c0, c1) in (
        ((0, 0), (0, 1), (1, 0)),
        ((1, 1), (1, 0), (0, 1)),
    ):
        ka, kb, kc = (
            (lead + (rows + a0, cols + a1)),
            (lead + (rows + b0, cols + b1)),
            (lead + (rows + c0, cols + c1)),
        )
        pax, pay = bx[ka], by[ka]
        pbx, pby = bx[kb], by[kb]
        pcx, pcy = bx[kc], by[kc]
        det = (pby - pcy) * (pax - pcx) + (pcx - pbx) * (pay - pcy)
        with np.errstate(divide="ignore", invalid="ignore"):
            l1 = ((pby - pcy) * (beta_x - pcx) + (pcx - pbx) * (beta_y - pcy)) / det
            l2 = ((pcy - pay) * (beta_x - pcx) + (pax - pcx) * (beta_y - pcy)) / det
        hit = (l1 >= 0) & (l2 >= 0) & (1 - l1 - l2 >= 0)
        w1, w2 = l1[hit], l2[hit]
        w3 = 1 - w1 - w2
        seeds_x.append(w1 * gx[ka][hit] + w2 * gx[kb][hit] + w3 * gx[kc][hit])
        seeds_y.append(w1 * gy[ka][hit] + w2 * gy[kb][hit] + w3 * gy[kc][hit])
    return np.concatenate(seeds_x), np.concatenate(seeds_y)


def _fold_cells(
    bx: np.ndarray, by: np.ndarray, beta_x: float, beta_y: float
) -> tuple[np.ndarray, np.ndarray]:
    """Grid cells next to a critical curve whose source-plane footprint may hide an image of beta.

    A critical curve is where the mapped triangles change orientation (the sign of their signed
    area). Near a fold, a merging image pair can sit inside one cell even when the cell's mapped
    corners do not bracket beta, so a cell qualifies when it, or a neighbour, contains an
    orientation change and beta lies within its mapped bounding box widened by the box's own size.
    """

    def widened(m: np.ndarray, value: float) -> np.ndarray:
        corners = (m[:-1, :-1], m[:-1, 1:], m[1:, :-1], m[1:, 1:])
        lo, hi = np.minimum.reduce(corners), np.maximum.reduce(corners)
        pad = hi - lo
        return (lo - pad <= value) & (value <= hi + pad)

    cand = widened(bx, beta_x) & widened(by, beta_y)
    if not cand.any():
        return np.nonzero(cand)

    def dilate(mask: np.ndarray) -> np.ndarray:
        out = mask.copy()
        out[1:, :] |= mask[:-1, :]
        out[:-1, :] |= mask[1:, :]
        out[:, 1:] |= mask[:, :-1]
        out[:, :-1] |= mask[:, 1:]
        return out

    # Orientation signs only where the test below can look (two cells around the candidates).
    band = dilate(dilate(cand))
    rr, cc = np.nonzero(band)
    x00, x01, x10, x11 = bx[rr, cc], bx[rr, cc + 1], bx[rr + 1, cc], bx[rr + 1, cc + 1]
    y00, y01, y10, y11 = by[rr, cc], by[rr, cc + 1], by[rr + 1, cc], by[rr + 1, cc + 1]
    sa = np.zeros(cand.shape, np.int8)
    sb = np.zeros(cand.shape, np.int8)
    with np.errstate(invalid="ignore"):  # map models are NaN outside their coverage: sign 0
        sa[rr, cc] = np.nan_to_num(np.sign((x01 - x00) * (y10 - y00) - (x10 - x00) * (y01 - y00)))
        sb[rr, cc] = np.nan_to_num(np.sign((x10 - x11) * (y01 - y11) - (x01 - x11) * (y10 - y11)))
    crit = sa != sb
    crit[:, 1:] |= sa[:, 1:] != sa[:, :-1]
    crit[1:, :] |= sa[1:, :] != sa[:-1, :]
    # dilate by one cell so both sides of the curve are refined; signs outside the two-cell
    # band are unset, but only cells within one cell of a candidate are kept
    return np.nonzero(dilate(crit & band) & cand)


def find_images(
    model: LensModel,
    grid: DeflectionGrid,
    beta_x: float,
    beta_y: float,
    z_s: float,
    newton_steps: int = 12,
    tol_arcsec: float = 1e-5,
    refine_arcsec: float = 0.02,
) -> Table:
    """All image positions of a point source at ``(beta_x, beta_y)`` (model frame, arcsec).

    Every grid cell is split into two triangles and mapped to the source plane; a triangle that
    contains the source seeds a Newton iteration on the lens equation with the analytic model.
    Cells next to a critical curve whose source-plane footprint is near the source are subdivided
    into sub-cells of at most ``refine_arcsec`` (below the 0.05" merge radius) with the model
    evaluated directly, so a merging pair inside one grid cell is still found (D-040);
    ``refine_arcsec <= 0`` disables this. Converged solutions closer than 0.05" are merged.
    Columns: ``x, y`` (arcsec), ``ra, dec``, ``magnification`` (signed; negative = odd parity),
    ``residual_arcsec`` (source-plane misfit). Images outside the grid can be missed.
    Provenance ``model_prediction``.
    """
    if grid.model_sha256 and model.sha256 and grid.model_sha256 != model.sha256:
        raise ValueError("find_images: the deflection grid was computed for another model")
    g = grid.x
    bx, by = model.source_grid(grid, z_s)
    gx, gy = np.broadcast_to(g[None, :], bx.shape), np.broadcast_to(g[:, None], by.shape)
    x0, y0 = _triangle_seeds(bx, by, gx, gy, beta_x, beta_y)
    seeds_x, seeds_y = [x0], [y0]
    step = g[1] - g[0]
    n_sub = int(np.ceil(step / refine_arcsec - 1e-9)) if refine_arcsec > 0 else 1
    if n_sub > 1:
        rows, cols = _fold_cells(bx, by, beta_x, beta_y)
        if len(rows):
            t = np.linspace(0.0, step, n_sub + 1)
            fx = g[cols][:, None, None] + t[None, None, :]
            fy = g[rows][:, None, None] + t[None, :, None]
            fx, fy = np.broadcast_arrays(fx, fy)
            fbx, fby = model.source_points(fx.ravel(), fy.ravel(), z_s)
            fbx, fby = fbx.reshape(fx.shape), fby.reshape(fy.shape)
            x1, y1 = _triangle_seeds(fbx, fby, fx, fy, beta_x, beta_y)
            seeds_x.append(x1)
            seeds_y.append(y1)
    # Newton steps for every seed at once: one lens_map call per step over all potentials.
    x, y = np.concatenate(seeds_x), np.concatenate(seeds_y)
    for _ in range(newton_steps):
        m = model.lens_map(x, y, z_s)
        rx, ry = beta_x - m["beta_x"], beta_y - m["beta_y"]
        det = m["a11"] * m["a22"] - m["a12"] * m["a21"]
        x = x + (m["a22"] * rx - m["a12"] * ry) / det  # a seed on a critical curve diverges and
        y = y + (-m["a21"] * rx + m["a11"] * ry) / det  # is rejected below
    m = model.lens_map(x, y, z_s)
    res = np.hypot(beta_x - m["beta_x"], beta_y - m["beta_y"])
    det = m["a11"] * m["a22"] - m["a12"] * m["a21"]
    found: list[tuple[float, float, float, float]] = []
    # same acceptance and merge order as one seed at a time (seed order is preserved)
    for k in np.flatnonzero(res < tol_arcsec):
        mu = float(1.0 / det[k])
        if all(np.hypot(x[k] - u_, y[k] - v_) > 0.05 for u_, v_, _, _ in found):
            found.append((float(x[k]), float(y[k]), mu, float(res[k])))
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


def imageplane_residuals(model: LensModel, grid: DeflectionGrid, backtrace: Table) -> Table:
    """Exact image-plane residual of each catalogued image (Lenstool "image plane optimization").

    For each family with two or more traced images, the source is the mean back-traced position
    (``beta_x, beta_y`` of :func:`backtrace_images`); :func:`find_images` solves for all of its
    images, and each catalogued image is matched to the nearest predicted one (``dtheta_arcsec``;
    NaN when the family has no prediction). ``n_predicted`` counts the family's predicted images;
    ``shared_match`` marks catalogued images whose nearest predicted image is also another
    catalogued image's nearest (a predicted image is missing or merged). Provenance
    ``model_prediction``.
    """
    systems = np.asarray(backtrace["system"]).astype(str)
    dtheta = np.full(len(backtrace), np.nan)
    mu = np.full(len(backtrace), np.nan)
    n_pred = np.zeros(len(backtrace), int)
    shared = np.zeros(len(backtrace), bool)
    for sys_id in dict.fromkeys(systems):
        sel = np.where((systems == sys_id) & np.isfinite(backtrace["beta_x"]))[0]
        if len(sel) < 2:
            continue
        obs = backtrace[sel]
        pred = find_images(
            model,
            grid,
            float(np.mean(obs["beta_x"])),
            float(np.mean(obs["beta_y"])),
            float(obs["z_used"][0]),
        )
        n_pred[sel] = len(pred)
        if not len(pred):
            continue
        ox, oy = model.to_frame(np.asarray(obs["ra"], float), np.asarray(obs["dec"], float))
        d = np.hypot(ox[:, None] - pred["x"][None, :], oy[:, None] - pred["y"][None, :])
        k = np.argmin(d, axis=1)
        dtheta[sel] = d[np.arange(len(sel)), k]
        mu[sel] = np.asarray(pred["magnification"])[k]
        shared[sel] = np.bincount(k, minlength=len(pred))[k] > 1
    out = Table(
        {
            "image_id": backtrace["image_id"],
            "system": systems,
            "z_used": backtrace["z_used"],
            "dtheta_arcsec": dtheta,
            "magnification": mu,
            "n_predicted": n_pred,
            "shared_match": shared,
        }
    )
    out.meta.update(model._meta())
    out.meta["source"] = f"{backtrace.meta.get('source', 'images')}, image-plane solve"
    return out


class MapLensModel(LensModel):
    """A lens model given as published deflection maps (RELICS / HFF / UNCOVER convention).

    ``alpha_x``, ``alpha_y`` are deflection maps in arcsec at D_LS/D_S = 1 along the image's pixel
    axes, on a north-up, east-left TAN grid (``wcs``; rotated grids raise
    :class:`UnsupportedModelError`): +i runs West and +j North, as the model frame's +x and +y.
    The model frame is centred on ``centre`` (RA, Dec; default: the map's reference pixel), with
    x = West and y = North in arcsec as for :class:`LensModel`. Positions are mapped to pixels
    through the WCS. Deflection is interpolated bilinearly; the Hessian comes from centred finite
    differences in float64 (so κ, γ and μ are resolution-limited at critical curves). Positions
    outside the maps give NaN. Provenance of every prediction: ``model_prediction``.
    """

    def __init__(
        self,
        alpha_x: np.ndarray,
        alpha_y: np.ndarray,
        wcs: Any,
        z_lens: float,
        cosmology: FlatLambdaCDM,
        source: str = "",
        sha256: str = "",
        centre: tuple[float, float] | None = None,
    ) -> None:
        cd = np.asarray(wcs.pixel_scale_matrix, float) * 3600.0
        if abs(cd[0, 1]) > 1e-9 * abs(cd[0, 0]) or abs(cd[1, 0]) > 1e-9 * abs(cd[1, 1]):
            raise UnsupportedModelError(f"{source}: rotated deflection maps are not supported")
        if not (cd[0, 0] < 0 < cd[1, 1]) or not np.isclose(-cd[0, 0], cd[1, 1], rtol=1e-6):
            raise UnsupportedModelError(f"{source}: maps must be north-up, east-left, square")
        if alpha_x.shape != alpha_y.shape:
            raise ValueError(f"{source}: deflection maps differ in shape")
        self.components = ()
        self.z_lens = float(z_lens)
        self.wcs = wcs
        if centre is None:
            ref = (float(wcs.wcs.crpix[0]) - 1.0, float(wcs.wcs.crpix[1]) - 1.0)
            centre = tuple(float(v) for v in wcs.pixel_to_world_values(*ref))
        self.ra0, self.dec0 = float(centre[0]), float(centre[1])
        self.cosmology = cosmology
        self.source = source
        self.sha256 = sha256
        self._cos0 = np.cos(np.deg2rad(self.dec0))
        self.pixel_arcsec = float(cd[1, 1])
        step = self.pixel_arcsec
        ax = np.asarray(alpha_x, np.float64)
        ay = np.asarray(alpha_y, np.float64)
        dax_dy, dax_dx = np.gradient(ax, step)
        day_dy, day_dx = np.gradient(ay, step)
        self._maps = {
            "alpha_x": ax.astype(np.float32),
            "alpha_y": ay.astype(np.float32),
            "psi_xx": dax_dx.astype(np.float32),
            "psi_yy": day_dy.astype(np.float32),
            "psi_xy": (0.5 * (dax_dy + day_dx)).astype(np.float32),
        }
        del ax, ay, dax_dy, dax_dx, day_dy, day_dx
        self.shape = alpha_x.shape

    @classmethod
    def from_fits(
        cls,
        alpha_x_path: str | Path,
        alpha_y_path: str | Path,
        z_lens: float,
        cosmology: FlatLambdaCDM,
        source: str = "",
        centre: tuple[float, float] | None = None,
    ) -> MapLensModel:
        """Build from two FITS deflection maps (arcsec, D_LS/D_S = 1); ``sha256`` hashes both
        files (read in chunks)."""
        import warnings

        from astropy.io import fits
        from astropy.wcs import WCS, FITSFixedWarning

        digest = hashlib.sha256()
        arrays = []
        wcs = None
        for p in (alpha_x_path, alpha_y_path):
            with open(p, "rb") as fh:
                while chunk := fh.read(1 << 22):
                    digest.update(chunk)
            with fits.open(p, memmap=True) as hdul:
                arrays.append(np.asarray(hdul[0].data, np.float32))
                if wcs is None:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", FITSFixedWarning)
                        wcs = WCS(hdul[0].header)
        return cls(
            arrays[0],
            arrays[1],
            wcs,
            z_lens,
            cosmology,
            source=source or f"{alpha_x_path} + {alpha_y_path}",
            sha256=digest.hexdigest(),
            centre=centre,
        )

    def shift_frame(self, dra_arcsec: float, ddec_arcsec: float) -> None:
        """Move the model and its maps on the sky (the maps are looked up by sky position, so
        moving only the frame anchor would leave them in place; D-040)."""
        crval = self.wcs.wcs.crval
        self.wcs.wcs.crval = [
            crval[0] + dra_arcsec / 3600.0 / self._cos0,
            crval[1] + ddec_arcsec / 3600.0,
        ]
        self.wcs.wcs.set()
        super().shift_frame(dra_arcsec, ddec_arcsec)

    def _interp(self, x: np.ndarray, y: np.ndarray, keys) -> dict[str, np.ndarray]:
        from scipy.ndimage import map_coordinates

        ra, dec = self.to_sky(x, y)
        i, j = (np.asarray(v, float) for v in self.wcs.world_to_pixel_values(ra, dec))
        inside = (i >= 0) & (i <= self.shape[1] - 1) & (j >= 0) & (j <= self.shape[0] - 1)
        coords = np.vstack([np.ravel(j), np.ravel(i)])
        out = {}
        for k in keys:
            v = map_coordinates(self._maps[k], coords, order=1, mode="nearest").reshape(x.shape)
            out[k] = np.where(inside, v.astype(float), np.nan)
        return out

    def fields_xy(self, x: Any, y: Any) -> dict[str, np.ndarray]:
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        return self._interp(x, y, ("alpha_x", "alpha_y", "psi_xx", "psi_xy", "psi_yy"))

    def deflection_xy(self, x: Any, y: Any) -> tuple[np.ndarray, np.ndarray]:
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        f = self._interp(x, y, ("alpha_x", "alpha_y"))
        return f["alpha_x"], f["alpha_y"]

    def kappa_xy(self, x: Any, y: Any) -> np.ndarray:
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        f = self._interp(x, y, ("psi_xx", "psi_yy"))
        return 0.5 * (f["psi_xx"] + f["psi_yy"])

    def _meta(self) -> dict[str, Any]:
        meta = super()._meta()
        meta.pop("n_potentials", None)
        meta.update(
            source=f"deflection maps {self.source}",
            maps={
                "shape": list(self.shape),
                "pixel_arcsec": self.pixel_arcsec,
                "interpolation": "bilinear; Hessian from centred finite differences",
            },
        )
        return meta


def _cosmology_key(cosmo: FlatLambdaCDM) -> str:
    """The cosmological parameters (not the name or meta) as text, for content hashes."""

    def lossless(v: Any) -> str:
        return f"{np.asarray(getattr(v, 'value', v)).tolist()!r} {getattr(v, 'unit', '')}"

    return repr(sorted((k, lossless(v)) for k, v in cosmo.parameters.items()))


def _detached(plane: LensModel) -> LensModel:
    """A shallow copy whose frame can move without moving the caller's object (or its WCS)."""
    q = copy.copy(plane)
    if getattr(q, "wcs", None) is not None:
        q.wcs = q.wcs.deepcopy()
    return q


class MultiPlaneLensModel:
    """Several lens planes, each a :class:`LensModel` at its own ``z_lens``, in one frame (D-046).

    Typical use: a single-plane cluster model with a foreground or background galaxy that the
    model treats as a cluster member, moved to its spectroscopic redshift with
    :meth:`LensModel.split_planes`. The planes share the reference point and cosmology. The lens
    equation is the standard multi-plane recursion (e.g. Schneider, Ehlers & Falco 1992, ch. 9):

    ``theta_j = theta_1 - sum_{i<j} (D_ij / D_j) alpha_i(theta_i)``, the source position
    ``beta = theta_1 - sum_i (D_is / D_s) alpha_i(theta_i)``, and
    ``A_j = I - sum_{i<j} (D_ij / D_j) H_i(theta_i) A_i`` (``A = I - sum_i (D_is / D_s) H_i A_i``),

    with ``alpha_i`` and ``H_i`` each plane's deflection and Hessian at D_LS/D_S = 1 and, in a flat
    cosmology, ``D_ij / D_j = 1 - D_M(z_i) / D_M(z_j)``. The ray positions ``theta_i`` do not depend
    on the source redshift, only the weights ``D_is / D_s`` do, so a :class:`DeflectionGrid` stores
    ``alpha_i(theta_i)`` for every plane (shape ``(n_planes, n, n)``) and serves every source
    redshift. It provides what :func:`backtrace_images`, :func:`find_images` and
    :func:`imageplane_residuals` need; it has no ``evaluate`` (kappa/shear tables) and no scalar
    ``z_lens`` (see ``z_planes``), so single-plane-only code fails loudly on it.
    """

    def __init__(self, planes: Sequence[LensModel], source: str = "", sha256: str = "") -> None:
        planes = sorted(planes, key=lambda p: p.z_lens)
        if not planes:
            raise ValueError("a multi-plane model needs at least one plane")
        z = [p.z_lens for p in planes]
        if len(set(z)) != len(z):
            raise ValueError(f"lens planes must have distinct redshifts (z_lens {z})")
        p0 = planes[0]
        for p in planes[1:]:
            if (p.ra0, p.dec0) != (p0.ra0, p0.dec0) or not p.cosmology.is_equivalent(p0.cosmology):
                raise ValueError("lens planes must share the reference point and cosmology")
        self.planes = tuple(_detached(p) for p in planes)
        self.components = tuple(c for p in planes for c in p.components)
        self.ra0, self.dec0 = p0.ra0, p0.dec0
        self._cos0 = p0._cos0
        self.cosmology = p0.cosmology
        self.source = source or " + ".join(p.source for p in planes)
        # identity for caches and provenance (see _identity)
        self.sha256 = sha256 or self._identity()
        self._d_m = [
            float(self.cosmology.comoving_transverse_distance(p.z_lens).to_value(u.Mpc))
            for p in planes
        ]
        # D_ij / D_j for every pair of planes i < j
        self._ratio = {
            (i, j): 1.0 - self._d_m[i] / self._d_m[j] for j in range(len(planes)) for i in range(j)
        }

    def _identity(self) -> str:
        """Each plane's file hash, redshift and potentials, the frame anchor and the cosmology.
        Parametric planes are fully described by these. A plane with neither a hash nor
        potentials (an in-memory map) is not, so then the model has no identity: "" is never
        cached, and ``find_images`` cannot check that a grid belongs to it."""
        if any(not p.sha256 and not p.components for p in self.planes):
            return ""
        ident = repr(
            [(p.sha256, p.z_lens, p.components) for p in self.planes]
            + [self.ra0, self.dec0, _cosmology_key(self.cosmology)]
        )
        return hashlib.sha256(ident.encode()).hexdigest()

    to_frame = LensModel.to_frame
    to_sky = LensModel.to_sky
    frame_angle_to_pa = staticmethod(LensModel.frame_angle_to_pa)

    @property
    def z_planes(self) -> tuple[float, ...]:
        return tuple(p.z_lens for p in self.planes)

    @property
    def z_lens(self) -> float:
        raise AttributeError(
            "MultiPlaneLensModel has no single z_lens (see z_planes): single-plane code path"
        )

    def shift_frame(self, dra_arcsec: float, ddec_arcsec: float) -> None:
        """Move every plane on the sky by (dRA cos dec, dDec) arcsec, in place (D-034)."""
        for p in self.planes:
            p.shift_frame(dra_arcsec, ddec_arcsec)
        p0 = self.planes[0]
        self.ra0, self.dec0, self._cos0 = p0.ra0, p0.dec0, p0._cos0
        # model-frame deflections do not change (a map plane's WCS moves too), so the identity,
        # and any grid computed before the shift, stay valid, as for one plane

    def _rays(
        self, x: np.ndarray, y: np.ndarray, n_planes: int, jacobian: bool
    ) -> tuple[list[tuple[np.ndarray, np.ndarray]], list[tuple[np.ndarray, ...]]]:
        """``alpha_i(theta_i)`` (and ``H_i A_i``) for the first ``n_planes`` planes; no z_s."""
        one, zero = np.ones(x.shape), np.zeros(x.shape)
        alphas: list[tuple[np.ndarray, np.ndarray]] = []
        ha: list[tuple[np.ndarray, ...]] = []
        for j, plane in enumerate(self.planes[:n_planes]):
            px, py = x.copy(), y.copy()
            a = [one.copy(), zero.copy(), zero.copy(), one.copy()] if jacobian else []
            for i in range(j):
                r = self._ratio[(i, j)]
                px -= r * alphas[i][0]
                py -= r * alphas[i][1]
                for k in range(len(a)):
                    a[k] -= r * ha[i][k]
            if jacobian:
                f = plane.fields_xy(px, py)
                hxx, hxy, hyy = f["psi_xx"], f["psi_xy"], f["psi_yy"]
                alphas.append((f["alpha_x"], f["alpha_y"]))
                ha.append(
                    (
                        hxx * a[0] + hxy * a[2],
                        hxx * a[1] + hxy * a[3],
                        hxy * a[0] + hyy * a[2],
                        hxy * a[1] + hyy * a[3],
                    )
                )
            else:
                ax, ay = plane.deflection_xy(px, py)
                alphas.append((np.asarray(ax), np.asarray(ay)))
        return alphas, ha

    def _weights(self, z_s: Any, shape: tuple[int, ...]) -> list[np.ndarray]:
        """``D_is / D_s`` per plane (0 for a plane at or behind the source), one distance
        evaluation per distinct redshift."""
        zs = np.broadcast_to(np.asarray(z_s, float), shape)
        u_z, inv = np.unique(zs.ravel(), return_inverse=True)
        d_s = np.full(u_z.shape, np.nan)
        fin = np.isfinite(u_z)
        lensed = fin & (u_z > self.planes[0].z_lens)  # only sources behind some plane
        d_s[lensed] = self.cosmology.comoving_transverse_distance(u_z[lensed]).to_value(u.Mpc)
        out = []
        for p, d_l in zip(self.planes, self._d_m, strict=True):  # as LensModel.dls_ds
            w = np.where(fin, 0.0, np.nan)
            behind = fin & (u_z > p.z_lens)
            w[behind] = 1.0 - d_l / d_s[behind]
            out.append(w[inv].reshape(shape))
        return out

    def _trace(self, x: Any, y: Any, z_s: Any, jacobian: bool) -> dict[str, np.ndarray]:
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        w = self._weights(z_s, x.shape)
        # planes behind every source never deflect (and neither do the ones behind them)
        n = len(self.planes) - int(sum(np.all(wi == 0) for wi in w))
        alphas, ha = self._rays(x, y, n, jacobian)
        out = {"beta_x": x.copy(), "beta_y": y.copy()}
        if jacobian:
            out.update(a11=np.ones(x.shape), a12=np.zeros(x.shape))
            out.update(a21=np.zeros(x.shape), a22=np.ones(x.shape))
        for i in range(n):
            on = w[i] != 0  # a plane behind this source adds nothing, even a NaN off its map
            out["beta_x"] = out["beta_x"] - np.where(on, w[i] * alphas[i][0], 0.0)
            out["beta_y"] = out["beta_y"] - np.where(on, w[i] * alphas[i][1], 0.0)
            if jacobian:
                for k, key in enumerate(("a11", "a12", "a21", "a22")):
                    out[key] = out[key] - np.where(on, w[i] * ha[i][k], 0.0)
        return out

    def lens_map(self, x: Any, y: Any, z_s: Any) -> dict[str, np.ndarray]:
        """As :meth:`LensModel.lens_map`; ``A`` is not symmetric in general (``a12 != a21``)."""
        return self._trace(x, y, z_s, jacobian=True)

    def source_points(self, x: Any, y: Any, z_s: Any) -> tuple[np.ndarray, np.ndarray]:
        out = self._trace(x, y, z_s, jacobian=False)
        return out["beta_x"], out["beta_y"]

    def plane_deflections(self, x: Any, y: Any) -> tuple[np.ndarray, np.ndarray]:
        """``alpha_i(theta_i)`` of every plane along the rays through ``(x, y)``, stacked on a
        leading plane axis (independent of the source redshift; for :class:`DeflectionGrid`)."""
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        alphas, _ = self._rays(x, y, len(self.planes), jacobian=False)
        return np.stack([a[0] for a in alphas]), np.stack([a[1] for a in alphas])

    def source_grid(self, grid: DeflectionGrid, z_s: float) -> tuple[np.ndarray, np.ndarray]:
        """Source-plane position of every node of ``grid`` (computed for this model)."""
        if grid.alpha_x.ndim != 3 or grid.alpha_x.shape[0] != len(self.planes):
            raise ValueError("the deflection grid was not computed for this multi-plane model")
        g, shape = grid.x, grid.alpha_x.shape[1:]
        bx = np.broadcast_to(g[None, :], shape).astype(float)
        by = np.broadcast_to(g[:, None], shape).astype(float)
        for i, w in enumerate(self._weights(float(np.squeeze(z_s)), ())):
            r = float(w)
            if r:
                bx -= r * grid.alpha_x[i]
                by -= r * grid.alpha_y[i]
        return bx, by

    def _meta(self) -> dict[str, Any]:
        meta = self.planes[0]._meta()
        meta.update(
            source=f"Lenstool model {self.source}",
            model_sha256=self.sha256,
            n_potentials=len(self.components),
            z_lens=list(self.z_planes),
            lens_planes=[
                {
                    "z_lens": p.z_lens,
                    "n_potentials": len(p.components),
                    **(
                        {"potentials": [c.name for c in p.components]}
                        if len(p.components) <= 10
                        else {}
                    ),
                }
                for p in self.planes
            ],
        )
        return meta


# --- Lenstool MCMC posterior (bayes.dat) -------------------------------------------------------
# bayes.dat column "O<i> : <key> (unit)" -> potential i (1-based, best.par order) and keyword.
_BAYES_KEYS = {
    "x": "x_centre",
    "y": "y_centre",
    "emass": "ellipticity",
    "theta": "angle_pos",
    "rc": "core_radius",
    "rcut": "cut_radius",
    "sigma": "v_disp",
}
_BAYES_UNITS = {"x": "arcsec", "y": "arcsec", "emass": "", "theta": "deg", "rc": "arcsec",
                "rcut": "arcsec", "sigma": "km/s"}  # fmt: skip
_BAYES_POT = re.compile(r"^O(\d+)\s*:\s*(\w+)\s*(?:\((.*)\))?\s*$")
_BAYES_POTFILE = re.compile(r"^Pot0\s*:?\s*(rcut|sigma)\s*(?:\((.*)\))?\s*$")
_BAYES_REDSHIFT = re.compile(r"^Redshift of (\S+)")


def read_lenstool_bayes(path: str | Path) -> Table:
    """Read a Lenstool ``bayes.dat`` (MCMC samples) into a Table.

    The header has one ``#<name>`` line per column; the names are kept as written
    (``"O1 : x (arcsec)"``, ``"Pot0 sigma (km/s)"``, ``"Redshift of 7c"``, ``"Chi2"``).
    ``meta["sha256"]`` identifies the file. Provenance ``model_prediction`` (posterior samples of
    a published model).
    """
    path = Path(path)
    raw = path.read_bytes()
    text = raw.decode("latin-1")
    names = [line[1:].strip() for line in text.splitlines() if line.startswith("#")]
    data = np.loadtxt(io.StringIO(text), comments="#", ndmin=2)
    if not names or not data.size:
        raise ValueError(f"{path}: no header or no samples")
    if data.shape[1] != len(names):
        raise ValueError(f"{path}: {len(names)} header names but {data.shape[1]} columns")
    out = Table(data, names=names)
    out.meta.update(
        provenance="model_prediction", source=str(path), sha256=hashlib.sha256(raw).hexdigest()
    )
    return out


def _bayes_columns(bayes: Table) -> tuple[dict, dict, dict]:
    """``(potential, potfile, redshift)`` column maps of a bayes Table.

    potential: ``name -> (index0, keyword)``; potfile: ``name -> "sigma" | "rcut"``;
    redshift: ``name -> image family``. Unknown columns other than the bookkeeping ones raise
    :class:`UnsupportedModelError` (never silently ignored)."""
    pot, potfile, zcol = {}, {}, {}
    for name in bayes.colnames:
        if m := _BAYES_POT.match(name):
            key, unit = m.group(2), m.group(3) or ""
            if key not in _BAYES_KEYS or unit != _BAYES_UNITS[key]:
                raise UnsupportedModelError(f"bayes.dat column {name!r} is not interpreted")
            pot[name] = (int(m.group(1)) - 1, _BAYES_KEYS[key])
        elif m := _BAYES_POTFILE.match(name):
            if (m.group(2) or "") != _BAYES_UNITS[m.group(1)]:
                raise UnsupportedModelError(f"bayes.dat column {name!r} is not interpreted")
            potfile[name] = m.group(1)
        elif m := _BAYES_REDSHIFT.match(name):
            zcol[name] = image_family(m.group(1))
        elif name not in ("Nsample", "ln(Lhood)", "Chi2"):
            raise UnsupportedModelError(f"bayes.dat column {name!r} is not interpreted")
    return pot, potfile, zcol


def _check_potential_indices(par: dict[str, Any], pot: dict) -> None:
    for name, (i, _) in pot.items():
        if i >= len(par["potentials"]):
            raise UnsupportedModelError(
                f"bayes.dat column {name!r}: best.par has no potential {i + 1}"
            )


def best_sample_index(par: dict[str, Any], bayes: Table) -> tuple[int, float]:
    """The bayes row whose optimised potential parameters and family redshifts match ``par`` (a
    best.par), and that row's largest absolute difference from ``par`` over those parameters.

    A difference above best.par's print precision means best.par is not a row of this chain
    (e.g. a thinned chain); use :func:`potfile_reference` for the potfile scaling then."""
    pot, _, zcol = _bayes_columns(bayes)
    if not pot:
        raise ValueError("bayes.dat has no optimised potential columns")
    _check_potential_indices(par, pot)
    cols = [(name, par["potentials"][i][key]) for name, (i, key) in pot.items()]
    cols += [(name, par["z_m_limit"].get(fam, np.nan)) for name, fam in zcol.items()]
    diffs = np.array([np.asarray(bayes[name], float) - ref for name, ref in cols])
    diffs = diffs[np.all(np.isfinite(diffs), axis=1)]
    scale = np.maximum(np.std(diffs, axis=1, keepdims=True), 1e-9)
    k = int(np.argmin(np.sum((diffs / scale) ** 2, axis=0)))
    return k, float(np.max(np.abs(diffs[:, k])))


def potfile_reference(par: dict[str, Any], mag0: float) -> dict[str, float]:
    """Best-fit potfile normalisation (``{"sigma", "rcut"}``) from the member at ``mag0``.

    Lenstool scales every potfile member from the reference magnitude ``mag0`` (the input
    ``.par`` file's potfile ``mag0``), so the best.par member with that magnitude carries
    sigma* and rcut* themselves. Raises ``ValueError`` when no member has it."""
    refs = {
        (p["v_disp"], p["cut_radius"])
        for p in par["potentials"]
        if abs(p.get("mag", np.nan) - mag0) < 1e-4
    }
    if len(refs) != 1:
        raise ValueError(f"{par['source']}: {len(refs)} distinct potfile members with mag {mag0}")
    sigma, rcut = refs.pop()
    return {"sigma": sigma, "rcut": rcut}


def posterior_par(
    par: dict[str, Any], bayes: Table, row: int, reference: int | dict[str, float]
) -> dict[str, Any]:
    """A copy of ``par`` (from :func:`parse_lenstool_par`) with the parameters of bayes ``row``.

    * ``O<i>`` columns replace keywords of the i-th potential (best.par order). A replaced
      radius drops its ``_kpc`` twin, so the two cannot disagree.
    * ``Pot0 sigma`` / ``Pot0 rcut`` (the potfile scaling-relation normalisation) multiply
      ``v_disp`` / ``cut_radius`` of every potential with a ``mag`` (a potfile member) by their
      ratio to the best-fit values: ``reference`` is best.par's row (:func:`best_sample_index`)
      or ``{"sigma", "rcut"}`` (:func:`potfile_reference`). The member scaling relations are
      power laws in luminosity, so this ratio is exact.
    * ``Redshift of <id>`` columns replace the fixed ``z_m_limit`` redshift of every family on
      ``<id>``'s ``z_m_limit`` line (one line may list several ids, e.g. Bergamini+2023b's
      A200.1a B200.2a C200.3a).

    ``sha256`` gets a ``#<Nsample or row>`` suffix so caches keyed on it stay per sample.
    """
    pot, potfile, zcol = _bayes_columns(bayes)
    _check_potential_indices(par, pot)
    out = copy.deepcopy(par)
    r = bayes[row]
    sampled = set()
    for name, (i, key) in pot.items():
        p = out["potentials"][i]
        p[key] = float(r[name])
        sampled.add((i, key))
        if key in ("core_radius", "cut_radius"):
            p[f"{key}_kpc"] = float("nan")
    for name, kind in potfile.items():
        ref = reference[kind] if isinstance(reference, dict) else bayes[reference][name]
        ratio = float(r[name]) / float(ref)
        key = "v_disp" if kind == "sigma" else "cut_radius"
        for i, p in enumerate(out["potentials"]):
            # a member optimised on its own (an O<i> column for this key) is not rescaled
            if np.isfinite(p.get("mag", np.nan)) and (i, key) not in sampled:
                p[key] *= ratio
                if key == "cut_radius":
                    p["cut_radius_kpc"] = float("nan")
    groups = par.get("z_m_limit_groups") or [[f] for f in par["z_m_limit"]]
    for name, family in zcol.items():
        group = next((g for g in groups if family in g), None)
        if group is None:
            raise UnsupportedModelError(f"bayes.dat column {name!r}: no z_m_limit for {family}")
        for fam in group:
            out["z_m_limit"][fam] = float(r[name])
    tag = int(r["Nsample"]) if "Nsample" in bayes.colnames else row
    out["sha256"] = f"{par['sha256']}#{tag}-{row}"
    return out
