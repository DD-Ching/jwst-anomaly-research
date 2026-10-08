"""Published strong-lens catalogues as one survey, and the "no visible deflector" screen (D-056).

W1 (negative-mass lens) and W2 (Ellis wormhole) predict lensed images with no luminous deflector
(D-047). Before running any lens finder, this module reuses the published candidate and confirmed
lens lists (D-054): lenscat ``catalog.csv`` (Vujeva et al. 2024, MIT), the Euclid Q1 Strong
Lensing Discovery Engine catalogue and lens-model tables (Walmsley et al. 2025, CC-BY-4.0) and the
SuGOHI public list (HSC-SSP). Each reader returns one row per catalogue entry in
:data:`ENTRY_COLUMNS` (``observed``: values as published); :func:`merge` groups entries within a
radius into systems and keeps every entry's provenance; :class:`PublishedLensSurvey` exposes the
merged list through the ``signatures.CatalogueSurvey`` protocol.

The screen itself (:func:`classify_deflectors`) asks whether a deep-imaging catalogue (Legacy
Surveys DR10 Tractor in ``scripts/w12_lenscats.py``) has a galaxy at the lens position, and
:func:`required_lens_mag` estimates how bright an ordinary lens galaxy of the needed velocity
dispersion would be (singular isothermal sphere plus an empirical Faber-Jackson calibration on
the catalogues' own visible lenses; ASSUMPTION), so a non-detection can be compared with the
survey depth. A lens without a visible deflector is an anomaly to vet, never a discovery.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np
from astropy.cosmology import Planck18
from astropy.table import Table
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from jwst_anomaly import schema

C_KM_S = 299792.458

# Pinned inputs (SOURCES.md "Published lens catalogues (D-056)"); fetched with
# photometry.fetch_catalog, which verifies sha256 before caching.
CATALOGUES: dict[str, dict] = {
    "lenscat": {
        "url": "https://raw.githubusercontent.com/lenscat/lenscat/"
        "f531b8a8f3fa4bbd1ee8a58b2ddcd54e6936ea14/lenscat/data/catalog.csv",
        "sha256": "7de5111afb6486c119198cb2df868a2e6a15f79869e6cdc8109ec7d59b06a7cc",
        "licence": "MIT",
    },
    "euclid_q1": {
        "url": "https://zenodo.org/api/records/15025832/files/"
        "q1_discovery_engine_lens_catalog.csv/content",
        "sha256": "ee5e60cd507413eabf3da8ffa37a3c527212da00feb29083cf516a86d6f26877",
        "licence": "CC-BY-4.0",
    },
    "euclid_q1_mass": {
        "url": "https://zenodo.org/api/records/15025832/files/modeling_lens_mass.csv/content",
        "sha256": "f2a52616a1ac65137b34abe0892d251e18e2c60e9a64c459914102f44ec49db5",
        "licence": "CC-BY-4.0",
    },
    "euclid_q1_sersic_mag": {
        "url": "https://zenodo.org/api/records/15025832/files/"
        "modeling_sersic_magnitude.csv/content",
        "sha256": "7c7506af33f27d65cc618e998cbbb0868c36a106703f303d6b0b0a74240d94f2",
        "licence": "CC-BY-4.0",
    },
    "sugohi": {
        # Served by a PHP script: the list can be updated in place, so the pin may need refreshing.
        "url": "https://www-utap.phys.s.u-tokyo.ac.jp/~oguri/sugohi/"
        "download_list.php?file=list_ra_asc_public.csv",
        "sha256": "72fb96dc8d13851c20304b25b8889087d92df27cf1549e6522c8657befda4a55",
        "licence": "none stated",
    },
}

ENTRY_COLUMNS = (
    "entry_id",  # "<catalogue>:<name>"
    "catalogue",
    "name",
    "ra",  # deg, as published
    "dec",
    "scale",  # galaxy | group | cluster | unknown
    "grade",  # A | B | C | ? (normalised; native grade kept in grade_native)
    "grade_native",
    "z_lens",  # NaN when not given
    "lens_z_known",  # the catalogue says a lens redshift exists (even if no value is given)
    "z_source",
    "theta_e",  # arcsec, NaN when not given
    "lens_mag",  # lens-galaxy AB magnitude, NaN when not given
    "lens_mag_band",
    "system_type",  # e.g. GG (galaxy-galaxy), GQ (galaxy-quasar), as published
    "ref",
)

GRADE_ORDER = {"A": 0, "B": 1, "C": 2, "?": 3}
_MISSING = {"", "-", "not measured", "nan", "none"}


def _entries(cat: str, n: int, **cols) -> Table:
    out = Table()
    for c in ENTRY_COLUMNS:
        if c in cols:
            out[c] = cols[c]
        elif c in ("z_lens", "z_source", "theta_e", "lens_mag"):
            out[c] = np.full(n, np.nan)
        elif c == "lens_z_known":
            out[c] = np.zeros(n, bool)
        else:
            out[c] = np.full(n, "", dtype=object)
    out["catalogue"] = np.full(n, cat, dtype=object)
    out["entry_id"] = np.array([f"{cat}:{nm}" for nm in out["name"]], dtype=object)
    out.meta.update(provenance=schema.Provenance.OBSERVED.value, source=cat)
    return out


def _first_float(token: str) -> float:
    """First numeric value of a lenscat ``zlens`` field ("0.49 | 0.49" -> 0.49)."""
    for part in str(token).split("|"):
        try:
            v = float(part)
        except ValueError:
            continue
        if np.isfinite(v) and v > 0:  # 0 is a placeholder in a few rows
            return v
    return np.nan


# lenscat rows inherited from cluster surveys are typed "galaxy" but list a cluster position
# (titles checked on the arXiv API, 2026-10-08). Any of these references -> scale "cluster".
CLUSTER_SURVEY_REFS = (
    "doi:10.1086/423038",  # Lopes et al. 2004, NoSOCS
    "doi:10.1086/191426",  # Gioia et al. 1990, EMSS
    "astro-ph/0106055",  # Gonzalez et al. 2001, Las Campanas Distant Cluster Survey
    "astro-ph/0405546",  # Boehringer et al. 2004, REFLEX
    "astro-ph/0003219",  # Boehringer et al. 2000, NORAS
    "0411075",  # Gladders & Yee 2005, Red-Sequence Cluster Survey
    "astro-ph/0609815",  # Olsen et al. 2007, CFHTLS clusters
    "astro-ph/9812394",  # Ebeling et al. 1998, BCS
    "astro-ph/0003191",  # Ebeling et al. 2000, eBCS
    "1002.2226",  # Menanteau et al. 2010, Southern Cosmology Survey clusters
    "astro-ph/0403354",  # Popesso et al. 2004, RASS-SDSS clusters
    "astro-ph/0310009",  # Cypriano et al. 2004, weak lensing of Abell clusters
    "1210.4136",  # Furlanetto et al. 2013, SOAR Gravitational Arc Survey (cluster targets)
    "astro-ph/0608624",  # van Breukelen et al. 2006, UKIDSS UDS clusters at 0.6 < z < 1.4
)


def read_lenscat(path) -> Table:
    """lenscat ``catalog.csv``: name, RA, Dec, zlens, type, grading, ref (32,838 rows in 1.1.3).

    ``zlens`` is a value, several values joined by "|", "measured" (a value exists elsewhere),
    a flag such as "L", or "not measured" / "-". Anything but the last two counts as
    ``lens_z_known``. Grades: confident -> A, probable -> B (ASSUMPTION: lenscat merges graded
    lists into these two bins). Rows citing a cluster survey (:data:`CLUSTER_SURVEY_REFS`) get
    scale "cluster" whatever their type.
    """
    t = Table.read(path, format="ascii.csv", fill_values=[("", "")])
    z_tok = [str(v).strip() for v in t["zlens"]]
    grade_native = np.array([str(g) for g in t["grading"]], dtype=object)
    return _entries(
        "lenscat",
        len(t),
        name=np.array([str(v) for v in t["name"]], dtype=object),
        ra=np.asarray(t["RA [deg]"], float),
        dec=np.asarray(t["DEC [deg]"], float),
        scale=np.array(
            [
                "cluster" if any(k in str(r) for k in CLUSTER_SURVEY_REFS) else str(v).lower()
                for v, r in zip(t["type"], t["ref"], strict=True)
            ],
            dtype=object,
        ),
        grade=np.array(
            [{"confident": "A", "probable": "B"}.get(g, "?") for g in grade_native], dtype=object
        ),
        grade_native=grade_native,
        z_lens=np.array([_first_float(z) for z in z_tok]),
        lens_z_known=np.array([z.lower() not in _MISSING for z in z_tok]),
        ref=np.array([str(v) for v in t["ref"]], dtype=object),
    )


SUGOHI_COLUMNS = (
    "name ra dec zl_spec zs_spec zl_phot zs_phot rein lens_i src_i type method grade ref"
).split()


def read_sugohi(path) -> Table:
    """SuGOHI public list (headerless CSV; -99 = missing; column order of the query form).

    Type GG / GQ = galaxy-scale lens of a galaxy / quasar, CG = group or cluster scale.
    ``z_lens`` prefers the spectroscopic value; ``lens_mag`` is HSC i.
    """
    t = Table.read(path, format="ascii.no_header", delimiter=",", names=SUGOHI_COLUMNS)

    def col(k):
        v = np.asarray(t[k].filled(np.nan) if hasattr(t[k], "filled") else t[k], float)
        return np.where(v < -90, np.nan, v)

    zl = np.where(np.isfinite(col("zl_spec")), col("zl_spec"), col("zl_phot"))
    zl = np.where(zl > 0, zl, np.nan)
    zs = np.where(np.isfinite(col("zs_spec")), col("zs_spec"), col("zs_phot"))
    typ = np.array([str(v) for v in t["type"]], dtype=object)
    g = [str(v) if str(v) in ("A", "B", "C") else "?" for v in t["grade"]]
    return _entries(
        "sugohi",
        len(t),
        name=np.array([str(v) for v in t["name"]], dtype=object),
        ra=np.asarray(t["ra"], float),
        dec=np.asarray(t["dec"], float),
        scale=np.array(["group" if v.startswith("C") else "galaxy" for v in typ], dtype=object),
        grade=np.array(g, dtype=object),
        grade_native=np.array([str(v) for v in t["grade"]], dtype=object),
        z_lens=zl,
        lens_z_known=np.isfinite(zl),
        z_source=zs,
        theta_e=np.where(col("rein") > 0, col("rein"), np.nan),
        lens_mag=col("lens_i"),
        lens_mag_band=np.full(len(t), "hsc_i", dtype=object),
        system_type=typ,
        ref=np.array([str(v) for v in t["ref"]], dtype=object),
    )


def read_euclid_q1(path, mass_path=None, mag_path=None) -> Table:
    """Euclid Q1 Discovery Engine catalogue, optionally joined (on ``id_str``) with the PyAutoLens
    SIE Einstein radius (``einstein_radius_effective_median_pdf``, arcsec) and the Sersic lens
    VIS magnitude. Positions are of the host (lens) galaxy and may be off by up to 10''
    (Zenodo README). Every candidate was selected around a catalogued galaxy.
    """
    t = Table.read(path, format="ascii.csv")
    ids = np.array([str(v) for v in t["id_str"]], dtype=object)
    n = len(t)
    theta = np.full(n, np.nan)
    mag = np.full(n, np.nan)
    row = {k: i for i, k in enumerate(ids)}
    for p, col, out in (
        (mass_path, "einstein_radius_effective_median_pdf", theta),
        (mag_path, "vis_lens_magnitude_ab_median_pdf", mag),
    ):
        if p is None:
            continue
        m = Table.read(p, format="ascii.csv")
        vals = np.ma.filled(np.ma.asarray(m[col], float), np.nan)
        for k, v in zip(m["id_str"], vals, strict=True):
            if str(k) in row:
                out[row[str(k)]] = float(v)
    names = np.array(
        [i if i not in ("", "--") else f"row{j}" for j, i in enumerate(ids)], dtype=object
    )
    notes = np.array([str(v) if str(v) != "--" else "" for v in t["notes"]], dtype=object)
    return _entries(
        "euclid_q1",
        n,
        name=names,
        ra=np.asarray(t["right_ascension"], float),
        dec=np.asarray(t["declination"], float),
        scale=np.array(["group" if "cluster" in s else "galaxy" for s in notes], dtype=object),
        grade=np.array([str(g) if str(g) in "ABC" else "?" for g in t["grade"]], dtype=object),
        grade_native=np.array([str(g) for g in t["grade"]], dtype=object),
        theta_e=theta,
        lens_mag=mag,
        lens_mag_band=np.full(n, "euclid_vis", dtype=object),
        system_type=np.array([f"{s}" for s in t["subset"]], dtype=object),
        ref=notes,
    )


def _unit(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def _join_refs(refs, max_len: int = 400) -> str:
    """Unique references of a system's entries (lenscat repeats them), " | "-joined."""
    seen: list[str] = []
    for r in refs:
        for part in str(r).split(" | "):
            part = part.strip()
            if part and part not in seen:
                seen.append(part)
    return " | ".join(seen)[:max_len]


def merge(tables: list[Table], radius_arcsec: float = 3.0) -> Table:
    """Group entries of all catalogues into systems (friends-of-friends within ``radius_arcsec``).

    One row per system: the position, scale and name of its best-graded entry (catalogue order
    breaks ties), the first finite lens / source redshift, the median Einstein radius,
    ``lens_z_known`` and ``lens_mag_known`` if any entry has them, and ``entries``
    (";"-joined ``entry_id``) as the per-entry provenance, and ``refs`` (unique references).
    Entries of one catalogue closer than the radius merge too (SuGOHI already merges within
    3''). ``derived``.
    """
    from astropy.table import vstack

    allt = vstack([t[list(ENTRY_COLUMNS)] for t in tables], metadata_conflicts="silent")
    n = len(allt)
    xyz = _unit(np.asarray(allt["ra"], float), np.asarray(allt["dec"], float))
    chord = 2 * np.sin(np.radians(radius_arcsec / 3600) / 2)
    pairs = cKDTree(xyz).query_pairs(chord, output_type="ndarray")
    graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    _, label = connected_components(graph, directed=False)
    order = sorted(range(n), key=lambda i: (label[i], GRADE_ORDER.get(allt["grade"][i], 3), i))
    rows = []
    start = 0
    while start < n:
        stop = start
        while stop < n and label[order[stop]] == label[order[start]]:
            stop += 1
        idx = order[start:stop]
        best = idx[0]

        def first(col, idx=idx):
            v = np.asarray(allt[col][idx], float)
            v = v[np.isfinite(v)]
            return float(v[0]) if len(v) else np.nan

        th = np.asarray(allt["theta_e"][idx], float)
        th = th[np.isfinite(th)]
        rows.append(
            {
                "system_id": f"L{len(rows):05d}",
                "name": allt["name"][best],
                "ra": float(allt["ra"][best]),
                "dec": float(allt["dec"][best]),
                "scale": allt["scale"][best],
                "grade": allt["grade"][best],
                "n_entries": len(idx),
                "catalogues": ";".join(sorted({allt["catalogue"][i] for i in idx})),
                "entries": ";".join(allt["entry_id"][i] for i in idx),
                "any_group_scale": bool(any(allt["scale"][i] in ("group", "cluster") for i in idx)),
                "z_lens": first("z_lens"),
                "lens_z_known": bool(np.any(np.asarray(allt["lens_z_known"][idx], bool))),
                "z_source": first("z_source"),
                "theta_e": float(np.median(th)) if len(th) else np.nan,
                "lens_mag": first("lens_mag"),
                "lens_mag_known": bool(
                    np.any(np.isfinite(np.asarray(allt["lens_mag"][idx], float)))
                ),
                "system_type": ";".join(sorted({str(allt["system_type"][i]) for i in idx} - {""})),
                "refs": _join_refs(allt["ref"][i] for i in idx),
            }
        )
        start = stop
    out = Table(rows=rows)
    out.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        source="merge of " + ", ".join(sorted(set(allt["catalogue"]))),
        merge_radius_arcsec=radius_arcsec,
        n_entries=n,
    )
    return out


@dataclass
class PublishedLensSurvey:
    """Merged published lens catalogues as a ``signatures.CatalogueSurvey``.

    ``area`` is the footprint over which a screen was run (e.g. the Legacy Surveys DR10 area the
    deflector test covers); published lists have no single area, so it is supplied by the caller
    (NaN until a screen sets it).
    """

    tables: list[Table]
    name: str = "published-lens-catalogues"
    radius_arcsec: float = 3.0
    area: float = float("nan")

    def catalogue(self) -> Table:
        return merge(self.tables, self.radius_arcsec)

    def area_deg2(self) -> float:
        return self.area


# ----------------------------------------------------------------------------- selection


def select_no_lens_info(systems: Table, theta_e_max: float = 3.0) -> np.ndarray:
    """Galaxy-scale systems that no catalogue gives a lens redshift or lens magnitude for.

    Galaxy scale: no entry is group / cluster scale, and the Einstein radius (when given) is at
    most ``theta_e_max`` arcsec (ASSUMPTION: galaxy-scale lenses have theta_E < 3''). Missing
    information is not evidence of a missing lens: this only orders the deep-imaging test.
    """
    galaxy = ~np.asarray(systems["any_group_scale"], bool)
    th = np.asarray(systems["theta_e"], float)
    galaxy &= ~(th > theta_e_max)
    no_info = ~np.asarray(systems["lens_z_known"], bool) & ~np.asarray(
        systems["lens_mag_known"], bool
    )
    return galaxy & no_info


_JNAME = re.compile(r"(\d{2})(\d{2})(\d{2}(?:\.\d+)?)([+-])(\d{2})(\d{2})(\d{2}(?:\.\d+)?)?")
_JNAME_SHORT = re.compile(r"(\d{2})(\d{2})([+-])(\d{2})(\d{2})?$")
_JNAME_TENTHS = re.compile(r"(\d{2})(\d{2})([+-])(\d{2})(\d)$")  # HHMM+DDd (0302+006)


# Survey prefixes whose designations are J2000 even without a "J" (SPT0418-47, PS10042-1713).
J2000_PREFIXES = {"SPT", "SPT-S", "PS1", "HATLAS", "H-ATLAS", "HELMS", "HERS", "ACT", "AGEL"}


def name_position_offset(name: str, ra: float, dec: float) -> float:
    """Offset (arcsec) between a J2000-style designation and the catalogued position.

    Truncated designations (``JHHMMSS+DDMMSS``, ``JHHMM+DDMM`` or ``HHMM+DDd``) are compared
    with the position truncated the same way, so a consistent entry gives a small offset (a few
    arcsec to an arcmin for short names); a large one flags a catalogue position error. NaN when
    the name has no designation.
    """
    token = str(name).split()[-1]
    prefix = re.match(r"[A-Za-z][A-Za-z0-9-]*?(?=\d{4}[+-])", token)
    if "J" in token:
        token = token[token.rindex("J") + 1 :]
    elif prefix and prefix.group(0).rstrip("-").upper() in J2000_PREFIXES:
        token = token[prefix.end() :]
    elif not token[:1].isdigit():
        return np.nan  # B1950 or survey-specific names (HE0435-1223) are not J2000
    m = _JNAME.match(token)
    if m:
        h, mi, s, sg, d, dm, ds = m.groups()
        ra_n = 15 * (int(h) + int(mi) / 60 + float(s) / 3600)
        dec_n = int(d) + int(dm) / 60 + (float(ds) if ds else 0.0) / 3600
        res_ra, res_dec = 15 * 1.0, 60.0 if ds is None else 1.0  # truncation unit (arcsec)
    elif m := _JNAME_TENTHS.match(token):
        h, mi, sg, d, dt = m.groups()
        ra_n = 15 * (int(h) + int(mi) / 60)
        dec_n = int(d) + int(dt) / 10
        res_ra, res_dec = 15 * 60.0, 360.0
    else:
        m = _JNAME_SHORT.match(token)
        if not m:
            return np.nan
        h, mi, sg, d, dm = m.groups()
        ra_n = 15 * (int(h) + int(mi) / 60)
        dec_n = int(d) + (int(dm) if dm else 0) / 60
        res_ra, res_dec = 15 * 60.0, 60.0 if dm else 3600.0
    dec_n = -dec_n if sg == "-" else dec_n
    # distance beyond the truncation cell [name, name + unit)
    dra = ((ra - ra_n + 180) % 360 - 180) * 3600 * np.cos(np.radians(dec))
    ddec = (abs(dec) - abs(dec_n)) * 3600
    ex_ra = max(0.0, -dra, dra - res_ra * np.cos(np.radians(dec)))
    ex_dec = max(0.0, -ddec, ddec - res_dec)
    return float(np.hypot(ex_ra, ex_dec))


def position_quantum_arcsec(ra: float, dec: float, printed_decimals: int = 5) -> float:
    """Coarsest rounding step consistent with a catalogued position (arcsec on the sky).

    Positions converted from truncated sexagesimal (whole seconds of RA, whole arcsec) or
    rounded decimal degrees (0.01 deg) cannot locate a deflector to ~1''. The RA step is
    multiplied by cos(dec); the larger of the RA and Dec steps is returned (0 when finer than
    0.01''). A value counts as a multiple of a step when it is one to within the rounding of
    ``printed_decimals`` decimal degrees (lenscat prints 5). Steps finer than 3.6'' are not
    tested.
    """
    tol = 0.5 * 10.0 ** (-printed_decimals) * 1.01
    cosd = np.cos(np.radians(dec))

    def step(value: float, units: list[tuple[float, float]]) -> float:
        for per_unit, arcsec in units:  # coarsest first
            x = value * per_unit
            if abs(x - round(x)) <= tol * per_unit + 1e-9:
                return arcsec
        return 0.0

    ra_units = [(1e1, 360.0 * cosd), (1e2, 36.0 * cosd), (240.0, 15.0 * cosd)]
    ra_units += [(1e3, 3.6 * cosd)]
    dec_units = [(1e1, 360.0), (1e2, 36.0), (60.0, 60.0), (1e3, 3.6)]
    return float(max(step(ra, ra_units), step(abs(dec), dec_units)))


EXTENDED_TYPES = ("REX", "DEV", "EXP", "SER")


def classify_deflectors(
    systems: Table, sources: Table, radius_arcsec: float = 1.5, cover_radius_arcsec: float = 5.0
) -> Table:
    """Per system: is there a catalogued galaxy at the deflector position? ``derived``.

    ``sources``: deep-imaging catalogue rows with ``ra, dec, type`` (Tractor morphological types;
    extended = :data:`EXTENDED_TYPES`), ``mag_z`` and ``galdepth_z`` (inverse variance of a
    galaxy-model flux in nanomaggies). Returns ``system_id``, ``covered`` (any source within
    ``cover_radius_arcsec``), ``n_ext`` and ``n_psf`` within ``radius_arcsec``, ``sep_ext`` and
    ``mag_ext`` of the nearest extended source, ``depth_z`` (5 sigma galaxy depth, median of the
    sources within the cover radius) and ``visible_deflector`` (``n_ext > 0``). The lens position is
    the catalogue position (galaxy-scale finders list the lens galaxy; ASSUMPTION).
    """
    xyz_s = _unit(np.asarray(sources["ra"], float), np.asarray(sources["dec"], float))
    xyz = _unit(np.asarray(systems["ra"], float), np.asarray(systems["dec"], float))
    tree = cKDTree(xyz_s) if len(sources) else None
    typ = np.array([str(t).strip() for t in sources["type"]])
    ext = np.isin(typ, EXTENDED_TYPES)
    magz = np.asarray(sources["mag_z"], float) if len(sources) else np.zeros(0)
    gd = np.asarray(sources["galdepth_z"], float) if len(sources) else np.zeros(0)
    rows = []
    to_arcsec = 3600 * 180 / np.pi
    for sid, p in zip(systems["system_id"], xyz, strict=True):
        near = tree.query_ball_point(p, cover_radius_arcsec / to_arcsec) if tree else []
        near = np.asarray(near, int)
        sep = np.linalg.norm(xyz_s[near] - p, axis=1) * to_arcsec if len(near) else np.zeros(0)
        inner = near[sep <= radius_arcsec]
        e = inner[ext[inner]]
        depth = gd[near][np.isfinite(gd[near]) & (gd[near] > 0)]
        if len(e):
            se = np.linalg.norm(xyz_s[e] - p, axis=1) * to_arcsec
            k = int(np.argmin(se))
            sep_e, mag_e = float(se[k]), float(magz[e[k]])
        else:
            sep_e = mag_e = np.nan
        rows.append(
            {
                "system_id": sid,
                "covered": bool(len(near)),
                "n_ext": len(e),
                "n_psf": int(np.sum(typ[inner] == "PSF")),
                "sep_ext": sep_e,
                "mag_ext": mag_e,
                "depth_z": float(galdepth_to_mag(np.median(depth))) if len(depth) else np.nan,
                "visible_deflector": bool(len(e)),
            }
        )
    out = Table(rows=rows) if rows else Table(names=["system_id"], dtype=[object])
    out.meta.update(
        provenance=schema.Provenance.DERIVED.value,
        lens_radius_arcsec=radius_arcsec,
        cover_radius_arcsec=cover_radius_arcsec,
    )
    return out


def galdepth_to_mag(galdepth, nsigma: float = 5.0):
    """Legacy Surveys ``galdepth`` (flux inverse variance, nanomaggies^-2) -> AB n-sigma depth."""
    gd = np.asarray(galdepth, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return 22.5 - 2.5 * np.log10(nsigma / np.sqrt(gd))


# ----------------------------------------------------------------------------- required lens light


def sis_sigma(theta_e_arcsec, z_l, z_s, cosmo=Planck18):
    """Velocity dispersion (km/s) of a singular isothermal sphere with Einstein radius theta_E:
    theta_E = 4 pi (sigma/c)^2 D_ls / D_s. ``model_prediction`` for a given (z_l, z_s); NaN where
    z_l >= z_s."""
    th = np.radians(np.asarray(theta_e_arcsec, float) / 3600)
    z_l = np.asarray(z_l, float)
    z_s = np.asarray(z_s, float)
    ok = z_s > z_l
    zl = np.where(ok, z_l, 0.1)
    zs = np.where(ok, z_s, 1.0)
    ratio = (cosmo.angular_diameter_distance(zl, zs) / cosmo.angular_diameter_distance(zs)).value
    return np.where(ok, C_KM_S * np.sqrt(th / (4 * np.pi * ratio)), np.nan)


@dataclass(frozen=True)
class FJCalibration:
    """Apparent lens magnitude m = a + slope * log10(sigma / 200 km/s) + 5 log10(D_L(z)/D_L(0.5))
    + k (z - 0.5), fitted to visible lenses (``derived``); ``rms`` is the scatter. ``slope`` -10
    is Faber-Jackson L proportional to sigma^4 (ASSUMPTION when fixed)."""

    a: float
    k: float
    slope: float
    rms: float
    n: int
    band: str

    def mag(self, sigma, z_l, cosmo=Planck18):
        z_l = np.asarray(z_l, float)
        dl = cosmo.luminosity_distance(z_l).value / cosmo.luminosity_distance(0.5).value
        return (
            self.a
            + self.slope * np.log10(np.asarray(sigma, float) / 200.0)
            + 5 * np.log10(dl)
            + self.k * (z_l - 0.5)
        )


def fit_fj(
    theta_e, z_l, z_s, mag, band: str, slope: float = -10.0, cosmo=Planck18
) -> FJCalibration:
    """Least-squares fit of :class:`FJCalibration` ``a`` and ``k`` (slope fixed), 3-sigma
    clipped."""
    sig = sis_sigma(theta_e, z_l, z_s, cosmo)
    z_l = np.asarray(z_l, float)
    mag = np.asarray(mag, float)
    dl = cosmo.luminosity_distance(z_l).value / cosmo.luminosity_distance(0.5).value
    y = mag - slope * np.log10(sig / 200.0) - 5 * np.log10(dl)
    good = np.isfinite(y) & np.isfinite(z_l)
    a = k = rms = np.nan
    for _ in range(5):
        A = np.column_stack([np.ones(good.sum()), z_l[good] - 0.5])
        (a, k), *_ = np.linalg.lstsq(A, y[good], rcond=None)
        res = y - a - k * (z_l - 0.5)
        rms = float(np.std(res[good]))
        new = np.isfinite(y) & np.isfinite(z_l) & (np.abs(res) < 3 * rms)
        if new.sum() == good.sum():
            break
        good = new
    return FJCalibration(float(a), float(k), slope, rms, int(good.sum()), band)


def required_lens_mag(
    theta_e, z_s, calib: FJCalibration, z_l_grid=None, n_sigma_faint: float = 2.0, cosmo=Planck18
):
    """Faintest apparent magnitude an ordinary lens galaxy producing ``theta_e`` could have.

    Maximised over lens redshift ``z_l_grid`` (default 0.1-1.5, below z_s - 0.1) and made
    ``n_sigma_faint`` x rms fainter than the calibration (ASSUMPTION: a 2-sigma under-luminous
    lens). Returns (mag, z_l at the maximum). ``model_prediction``.
    """
    grid = np.arange(0.1, 1.51, 0.05) if z_l_grid is None else np.asarray(z_l_grid, float)
    th = np.atleast_1d(np.asarray(theta_e, float))
    zs = np.broadcast_to(np.asarray(z_s, float), th.shape)
    best = np.full(th.shape, np.nan)
    zbest = np.full(th.shape, np.nan)
    for zl in grid:
        ok = zl < zs - 0.1
        if not ok.any():
            continue
        m = calib.mag(sis_sigma(th, zl, zs, cosmo), zl, cosmo) + n_sigma_faint * calib.rms
        upd = ok & np.isfinite(m) & ~(m <= best)
        best = np.where(upd, m, best)
        zbest = np.where(upd, zl, zbest)
    return best, zbest
