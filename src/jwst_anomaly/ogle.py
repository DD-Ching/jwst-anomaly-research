"""OGLE-IV published microlensing samples as a :class:`signatures.LightCurveSurvey` (D-054, D-057).

Two homogeneous samples with published photometry, fits and detection efficiencies:

- ``bulge2019``: Mróz et al. 2019, ApJS 244, 29 (arXiv:1906.02210), 5,790 events in 112 low-cadence
  Galactic-bulge fields, 2010-06-29 … 2017-11-01 (t0 window 2455377–2458118, ΔT = 2741 d);
- ``disk2020``: Mróz et al. 2020, ApJS 249, 16 (arXiv:2004.07289), 460 events in the OGLE GVS
  Galactic-plane fields, ΔT = 2650 d (the 170 "possible" events of Table B2 are not used: they
  failed the published selection).

Files are pinned by sha256 (``FILES``) and fetched with ``photometry.fetch_catalog``. Light curves
are read from the published tarballs without unpacking them to disk. Everything here is ``observed``
as published (the fit parameters in ``events()`` are the authors' PSPL fits).
"""

from __future__ import annotations

import io
import re
import tarfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.table import Table

from jwst_anomaly import photometry, schema
from jwst_anomaly.signatures import standard_light_curve

OGLE_FLUX_ZP = 18.0  # Mróz et al.: flux 1 = I of 18 mag, so F_min = 0.1 is I = 20.5
F_MIN = 0.1  # maximum negative blend flux allowed by the published selection

_BULGE = "https://www.astrouw.edu.pl/ogle/ogle4/microlensing_maps/"
_DISK = "https://www.astrouw.edu.pl/ogle/ogle4/galactic_disk_microlensing/"

# url -> (sha256, size in bytes), retrieved 2026-10-08 (SOURCES.md "OGLE-IV microlensing samples").
FILES: dict[str, tuple[str, int]] = {
    _BULGE + "README": ("9e3e62038163881521d5a895b27ba0980edd502c20dcdb65f772ff4c0d9b2136", 2688),
    _BULGE + "table3.dat": (
        "ca47555840c808967e9a257dba9071acddd499dc078efb5f5c361b3ec8fef67f",
        1542323,
    ),
    _BULGE + "table6.dat": (
        "a3e27a63b597e49563a061480c0ea37ae475e569900549d04dd8a8afa3ad50da",
        7355,
    ),
    _BULGE + "table7.dat": (
        "56455aecadfea463c9ba623d7fe459a552a1017a2532601d3c610057fa24aa2d",
        10537,
    ),
    _BULGE + "eff.tar.gz": (
        "dd5ffa37e4860dfb137691f90e93e11278c16bef705c8950d5f4a5bacc133ac1",
        15845,
    ),
    _BULGE + "phot.tar.gz": (
        "5dafa6835b8456b00eb379d1803a6f6dae46bb5eac54fb8c2f809fe4685bfe6f",
        50731904,
    ),
    _DISK + "README": ("a456c9f2f0041b81aa1b3dbc21e0be0b7957e29dae7051cc1d88534270b85f33", 2571),
    _DISK + "table_A1.txt": (
        "60990a6d207171333c8771b9dcc1240217da80ad7f272d20c5faf3432271a954",
        139773,
    ),
    _DISK + "table_B1.txt": (
        "065fea96c0f0343f268c675b3132dd36c2c0243ba6c6a1690b784549fbeec410",
        131669,
    ),
    _DISK + "eff21.tar.gz": (
        "aacd45898288269573cb73ecd936c96fb0992563f65f6752b63701129d11f563",
        213175,
    ),
    _DISK + "data.tar.gz": (
        "74868d53863167ebca08cbbbc2df6214433bfc7382f8ee82705529ca206c93c7",
        497649,
    ),
}


@dataclass(frozen=True)
class SampleSpec:
    key: str
    base: str
    events: str
    phot: str
    eff: str
    fields: str  # number of monitored sources per field
    t_min: float  # t0 window of the published efficiencies (HJD)
    t_max: float
    delta_t_days: float  # ΔT the published efficiencies are normalised to
    reference: str


SAMPLES: dict[str, SampleSpec] = {
    "bulge2019": SampleSpec(
        "bulge2019",
        _BULGE,
        "table3.dat",
        "phot.tar.gz",
        "eff.tar.gz",
        "table7.dat",
        2455377.0,
        2458118.0,
        2741.0,
        "Mróz et al. 2019, ApJS 244, 29 (arXiv:1906.02210)",
    ),
    "disk2020": SampleSpec(
        "disk2020",
        _DISK,
        "table_B1.txt",
        "data.tar.gz",
        "eff21.tar.gz",
        "table_A1.txt",
        2456200.0,  # 2012-09-29, start of ΔT (Mróz et al. 2020, Sect. 5)
        2458850.0,  # 2020-01-01
        2650.0,
        "Mróz et al. 2020, ApJS 249, 16 (arXiv:2004.07289)",
    ),
}

# Mróz et al. 2019 used Mróz et al. 2017's events in these nine high-cadence fields; their light
# curves are not in phot.tar.gz, so they are outside this sample (efficiency files exist for them).
HIGH_CADENCE_FIELDS = frozenset(
    ["BLG500", "BLG501", "BLG504", "BLG505", "BLG506", "BLG511", "BLG512", "BLG534", "BLG611"]
)

_BYTES = re.compile(r"^#\s+(\d+)-\s*(\d+)\s+(\S+)\s+(\S+)\s")


def parse_byte_table(text: str) -> Table:
    """Parse a table whose ``#`` header lists ``Bytes Format Label`` per column.

    Values are split on whitespace (labels and formats from the header): the published rows
    overflow some byte ranges (e.g. t_E > 999 d), but every row has one token per column
    once sexagesimal colons are split.
    """
    cols = []
    for line in text.splitlines():
        if not line.startswith("#"):
            continue
        m = _BYTES.match(line)
        if m:
            cols.append((m.group(3), m.group(4)))
    if not cols:
        raise ValueError("no byte-by-byte header found")
    # sexagesimal "17:13:08.00" spans three header columns (h, m, s)
    rows = [
        ln.replace(":", " ").split()
        for ln in text.splitlines()
        if ln.strip() and not ln.startswith("#")
    ]
    bad = [i for i, r in enumerate(rows) if len(r) != len(cols)]
    if bad:
        raise ValueError(f"rows {bad[:5]} do not have {len(cols)} tokens")
    out = {}
    for j, (fmt, label) in enumerate(cols):
        raw = [r[j] for r in rows]
        if fmt[0] in "FE":
            out[label] = np.array([float(v) if v not in ("", "--") else np.nan for v in raw])
        elif fmt[0] == "I":
            out[label] = np.array([int(v) if v not in ("", "--") else -1 for v in raw])
        else:
            out[label] = np.array(raw, dtype=str)
    return Table(out)


def parse_whitespace_table(text: str, labels: list[str]) -> Table:
    rows = [ln.split() for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    out = {}
    for i, label in enumerate(labels):
        vals = [r[i] for r in rows]
        try:
            out[label] = np.array(vals, dtype=float)
        except ValueError:
            out[label] = np.array(vals, dtype=str)
    return Table(out)


def parse_efficiency(text: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(log_tE_min, log_tE_max, efficiency)`` from one published efficiency file."""
    arr = np.loadtxt(io.StringIO(text), comments="#", ndmin=2)
    return arr[:, 0], arr[:, 1], arr[:, 2]


def parse_photometry(raw: bytes | str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(HJD, I, sigma_I)`` from one published photometry file."""
    if isinstance(raw, bytes):
        raw = raw.decode()
    arr = np.loadtxt(io.StringIO(raw), ndmin=2)
    if arr.size == 0:
        return np.empty(0), np.empty(0), np.empty(0)
    return arr[:, 0], arr[:, 1], arr[:, 2]


def eff_at(log_lo: np.ndarray, log_hi: np.ndarray, eff: np.ndarray, t_e_days) -> np.ndarray:
    """Published efficiency of the log t_E bin holding ``t_e_days``; 0 outside the grid."""
    lt = np.log10(np.atleast_1d(np.asarray(t_e_days, float)))
    i = np.searchsorted(log_hi, lt, side="left")
    ok = (lt >= log_lo[0]) & (lt <= log_hi[-1]) & (i < len(eff))
    return np.where(ok, eff[np.clip(i, 0, len(eff) - 1)], 0.0)


class OgleMrozSample:
    """Adapter for one Mróz et al. OGLE-IV sample (``bulge2019`` or ``disk2020``)."""

    def __init__(self, sample: str = "bulge2019", cache_dir: Path | None = None):
        if sample not in SAMPLES:
            raise KeyError(f"unknown sample {sample!r}; one of {sorted(SAMPLES)}")
        self.spec = SAMPLES[sample]
        self.name = f"ogle4-{sample}"
        self._cache_dir = cache_dir
        self._phot: dict[str, bytes] | None = None
        self._eff: dict[str, tuple] | None = None

    # ------------------------------------------------------------------ files
    def path(self, name: str) -> Path:
        url = self.spec.base + name
        sha, _size = FILES[url]
        return photometry.fetch_catalog(url, sha, cache_dir=self._cache_dir)

    def _text(self, name: str) -> str:
        return self.path(name).read_text()

    # ------------------------------------------------------------------ protocol
    def events(self) -> Table:
        t = parse_byte_table(self._text(self.spec.events))
        out = Table(
            {
                "event_id": t["name"],
                "field": np.array([f.split(".")[0] for f in t["field"]]),
                "subfield": t["field"],
                "ra": t["ra"],
                "dec": t["dec"],
                "glon": t["glon"],
                "glat": t["glat"],
                "t0_pub": t["t0_best"],
                "tE_pub": t["tE_best"],
                "u0_pub": t["u0_best"],
                "Is_pub": t["Is_best"],
                "fs_pub": t["fs_best"],
                "tE_err1_pub": t["tE_err1"],
                "tE_err2_pub": t["tE_err2"],
                "weight_pub": t["weight"],
                "alt_id": t["ews_id"] if "ews_id" in t.colnames else t["other_id"],
            }
        )
        out.meta.update(
            provenance=schema.Provenance.OBSERVED.value,
            source=f"{self.spec.base}{self.spec.events} ({self.spec.reference}); "
            "*_pub columns are the published PSPL fits",
        )
        return out

    def _load_phot(self) -> dict[str, bytes]:
        if self._phot is None:
            phot = {}
            with tarfile.open(self.path(self.spec.phot), "r:gz") as tar:
                for member in tar:  # one streaming pass; nothing is written to disk
                    if member.isfile() and member.name.endswith(".dat"):
                        key = member.name.rsplit("/", 1)[-1][: -len(".dat")]
                        phot[key] = tar.extractfile(member).read()
            self._phot = phot
        return self._phot

    def light_curve(self, event_id: str) -> Table:
        raw = self._load_phot().get(event_id)
        if raw is None:
            raise KeyError(f"{self.name}: no photometry for {event_id!r}")
        hjd, mag, err = parse_photometry(raw)
        return standard_light_curve(
            hjd,
            mag,
            err,
            "I",
            source=f"{self.spec.base}{self.spec.phot}:{event_id}.dat ({self.spec.reference})",
            time_system="HJD",
        )

    def _load_eff(self) -> dict[str, tuple]:
        if self._eff is None:
            eff = {}
            with tarfile.open(self.path(self.spec.eff), "r:gz") as tar:
                for member in tar:
                    if member.isfile():
                        key = member.name.rsplit("/", 1)[-1].split(".")[0]
                        eff[key] = parse_efficiency(tar.extractfile(member).read().decode())
            self._eff = eff
        return self._eff

    def fields(self) -> Table:
        """Fields of this sample with ``n_sources`` (monitored sources) and efficiency coverage.

        bulge2019: N_s of sources brighter than I = 21 (Table 7), the nine high-cadence fields
        excluded. disk2020: number of stars in the database (Table A1, millions); Mróz et al.
        2020 count sources brighter than I = 21 separately, so this is an ASSUMPTION (≈ N_s).
        """
        if self.spec.key == "bulge2019":
            t = parse_byte_table(self._text(self.spec.fields))
            keep = np.array([f not in HIGH_CADENCE_FIELDS for f in t["field"]])
            out = Table({"field": t["field"][keep], "n_sources": t["N_stars"][keep] * 1e6})
            note = "N_s (I < 21) from Table 7"
        else:
            t = parse_whitespace_table(
                self._text(self.spec.fields),
                ["field", "ra", "dec", "glon", "glat", "n_stars", "n_epochs", "delta_t"],
            )
            out = Table({"field": t["field"], "n_sources": t["n_stars"] * 1e6})
            note = "database stars from Table A1 (ASSUMPTION: equal to N_s)"
        eff = self._load_eff()
        out["has_efficiency"] = [f in eff for f in out["field"]]
        out.meta.update(provenance=schema.Provenance.OBSERVED.value, source=note)
        return out

    def field_efficiency(self, field: str, t_e_days) -> np.ndarray:
        lo, hi, eff = self._load_eff()[field]
        return eff_at(lo, hi, eff, t_e_days)

    def efficiency(self, t_e_days: float) -> float | None:
        """N_s-weighted mean published efficiency over the sample's fields at ``t_e_days``."""
        f = self.fields()
        f = f[f["has_efficiency"]]
        w = np.asarray(f["n_sources"], float)
        e = np.array([self.field_efficiency(x, t_e_days)[0] for x in f["field"]])
        return float(np.sum(w * e) / np.sum(w))

    def exposure_star_years(self, t_e_days: float) -> float:
        """Σ_fields N_s ε(t_E) ΔT in star-years (the denominator of the published event rate)."""
        f = self.fields()
        f = f[f["has_efficiency"]]
        e = np.array([self.field_efficiency(x, t_e_days)[0] for x in f["field"]])
        return float(
            np.sum(np.asarray(f["n_sources"], float) * e) * self.spec.delta_t_days / 365.25
        )
