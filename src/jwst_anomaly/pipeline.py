"""Config-driven end-to-end runner chaining the stage modules. Owner: bootstrap unit 6.

For each sample in the YAML config (e.g. ``configs/reference_sample.yaml``)::

    query.query_observations -> query.list_products -> acquire.fetch_products (CAT only)
    -> catalog.load_pipeline_catalog (per band) -> catalog.merge_bands
    -> features.build_features -> quality.assess_sources (optional gate, D-011)
    -> classify.classify_sources (optional, D-012: stars ranked as stratum "<sample>-stars")
    -> rank.score_anomalies (quality_ok sources, per stratum) -> top-k targets
    -> cutouts.make_cutouts (optional, per band, S3 i2d URI) -> crossmatch.crossmatch (optional)

then the top-ranked candidates go into the :class:`~jwst_anomaly.candidates.CandidateStore`
and a ``report.md`` is written. Every stage output is checked with ``schema.validate``.
A failing required stage aborts the run (recorded as ``failed``); a failing optional stage
(quality gate, cutouts, crossmatch) is recorded and the run continues (an unusable gate means
all sources are ranked, marked "ungated"). Image products (``_i2d.fits``,
~1.8 GB) are never downloaded here: only catalogs are fetched, cutouts read from S3.

Config keys read: ``name``, ``archive.{collection, calib_level, data_rights,
product_subgroups}``, ``samples[].{id, role, description, proposal_id, instrument_name,
ref_band, obs_ids, query}`` (``query``: extra MAST criteria, optional),
``cloud.{s3_bucket, l3_key_pattern}``, ``stages.catalog.merge_radius_arcsec``,
``stages.quality.{enabled, grid_arcsec, step, min_rel_weight, min_edge_arcsec, max_artifact_ci,
min_ranked, min_detection_snr, require_multiband}`` (gate skipped when the block is absent),
``stages.classify.{enabled, services, radius_arcsec, gaia_radius_arcsec, star_top_k, min_stars,
stellar_locus}``
(skipped when absent), ``stages.features.daofind_max_ci``,
``samples[].matched_photometry.{url, sha256, label, aperture, radius_arcsec, max_bytes}`` (D-013),
``stages.rank.{methods, random_state}``,
``stages.cutouts.{enabled, top_k, size_arcsec, bands, spike, screen_low_weight}``
(spike: D-018-020; screen_low_weight: D-021),
``stages.crossmatch.{enabled, top_k, radius_arcsec, services}`` and the optional top-level
``outputs`` block added by this unit:

* ``outputs.candidates_top_k`` -- candidates stored/reported per sample
  (default: the larger cutouts/crossmatch ``top_k``, else 20);
* ``outputs.table_format`` -- ``ecsv`` (default) or ``parquet`` for intermediate tables.

Run outputs live in ``<outputs>/runs/<run_id>/`` (``paths.outputs_dir()`` by default):
``config.yaml`` (verbatim copy), ``run_context.json``, ``run_record.json``, ``report.md`` and
``<sample_id>/<table>.<ext>``. The store defaults to ``<outputs>/candidates.sqlite``.
"""

from __future__ import annotations

import json
import logging
import re
import shutil
import time
from collections.abc import Callable, Mapping, Sequence
from dataclasses import asdict, dataclass, field
from functools import partial
from pathlib import Path
from typing import Any

import numpy as np
import yaml
from astropy.table import Table, vstack

from jwst_anomaly import (
    acquire,
    catalog,
    classify,
    crossmatch,
    cutouts,
    features,
    paths,
    photometry,
    provenance,
    quality,
    query,
    rank,
    schema,
    viz,
)
from jwst_anomaly.candidates import CandidateStore, to_python
from jwst_anomaly.features import column_as_float

log = logging.getLogger(__name__)

DEFAULT_DB_NAME = "candidates.sqlite"
DEFAULT_TOP_K = 20
TABLE_FORMATS = {"ecsv": "ascii.ecsv", "parquet": "parquet"}

DISCLAIMER = (
    "Anomaly scores rank sources by how unusual they are *relative to this sample, under this "
    "feature set and these models* (provenance: `model_prediction`). A high score is not a "
    "discovery and is not evidence of new physics. Every candidate below has status `new` "
    "(not vetted). Expected explanations, in order: processing or instrument artifacts, catalog "
    "effects, known but rare astrophysical populations; only then 'unexplained under tests X, "
    "Y, Z' (docs/methodology.md)."
)

LIMITATIONS = (
    "Catalog-first: JWST level-3 pipeline catalogs detect each band independently (no forced "
    "photometry), so cross-band colors are approximate and cross-band mismatches can create "
    "spurious outliers (DECISIONS.md D-001).",
    "Each sample is ranked on its own; scores and ranks are not comparable across samples, "
    "feature sets, methods or runs with different configs.",
    "Image quality flags are automatic heuristics on cutouts; they do not replace visual "
    "inspection.",
    "Cross-match: no counterpart within the search radius does not make a source unknown "
    "(catalog depth and completeness vary); a counterpart does not prove association.",
    "Optional stages that failed or were skipped are listed under Stages; their columns read "
    "'not run'.",
)

_BAND_RE = re.compile(r"f\d{3,4}[wmnc]\d?", re.IGNORECASE)
_OBS_ID_RE = re.compile(r"^jw(\d{5})-o(\d{3})_")
_SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")  # used as a directory name


class ConfigError(ValueError):
    """The run configuration is missing or malformed."""


class PipelineError(RuntimeError):
    """A required stage failed; the run was recorded as ``failed``."""

    def __init__(
        self, message: str, *, run_id: str | None = None, report_path: Path | None = None
    ) -> None:
        super().__init__(message)
        self.run_id = run_id
        self.report_path = report_path


class _StageFailure(Exception):
    def __init__(self, sample_id: str, stage: str, message: str) -> None:
        super().__init__(f"required stage {stage} failed for sample {sample_id!r}: {message}")


@dataclass
class StageRecord:
    """Outcome of one stage call, as stored in the run record and report."""

    sample: str
    stage: str
    required: bool
    status: str  # ok | failed | skipped
    n_rows: int | None = None
    seconds: float = 0.0
    message: str = ""
    output: str | None = None


@dataclass
class SampleSummary:
    """Per-sample facts for the report (counts, bands, methods)."""

    id: str
    role: str
    description: str
    ref_band: str
    n_observations: int = 0
    bands: list[str] = field(default_factory=list)
    n_sources: int = 0
    feature_spec: dict[str, str] = field(default_factory=dict)
    methods: list[str] = field(default_factory=list)
    n_scored: int = 0
    cutout_bands: list[str] = field(default_factory=list)
    tables: list[str] = field(default_factory=list)  # saved files, relative to the run dir
    contact_sheet: str | None = None  # PNG of the top-k cutouts, relative to the run dir
    quality: dict[str, int] = field(default_factory=dict)  # quality-gate counts (D-011)
    quality_text: str = ""  # replaces the counts line (e.g. strata of a gated sample)
    topk: dict[str, int] = field(default_factory=dict)  # top-k composition (tracked metric)
    notes: list[str] = field(default_factory=list)


# -- config ------------------------------------------------------------------------------------


def _check_positive_int(value: Any, where: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ConfigError(f"{where} must be a positive integer, got {value!r}")


def load_config(path: str | Path) -> dict[str, Any]:
    """Read and minimally validate a run config (see module docstring for the keys)."""
    path = Path(path)
    if not path.is_file():
        raise ConfigError(f"config file not found: {path}")
    try:
        config = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise ConfigError(f"{path}: invalid YAML: {exc}") from exc
    if not isinstance(config, dict):
        raise ConfigError(f"{path}: top level must be a mapping")
    if not isinstance(config.get("name"), str) or not config["name"]:
        raise ConfigError(f"{path}: 'name' must be a non-empty string")
    samples = config.get("samples")
    if not isinstance(samples, list) or not samples:
        raise ConfigError(f"{path}: 'samples' must be a non-empty list")
    seen: set[str] = set()
    for i, sample in enumerate(samples):
        if not isinstance(sample, dict) or not sample.get("id"):
            raise ConfigError(f"{path}: samples[{i}] needs an 'id'")
        sid = str(sample["id"])
        if not _SAMPLE_ID_RE.fullmatch(sid):
            raise ConfigError(f"{path}: sample id {sid!r} may only use letters, digits, _ . -")
        if sid in seen:
            raise ConfigError(f"{path}: duplicate sample id {sid!r}")
        seen.add(sid)
        if not sample.get("ref_band"):
            raise ConfigError(f"{path}: sample {sid!r} needs a 'ref_band'")
        if not sample.get("obs_ids") and not sample.get("proposal_id"):
            raise ConfigError(f"{path}: sample {sid!r} needs 'obs_ids' or 'proposal_id'")
        if sample.get("matched_photometry") is not None:
            _check_matched_photometry(sample["matched_photometry"], f"{path}: sample {sid!r}")
    for key in ("archive", "stages", "cloud", "outputs"):
        if config.get(key) is not None and not isinstance(config[key], dict):
            raise ConfigError(f"{path}: '{key}' must be a mapping")
    for name, block in (config.get("stages") or {}).items():
        if block is None:
            continue
        if not isinstance(block, dict):
            raise ConfigError(f"{path}: stages.{name} must be a mapping")
        if "enabled" in block and not isinstance(block["enabled"], bool):
            raise ConfigError(f"{path}: stages.{name}.enabled must be true or false")
        if "top_k" in block:
            _check_positive_int(block["top_k"], f"{path}: stages.{name}.top_k")
    qcfg = (config.get("stages") or {}).get("quality") or {}
    if "require_multiband" in qcfg and not isinstance(qcfg["require_multiband"], bool):
        raise ConfigError(f"{path}: stages.quality.require_multiband must be true or false")
    for key in ("min_detection_snr", "min_rel_weight", "min_edge_arcsec", "max_artifact_ci"):
        value = qcfg.get(key)
        if value is not None and not (
            isinstance(value, int | float) and not isinstance(value, bool) and value > 0
        ):
            raise ConfigError(f"{path}: stages.quality.{key} must be a positive number")
    locus = ((config.get("stages") or {}).get("classify") or {}).get("stellar_locus")
    if locus is not None and not isinstance(locus, bool | dict):
        raise ConfigError(f"{path}: stages.classify.stellar_locus must be a mapping or true/false")
    if isinstance(locus, dict):
        unknown = set(locus) - set(classify.DEFAULT_LOCUS) - {"enabled"}
        if unknown:
            raise ConfigError(
                f"{path}: unknown stages.classify.stellar_locus keys {sorted(unknown)}"
            )
    veto = ((config.get("stages") or {}).get("classify") or {}).get("extended_veto")
    if veto is not None and not isinstance(veto, bool | dict):
        raise ConfigError(f"{path}: stages.classify.extended_veto must be a mapping or true/false")
    if isinstance(veto, dict):
        unknown = set(veto) - set(classify.DEFAULT_EXTENDED_VETO) - {"enabled"}
        if unknown:
            raise ConfigError(
                f"{path}: unknown stages.classify.extended_veto keys {sorted(unknown)}"
            )
    screen_lw = ((config.get("stages") or {}).get("cutouts") or {}).get("screen_low_weight")
    if screen_lw is not None and not isinstance(screen_lw, bool):
        raise ConfigError(f"{path}: stages.cutouts.screen_low_weight must be true or false")
    spike = ((config.get("stages") or {}).get("cutouts") or {}).get("spike")
    if spike is not None:
        _check_spike(spike, f"{path}: stages.cutouts.spike")
    daofind = ((config.get("stages") or {}).get("features") or {}).get("daofind_max_ci")
    if daofind is not None and not (
        isinstance(daofind, int | float) and not isinstance(daofind, bool) and daofind > 0
    ):
        raise ConfigError(f"{path}: stages.features.daofind_max_ci must be a positive number")
    outputs = config.get("outputs") or {}
    if "candidates_top_k" in outputs:
        _check_positive_int(outputs["candidates_top_k"], f"{path}: outputs.candidates_top_k")
    if outputs.get("table_format", "ecsv") not in TABLE_FORMATS:
        raise ConfigError(f"{path}: outputs.table_format must be one of {sorted(TABLE_FORMATS)}")
    return config


def _check_matched_photometry(cfg: Any, where: str) -> None:
    """Validate a sample's ``matched_photometry`` block (D-013)."""
    if not isinstance(cfg, dict):
        raise ConfigError(f"{where}: matched_photometry must be a mapping")
    url = cfg.get("url")
    if not isinstance(url, str) or not url.startswith(("https://", "http://")):
        raise ConfigError(f"{where}: matched_photometry.url must be an http(s) URL")
    sha = cfg.get("sha256")
    if not (isinstance(sha, str) and re.fullmatch(r"[0-9a-fA-F]{64}", sha)):
        raise ConfigError(f"{where}: matched_photometry.sha256 must be 64 hex characters")
    label = cfg.get("label")
    if not (isinstance(label, str) and label.isalnum() and label[:1].isalpha()):
        raise ConfigError(f"{where}: matched_photometry.label must be alphanumeric, e.g. dja05")
    for key in ("aperture", "max_bytes"):
        if key in cfg:
            if key == "aperture" and cfg[key] == 0:
                continue
            _check_positive_int(cfg[key], f"{where}: matched_photometry.{key}")
    radius = cfg.get("radius_arcsec")
    if radius is not None and not (
        isinstance(radius, int | float) and not isinstance(radius, bool) and radius > 0
    ):
        raise ConfigError(f"{where}: matched_photometry.radius_arcsec must be positive")


def default_db_path(outputs_dir: str | Path | None = None) -> Path:
    """Default candidate store: ``<outputs>/candidates.sqlite``."""
    base = Path(outputs_dir) if outputs_dir is not None else paths.outputs_dir()
    return base / DEFAULT_DB_NAME


def band_from_name(name: str) -> str | None:
    """Last JWST filter name in ``name`` (``...nircam_clear-f200w_cat.ecsv`` -> ``F200W``)."""
    matches = _BAND_RE.findall(name)
    return matches[-1].upper() if matches else None


def l3_image_uri(cloud: Mapping[str, Any] | None, obs_id: str, suffix: str = "i2d.fits") -> str:
    """S3 URI of a level-3 product from ``cloud.s3_bucket`` and ``cloud.l3_key_pattern``."""
    if not cloud or not cloud.get("s3_bucket") or not cloud.get("l3_key_pattern"):
        raise ConfigError("config 'cloud' needs 's3_bucket' and 'l3_key_pattern' for cutouts")
    match = _OBS_ID_RE.match(obs_id)
    if not match:
        raise ValueError(f"cannot parse program/observation from obs_id {obs_id!r}")
    key = cloud["l3_key_pattern"].format(
        proposal=int(match[1]), observation=int(match[2]), obs_id=obs_id, suffix=suffix
    )
    return f"s3://{cloud['s3_bucket']}/{key}"


def _band_of(table: Table, *names: str) -> str | None:
    """Band of a catalog: explicit ``meta['band']``, else the obs_id/file name, else meta filter.

    The name wins over a FILTER keyword because NIRCam pupil-wheel filters (``f444w-f470n``)
    are named last in obs ids, while FILTER holds the filter-wheel element only.
    """
    candidates: list[Any] = [table.meta.get("band"), table.meta.get("BAND")]
    candidates += [band_from_name(n) for n in names]
    candidates += [table.meta.get("filter"), table.meta.get("FILTER")]
    for value in candidates:
        if isinstance(value, str) and _BAND_RE.fullmatch(value.strip()):
            return value.strip().upper()
    return None


def _describe(exc: BaseException) -> str:
    if isinstance(exc, NotImplementedError):
        return f"stage not implemented yet ({exc})"
    return f"{type(exc).__name__}: {exc}"


def _subset(table: Table, mask: np.ndarray, what: str) -> Table:
    """Rows of ``table`` under ``mask``, with ``meta['source']`` naming the selection."""
    out = table[mask]
    out.meta = dict(table.meta)
    out.meta["source"] = (
        f"{what}: {int(mask.sum())} of {len(table)} rows of: {table.meta.get('source', '')}"
    )
    return out


def _topk_metrics(
    cands: Sequence[Mapping[str, Any]],
    cutout_rows: Mapping[str, Mapping[str, Mapping[str, Any]]],
    xmatch_rows: Mapping[str, Mapping[str, Any]],
) -> dict[str, int]:
    """Composition of a stratum's top k: a tracked contamination and recovery metric.

    Counts candidates whose cutouts carry any quality flag, the ``spikes`` flag (D-018), a
    cross-match to a known object, to a star, or to a lens-related object. Candidates without
    cutouts or cross-match are counted as ``no_cutout`` / ``no_xmatch``.
    """
    out = dict.fromkeys(
        ("n", "flagged", "spikes", "no_cutout", "known", "star", "lens_related", "no_xmatch"), 0
    )
    for rec in cands:
        uid = rec["source_uid"]
        out["n"] += 1
        tokens = {
            t.strip()
            for q in (cutout_rows.get(uid) or {}).values()
            for t in str(q.get("quality_flag", "")).split(",")
        }
        if uid not in cutout_rows:
            out["no_cutout"] += 1
        elif tokens - {"ok", ""}:
            out["flagged"] += 1
        out["spikes"] += "spikes" in tokens
        xm = xmatch_rows.get(uid)
        if xm is None:
            out["no_xmatch"] += 1
            continue
        out["known"] += bool(xm.get("is_known_object"))
        out["star"] += bool(xm.get("is_star"))
        out["lens_related"] += bool(xm.get("is_lens_related"))
    return out


def _topk_line(m: Mapping[str, int]) -> str:
    """Report text; each fraction uses the candidates that have that evidence (cut out or
    cross-matched), since the evidence stages may cover fewer than the stored top k."""
    n = m.get("n", 0)
    if not n:
        screened = m.get("screened", 0)
        return f"none ({screened} screened out, D-019 to D-021)" if screened else "none"
    n_xm = n - m.get("no_xmatch", 0)
    n_cut = n - m.get("no_cutout", 0)

    def part(key: str, label: str, denom: int) -> str:
        frac = f" ({m.get(key, 0) / denom:.0%})" if denom else ""
        return f"{label} {m.get(key, 0)}{frac}"

    text = ", ".join(
        [
            part("known", "known object", n_xm),
            part("lens_related", "lens-related", n_xm),
            part("star", "catalogued star", n_xm),
            part("flagged", "cutout-flagged", n_cut),
            part("spikes", "spikes", n_cut),
        ]
    )
    extra = [f"{n_xm} cross-matched, {n_cut} cut out"] if (n_xm, n_cut) != (n, n) else []
    if m.get("screened"):
        extra.append(f"{m['screened']} screened out before selection (D-019 to D-021)")
    return f"n = {n}: {text}" + (f"; {', '.join(extra)}" if extra else "")


def _has_flag(per_band: Mapping[str, Mapping[str, Any]], token: str) -> bool:
    """Whether any band's cutout ``quality_flag`` carries ``token``."""
    return any(
        token in {t.strip() for t in str(q.get("quality_flag", "")).split(",")}
        for q in per_band.values()
    )


def _spike_bands(per_band: Mapping[str, Mapping[str, Any]]) -> dict[str, Mapping[str, Any]]:
    """The bands whose cutout carries the ``spikes`` flag (the host test runs only there)."""
    return {
        band: q
        for band, q in per_band.items()
        if "spikes" in {t.strip() for t in str(q.get("quality_flag", "")).split(",")}
    }


def _finite_extreme(
    per_band: Mapping[str, Mapping[str, Any]], key: str, pick: Any = max
) -> float | None:
    """``pick`` (max or min) of the finite ``key`` values over a source's bands; None if none."""
    values = [_float(q.get(key)) for q in per_band.values()]
    finite = [v for v in values if v is not None and np.isfinite(v)]
    return pick(finite) if finite else None


def _check_spike(spike: Any, where: str) -> None:
    """``stages.cutouts.spike`` (D-018): ``radii_arcsec`` [r_in, r_out], ``threshold``,
    optional ``search_arcsec`` and ``screen`` (D-019)."""
    if not isinstance(spike, dict):
        raise ConfigError(f"{where} must be a mapping")
    unknown = set(spike) - {
        "radii_arcsec",
        "threshold",
        "search_arcsec",
        "screen",
        "host_annulus_arcsec",
        "host_ratio_max",
        "host_exempt_colour",
    }
    if unknown:
        raise ConfigError(f"{where}: unknown keys {sorted(unknown)}")

    def number(v: Any) -> bool:
        return isinstance(v, int | float) and not isinstance(v, bool)

    radii = spike.get("radii_arcsec")
    if not (
        isinstance(radii, list | tuple)
        and len(radii) == 2
        and all(number(r) for r in radii)
        and 0 <= radii[0] < radii[1]
    ):
        raise ConfigError(f"{where}.radii_arcsec must be [r_in, r_out] with 0 <= r_in < r_out")
    if not (number(spike.get("threshold")) and spike["threshold"] > 0):
        raise ConfigError(f"{where}.threshold must be a positive number")
    if "search_arcsec" in spike and not (
        number(spike["search_arcsec"]) and spike["search_arcsec"] >= 0
    ):
        raise ConfigError(f"{where}.search_arcsec must be a non-negative number")
    if "screen" in spike and not isinstance(spike["screen"], bool):
        raise ConfigError(f"{where}.screen must be true or false")
    ann = spike.get("host_annulus_arcsec")
    if ann is not None and not (
        isinstance(ann, list | tuple)
        and len(ann) == 2
        and all(number(r) for r in ann)
        and 0 <= ann[0] < ann[1]
    ):
        raise ConfigError(
            f"{where}.host_annulus_arcsec must be [r_in, r_out] with 0 <= r_in < r_out"
        )
    hec = spike.get("host_exempt_colour")
    if hec is not None and not (
        isinstance(hec, list | tuple)
        and len(hec) == 3
        and number(hec[2])
        and all(isinstance(b, str) for b in hec[:2])
    ):
        raise ConfigError(f"{where}.host_exempt_colour must be [blue_band, red_band, min_colour]")
    if "host_ratio_max" in spike and not (
        number(spike["host_ratio_max"]) and spike["host_ratio_max"] > 0
    ):
        raise ConfigError(f"{where}.host_ratio_max must be a positive number")


def _stage_kwargs(cfg: Mapping[str, Any], keys: Sequence[str]) -> dict[str, Any]:
    """Stage parameters present in the config; lists become tuples (stage defaults otherwise)."""
    return {k: tuple(cfg[k]) if isinstance(cfg[k], list) else cfg[k] for k in keys if k in cfg}


def _text(value: Any) -> str:
    """Table cell as ``str`` (bytes decoded; masked -> empty string)."""
    value = to_python(value)
    return "" if value is None else str(value)


def _float(value: Any) -> float | None:
    value = to_python(value)
    return float(value) if isinstance(value, int | float) and not isinstance(value, bool) else None


# -- runner ------------------------------------------------------------------------------------


class _Runner:
    def __init__(
        self,
        config: dict[str, Any],
        run_id: str,
        run_dir: Path,
        manifest_path: Path | None = None,
    ) -> None:
        self.config = config
        # Shared with scripts/fetch_reference_sample.py: one tracked manifest per config file.
        self.manifest_path = manifest_path
        self.run_id = run_id
        self.run_dir = run_dir
        self.stages_cfg: dict[str, Any] = config.get("stages") or {}
        outputs = config.get("outputs") or {}
        self.table_format: str = outputs.get("table_format", "ecsv")
        cut_cfg = self.stages_cfg.get("cutouts") or {}
        xm_cfg = self.stages_cfg.get("crossmatch") or {}
        self.cutouts_top_k = int(cut_cfg.get("top_k", DEFAULT_TOP_K))
        self.crossmatch_top_k = int(xm_cfg.get("top_k", DEFAULT_TOP_K))
        default_k = max(self.cutouts_top_k, self.crossmatch_top_k)
        self.candidates_top_k = int(outputs.get("candidates_top_k", default_k))
        self.records: list[StageRecord] = []
        # Reference weight of each sample's gate map, by (sample id, image URI), shared with
        # cutouts of that image in the sample's strata.
        self.weight_refs: dict[tuple[str, str], float] = {}
        # D-019: per sample, the source_uids with stellar matched-photometry colours (D-015 box).
        self.stellar_colour_uids: dict[str, set[str]] = {}
        # D-025: the matched-photometry sources (and label) for colour checks during screening.
        self.matched_sources: dict[str, tuple[Table, str]] = {}
        self.summaries: dict[str, SampleSummary] = {}
        self.candidates: list[dict[str, Any]] = []

    def settings(self) -> dict[str, Any]:
        return {
            "candidates_top_k": self.candidates_top_k,
            "cutouts_top_k": self.cutouts_top_k,
            "crossmatch_top_k": self.crossmatch_top_k,
            "table_format": self.table_format,
        }

    # -- stage plumbing ------------------------------------------------------------------------

    def stage(
        self,
        sample_id: str,
        name: str,
        call: Callable[[], Any],
        *,
        columns: Sequence[str],
        required: bool,
        save_as: str | None = None,
    ) -> Table | None:
        """Run one stage call, validate its table, persist it and record the outcome."""
        t0 = time.perf_counter()
        try:
            table = call()
            if not isinstance(table, Table):
                raise TypeError(f"returned {type(table).__name__}, expected astropy Table")
            schema.validate(table, columns, name=name)
        except Exception as exc:  # stage code is outside this unit; report any failure
            message = _describe(exc)
            self.records.append(
                StageRecord(
                    sample_id, name, required, "failed", None, time.perf_counter() - t0, message
                )
            )
            if required:
                raise _StageFailure(sample_id, name, message) from exc
            log.warning("[%s] optional stage %s failed: %s", sample_id, name, message)
            return None
        output = self.save(table, sample_id, save_as) if save_as else None
        record = StageRecord(
            sample_id, name, required, "ok", len(table), time.perf_counter() - t0, "", output
        )
        self.records.append(record)
        log.info("[%s] %s: ok (%d rows, %.2f s)", sample_id, name, len(table), record.seconds)
        return table

    def record(
        self, sample_id: str, name: str, status: str, message: str, *, required: bool
    ) -> None:
        self.records.append(StageRecord(sample_id, name, required, status, message=message))
        log.info("[%s] %s: %s %s", sample_id, name, status, message)

    def fail(self, sample_id: str, name: str, message: str) -> None:
        self.record(sample_id, name, "failed", message, required=True)
        raise _StageFailure(sample_id, name, message)

    def save(self, table: Table, sample_id: str, name: str) -> str | None:
        """Persist an intermediate table (best effort); return its path relative to the run dir."""
        rel = f"{sample_id}/{name}.{self.table_format}"
        summary = self.summaries[sample_id]
        try:
            path = self.run_dir / rel
            path.parent.mkdir(parents=True, exist_ok=True)
            out = table.copy(copy_data=False)
            out.meta = to_python(dict(table.meta))  # plain types, so any writer can serialize
            out.write(path, format=TABLE_FORMATS[self.table_format], overwrite=True)
        except Exception as exc:  # a missing intermediate file must not abort the run
            log.warning("[%s] could not write %s: %s", sample_id, rel, exc)
            summary.notes.append(f"could not write {rel}: {exc}")
            return None
        summary.tables.append(rel)
        return rel

    # -- one sample ----------------------------------------------------------------------------

    def run_sample(self, sample: Mapping[str, Any]) -> None:
        sid = str(sample["id"])
        ref_band = str(sample["ref_band"]).upper()
        summary = SampleSummary(
            sid, str(sample.get("role", "science")), str(sample.get("description", "")), ref_band
        )
        self.summaries[sid] = summary
        archive = self.config.get("archive") or {}

        # 1. Observations (observed). Public only, and only the configured obs_ids.
        criteria: dict[str, Any] = {}
        if archive.get("collection"):
            criteria["obs_collection"] = archive["collection"]
        for key in ("proposal_id", "instrument_name"):
            if sample.get(key):
                criteria[key] = str(sample[key])
        if archive.get("calib_level") is not None:
            criteria["calib_level"] = archive["calib_level"]
        obs_ids = [str(o) for o in sample.get("obs_ids") or []]
        if obs_ids:
            criteria["obs_id"] = obs_ids
        criteria.update(sample.get("query") or {})
        obs = self.stage(
            sid,
            "query.query_observations",
            partial(query.query_observations, **criteria),
            columns=schema.OBSERVATION_COLUMNS,
            required=True,
        )
        obs = self._filter_observations(sid, obs, obs_ids, archive.get("data_rights"))
        summary.n_observations = len(obs)
        self.save(obs, sid, "observations")

        # 2. Product list, then fetch catalogs only (images are read remotely for cutouts).
        product_kwargs = {"calib_level": archive.get("calib_level", 3)}
        if archive.get("product_subgroups"):
            product_kwargs["subgroups"] = tuple(archive["product_subgroups"])
        products = self.stage(
            sid,
            "query.list_products",
            partial(query.list_products, obs, **product_kwargs),
            columns=schema.PRODUCT_COLUMNS,
            required=True,
            save_as="products",
        )
        subgroup = np.array([_text(s).upper() for s in products["productSubGroupDescription"]])
        cat_products = products[subgroup == "CAT"]
        if len(cat_products) == 0:
            self.fail(sid, "select catalog products", "no CAT products in the product list")
        manifest = self.stage(
            sid,
            "acquire.fetch_products",
            partial(acquire.fetch_products, cat_products, manifest_path=self.manifest_path),
            columns=schema.MANIFEST_COLUMNS,
            required=True,
            save_as="manifest",
        )

        # 3. Per-band pipeline catalogs (observed) -> merged sources (derived).
        uri_to_obs = {_text(r["dataURI"]): _text(r["obs_id"]) for r in products}
        catalogs: dict[str, Table] = {}
        band_obs: dict[str, str] = {}
        for row in manifest:
            filename = _text(row["productFilename"])
            obs_id = uri_to_obs.get(_text(row["dataURI"]), "")
            local = paths.data_root() / _text(row["local_path"])
            table = self.stage(
                sid,
                f"catalog.load_pipeline_catalog[{filename}]",
                partial(catalog.load_pipeline_catalog, local),
                columns=schema.BAND_CATALOG_COLUMNS,
                required=True,
            )
            band = _band_of(table, obs_id, filename)
            if band is None:
                self.fail(sid, "catalog band", f"cannot determine the band of {filename}")
            if band in catalogs:
                self.fail(sid, "catalog band", f"two catalogs for band {band} ({filename})")
            catalogs[band] = table
            band_obs[band] = obs_id or filename.rsplit("_", 1)[0]
        summary.bands = sorted(catalogs)
        if ref_band not in catalogs:
            self.fail(sid, "catalog band", f"ref_band {ref_band} not among {summary.bands}")
        cat_cfg = self.stages_cfg.get("catalog") or {}
        merge_kwargs = {}
        if "merge_radius_arcsec" in cat_cfg:
            merge_kwargs["radius_arcsec"] = float(cat_cfg["merge_radius_arcsec"])
        sources = self.stage(
            sid,
            "catalog.merge_bands",
            partial(catalog.merge_bands, catalogs, ref_band=ref_band, **merge_kwargs),
            columns=schema.SOURCE_COLUMNS,
            required=True,
            save_as="sources",
        )
        summary.n_sources = len(sources)

        # 4. Features (derived) -> scores (model_prediction). Colours come from matched-aperture
        # photometry when the sample configures it (D-013), else from the pipeline catalogs.
        feat_sources, aperture = self._matched_photometry(sid, sample, sources)
        feat_kwargs: dict[str, Any] = {}
        if aperture:
            feat_kwargs["aperture"] = aperture
        fcfg = self.stages_cfg.get("features") or {}
        if fcfg.get("daofind_max_ci") is not None:
            feat_kwargs["daofind_max_ci"] = float(fcfg["daofind_max_ci"])
        feats = self.stage(
            sid,
            "features.build_features",
            partial(features.build_features, feat_sources, **feat_kwargs),
            columns=schema.FEATURE_ID_COLUMNS,
            required=True,
            save_as="features",
        )
        summary.feature_spec = {
            str(k): str(v) for k, v in (feats.meta.get("feature_spec") or {}).items()
        }
        # 4b. Quality gate (derived, D-011): rank only sources whose measurements can be trusted.
        to_rank = feats
        confirm = photometry.match_sep_column(aperture) if aperture else None
        keep = self._quality(sid, summary, feat_sources, feats, ref_band, band_obs, confirm)
        if keep is not None:
            to_rank = _subset(feats, keep, "quality-gated subset (D-011)")
        # 4c. Star/galaxy separation (D-012): stars become their own ranking stratum.
        split = self._classify(sid, feat_sources, to_rank, aperture)
        if split is not None:
            stars, star_top_k, min_stars = split
            n_stars, n_other = int(stars.sum()), int((~stars).sum())
            if n_stars < min_stars or n_other < 2:
                summary.notes.append(
                    f"star/galaxy strata not used ({n_stars} stars, {n_other} others; "
                    f"min_stars {min_stars}): one ranking"
                )
                split = None
        if split is None:
            self._rank_stratum(sid, summary, to_rank, sources, band_obs, required=True)
            return
        star_sid = f"{sid}-stars"
        summary.notes.append(f"{n_stars} stars ranked separately as {star_sid} (D-012)")
        galaxies = _subset(to_rank, ~stars, "non-star stratum (D-012)")
        self._rank_stratum(
            sid, summary, galaxies, sources, band_obs, required=True, galaxy_stratum=True
        )
        star_summary = SampleSummary(
            star_sid,
            f"{summary.role}/stars",
            f"stars of {sid} (Gaia/SIMBAD; stellar locus D-015 where applied), ranked among "
            "themselves (D-012)",
            summary.ref_band,
            n_observations=summary.n_observations,
            bands=list(summary.bands),
            n_sources=n_stars,
            feature_spec=dict(summary.feature_spec),
            quality_text=f"applied in {sid}; this stratum holds its {n_stars} passing stars",
        )
        self.summaries[star_sid] = star_summary
        self._rank_stratum(
            star_sid,
            star_summary,
            _subset(to_rank, stars, "star stratum (D-012)"),
            sources,
            band_obs,
            required=False,
            top_k=star_top_k,
            parent=sid,
        )

    def _rank_stratum(
        self,
        label: str,
        summary: SampleSummary,
        to_rank: Table,
        sources: Table,
        band_obs: Mapping[str, str],
        *,
        required: bool,
        top_k: int | None = None,
        parent: str | None = None,
        galaxy_stratum: bool = False,
    ) -> None:
        """Score one stratum, then select, cut out, cross-match and store its top k.

        Spike screening (D-019) runs only for ``galaxy_stratum`` (a successful D-012 split),
        never for a star stratum or a single mixed ranking.

        ``parent`` is the sample a stratum (``<sample>-stars``) belongs to; its gate's
        reference weights are reused for the stratum's cutouts.
        """
        rank_kwargs = _stage_kwargs(self.stages_cfg.get("rank") or {}, ("methods", "random_state"))
        scores = self.stage(
            label,
            "rank.score_anomalies",
            partial(rank.score_anomalies, to_rank, **rank_kwargs),
            columns=schema.SCORE_COLUMNS,
            required=required,
            save_as="scores",
        )
        if scores is None:
            return
        summary.n_scored = len(scores)
        summary.methods = [c[len("score_") :] for c in scores.colnames if c.startswith("score_")]

        def cap(k: int) -> int:
            return k if top_k is None else min(k, top_k)

        cand_k = cap(self.candidates_top_k)
        cut_k = cap(self.cutouts_top_k)
        xm_k = cap(self.crossmatch_top_k)

        # 5. Top-k targets, joined to positions.
        ranked = self._ranked(label, scores, sources)
        n_top = max(cand_k, cut_k, xm_k)

        # 6. Optional evidence stages. With screening, cut out 2 x top_k, drop sources the
        # cutouts show to be artifacts or stars and keep the next clean ones; the dropped sources
        # are saved as "screened" with a reason. Spike screening (D-019/D-020) runs on galaxy
        # strata only; low-weight screening (D-021) on any stratum that is not a star stratum.
        cut_cfg = self.stages_cfg.get("cutouts") or {}
        spike_cfg = cut_cfg.get("spike") or {}
        spike_screen = bool(spike_cfg.get("screen")) and galaxy_stratum
        weight_screen = bool(cut_cfg.get("screen_low_weight")) and parent is None
        screen = spike_screen or weight_screen
        pool = self._targets(label, ranked[: (2 * cut_k if screen else cut_k)])
        cutout_rows, cut_table = self._cutouts(
            label, summary, pool, band_obs, parent=parent or label, render=not screen
        )
        if screen:
            removed: list[dict[str, Any]] = []
            if weight_screen:
                low = {
                    uid
                    for uid, per_band in cutout_rows.items()
                    if _has_flag(per_band, "low_weight")
                }
                removed += [
                    {**r, "reason": "low cutout weight (D-021)"}
                    for r in ranked
                    if r["source_uid"] in low
                ]
                ranked = [r for r in ranked if r["source_uid"] not in low]
            if spike_screen:
                flagged = {
                    uid for uid, per_band in cutout_rows.items() if _has_flag(per_band, "spikes")
                }
                ranked = self._screen(label, summary, ranked, flagged, cutout_rows, removed)
            self._save_screened(label, summary, removed, ranked, cutout_rows, cut_k)
            if cut_table is not None:
                keep = {r["source_uid"] for r in ranked[:cut_k]}
                shown = cut_table[[_text(u) in keep for u in cut_table["source_uid"]]]
                self._contact_sheet(label, summary, shown, self._targets(label, ranked[:cut_k]))
        targets_all = self._targets(label, ranked[:n_top])
        self.save(targets_all, label, "targets")
        xmatch_rows = self._crossmatch(label, targets_all[:xm_k])

        for rec in ranked[:cand_k]:
            uid = rec["source_uid"]
            per_band = cutout_rows.get(uid, {})
            rec.update(
                sample_id=label,
                role=summary.role,
                flags={
                    band: {k: v for k, v in row.items() if k != "path"}
                    for band, row in per_band.items()
                },
                cutouts={band: row["path"] for band, row in per_band.items()},
                xmatch=xmatch_rows.get(uid),
            )
            self.candidates.append(rec)
        summary.topk = {
            **_topk_metrics(ranked[:cand_k], cutout_rows, xmatch_rows),
            "screened": summary.topk.get("screened", 0),
        }

    def _screen(
        self,
        sid: str,
        summary: SampleSummary,
        ranked: list[dict[str, Any]],
        flagged: set[str],
        cutout_rows: Mapping[str, Mapping[str, Mapping[str, Any]]],
        removed: list[dict[str, Any]],
    ) -> list[dict[str, Any]]:
        """``ranked`` without star-dominated spike-flagged sources (D-019, D-020).

        The removed sources are appended to ``removed`` with their reason.

        A spiky source is a star or a star-dominated blend when its colours are stellar (D-019)
        or, with ``host_ratio_max``, when no host light surrounds its peak (D-020; corrupted
        colours of star pairs or saturated stars do not matter then). A spiky nucleus inside a
        galaxy with non-stellar colours (e.g. an AGN) stays ranked and is noted. Without
        stellar-locus colours or a host test nothing is removed.
        """
        if not flagged:
            return ranked
        stellar = self.stellar_colour_uids.get(sid)
        colours = set(stellar) if stellar is not None else set()
        host_max = ((self.stages_cfg.get("cutouts") or {}).get("spike") or {}).get("host_ratio_max")
        if host_max is not None:
            hostless = {
                uid
                for uid in flagged
                if (
                    h := _finite_extreme(
                        _spike_bands(cutout_rows.get(uid) or {}), "host_ratio", min
                    )
                )
                is not None
                and h < float(host_max)
            }
            exempt = self._red_exempt(sid, hostless)
            if exempt:
                summary.notes.append(
                    "host test (D-020) exempts very red spiky source(s), not stars (D-025): "
                    + ", ".join("_".join(u.split("_")[-2:]) for u in sorted(exempt))
                )
            hostless -= exempt
            if stellar is None:
                summary.notes.append(
                    "spike screening: no stellar-locus colours for this sample; only the host "
                    "test (D-020) applies"
                )
            stellar = (stellar or set()) | hostless
        if stellar is None:
            summary.notes.append(
                "spike screening (D-019) skipped: no stellar-locus colours for this sample; "
                "spike-flagged sources stay ranked"
            )
            return ranked
        kept_flagged = [r for r in ranked if r["source_uid"] in flagged - stellar]
        if kept_flagged:
            summary.notes.append(
                "spike-flagged but non-stellar colours and host light, kept ranked (a galaxy "
                "with a bright nucleus such as an AGN, or a saturated star; inspect): "
                + ", ".join(
                    f"#{r['rank']} " + "_".join(r["source_uid"].split("_")[-2:])
                    for r in kept_flagged
                )
            )
        flagged = flagged & stellar
        removed += [
            {
                **r,
                "reason": (
                    "spikes, stellar colours (D-019)"
                    if r["source_uid"] in colours
                    else "spikes, no host light (D-020)"
                ),
            }
            for r in ranked
            if r["source_uid"] in flagged
        ]
        return [r for r in ranked if r["source_uid"] not in flagged]

    def _red_exempt(self, sid: str, uids: set[str]) -> set[str]:
        """Host-less spiky sources whose colour is too red for a star (D-025): kept ranked."""
        spike = (self.stages_cfg.get("cutouts") or {}).get("spike") or {}
        rule = spike.get("host_exempt_colour")
        if not rule or not uids or sid not in self.matched_sources:
            return set()
        sources, label = self.matched_sources[sid]
        blue, red, min_colour = str(rule[0]).lower(), str(rule[1]).lower(), float(rule[2])
        try:
            colour = column_as_float(sources, f"{blue}_{label}_abmag") - column_as_float(
                sources, f"{red}_{label}_abmag"
            )
        except KeyError:
            self.summaries[sid].notes.append(
                f"host-test exemption (D-025) off: no {blue}/{red} matched photometry here"
            )
            return set()
        uid_col = [_text(u) for u in sources["source_uid"]]
        with np.errstate(invalid="ignore"):
            return {
                u for u, c in zip(uid_col, colour, strict=True) if u in uids and c >= min_colour
            }

    def _save_screened(
        self,
        sid: str,
        summary: SampleSummary,
        removed: list[dict[str, Any]],
        kept: list[dict[str, Any]],
        cutout_rows: Mapping[str, Mapping[str, Mapping[str, Any]]],
        cut_k: int,
    ) -> None:
        """Save and report the sources screened out of the top k (D-019, D-020, D-021)."""
        if not removed:
            return
        removed = sorted(removed, key=lambda r: r["rank"] or 0)
        table = Table(
            {
                "source_uid": [r["source_uid"] for r in removed],
                "rank": np.array([r["rank"] or 0 for r in removed], dtype=int),
                "reason": [r["reason"] for r in removed],
                "spike_s6": np.array(
                    [
                        np.nan
                        if (
                            v := _finite_extreme(cutout_rows.get(r["source_uid"]) or {}, "spike_s6")
                        )
                        is None
                        else v
                        for r in removed
                    ],
                    float,
                ),
            }
        )
        table.meta.update(
            provenance=schema.Provenance.DERIVED.value,
            source=f"sources screened out of the top {cut_k} on cutout evidence; see 'reason'",
        )
        self.save(table, sid, "screened")
        summary.topk["screened"] = len(removed)
        summary.notes.append(
            f"screening removed {len(removed)} source(s) from the top {cut_k} (ranks keep their "
            "original numbers; see the screened table): "
            + ", ".join(
                f"#{r['rank']} " + "_".join(r["source_uid"].split("_")[-2:]) + f" ({r['reason']})"
                for r in removed
            )
        )
        if any(r["source_uid"] not in cutout_rows for r in kept[:cut_k]):
            summary.notes.append("screening: some backfilled sources have no cutout")

    def _classify(
        self, sid: str, sources: Table, to_rank: Table, label: str | None = None
    ) -> tuple[np.ndarray, int, int] | None:
        """``(star mask over to_rank, star_top_k, min_stars)`` (optional stage, D-012/D-015).

        None means one ranking: not configured, failed, or invalid parameters. With
        ``classify.stellar_locus`` and matched photometry ``label``, point-like sources with
        stellar colours join the catalogued stars (D-015); a locus failure keeps the catalogue
        classification and adds a note.
        """
        cfg = self.stages_cfg.get("classify")
        name = "classify.classify_sources"
        if cfg is None or cfg.get("enabled", True) is False:
            reason = "not configured" if cfg is None else "disabled in config"
            self.record(sid, name, "skipped", reason, required=False)
            return None
        params: dict[str, int] = {}

        def run_classify() -> Table:
            star_top_k = int(cfg.get("star_top_k", 10))
            min_stars = int(cfg.get("min_stars", 5))
            if star_top_k < 1 or min_stars < 2:
                raise ValueError("classify needs star_top_k >= 1 and min_stars >= 2")
            configured = {str(smp.get("id")) for smp in self.config.get("samples") or []}
            if f"{sid}-stars" in configured:
                raise ValueError(f"stratum id {sid}-stars collides with a configured sample")
            params.update(star_top_k=star_top_k, min_stars=min_stars)
            targets = Table(
                {"source_uid": sources["source_uid"], "ra": sources["ra"], "dec": sources["dec"]}
            )
            targets.meta.update(
                provenance=schema.Provenance.DERIVED.value, source=f"all merged sources of {sid}"
            )
            matches = crossmatch.query_matches(
                targets,
                radius_arcsec=float(cfg.get("radius_arcsec", classify.DEFAULT_QUERY_RADIUS_ARCSEC)),
                services=tuple(cfg.get("services") or classify.DEFAULT_SERVICES),
            )
            table = classify.classify_sources(
                sources, matches, **_stage_kwargs(cfg, ("gaia_radius_arcsec",))
            )
            if label:  # D-025: matched colours for the red exemption, independent of the locus
                self.matched_sources[sid] = (sources, label)
            veto_cfg = cfg.get("extended_veto")
            veto_on = bool(veto_cfg) and not (
                isinstance(veto_cfg, dict) and veto_cfg.get("enabled") is False
            )
            if veto_on and label:  # D-025: before the locus calibrates on these stars
                try:
                    table = classify.veto_extended_stars(table, sources, label, veto_cfg)
                except Exception as exc:  # noqa: BLE001  keys are validated in load_config
                    self.summaries[sid].notes.append(f"extended-star veto not applied: {exc}")
                else:
                    vetoed = table.meta["extended_veto"]["n_vetoed"]
                    if vetoed:
                        self.summaries[sid].notes.append(
                            f"extended-star veto (D-025): {vetoed} catalogue star(s) too "
                            "extended in the DJA detection image are ranked as galaxies"
                        )
            locus_cfg = cfg.get("stellar_locus")
            enabled = not (isinstance(locus_cfg, dict) and locus_cfg.get("enabled") is False)
            if locus_cfg is not None and locus_cfg is not False and enabled:
                if not label:
                    self.summaries[sid].notes.append(
                        "stellar locus not applied: no matched photometry for this sample"
                    )
                else:
                    try:  # an optional add-on: any failure keeps the catalogue classification
                        table = classify.apply_stellar_locus(table, sources, label, locus_cfg)
                        mask = classify.stellar_colour_mask(
                            sources, label, table.meta["stellar_locus"]
                        )
                        self.stellar_colour_uids[sid] = {
                            _text(u) for u, m in zip(sources["source_uid"], mask, strict=True) if m
                        }
                    except Exception as exc:  # noqa: BLE001
                        self.summaries[sid].notes.append(f"stellar locus not applied: {exc}")
            return table

        table = self.stage(
            sid,
            name,
            run_classify,
            columns=schema.CLASSIFY_COLUMNS,
            required=False,
            save_as="populations",
        )
        if table is None:
            self.summaries[sid].notes.append("star/galaxy separation failed: one ranking")
            return None
        star_uids = {
            _text(u)
            for u, p in zip(table["source_uid"], table["population"], strict=True)
            if _text(p) == "star"
        }
        mask = np.array([_text(u) in star_uids for u in to_rank["source_uid"]], dtype=bool)
        locus = table.meta.get("stellar_locus")
        if locus:
            locus_uids = {
                _text(u)
                for u, b in zip(table["source_uid"], table["star_basis"], strict=True)
                if _text(b) == "stellar_locus"
            }
            n_ranked = sum(_text(u) in locus_uids for u in to_rank["source_uid"])
            self.summaries[sid].notes.append(
                f"stellar locus (D-015/D-016) added {locus['n_added']} stars, "
                f"{n_ranked} of them past "
                f"the quality gate (r50_psf {locus['r50_psf_pix']} px from "
                f"{locus['n_calibration_stars']} catalogued stars)"
            )
        return mask, params["star_top_k"], params["min_stars"]

    def _filter_observations(
        self, sid: str, obs: Table, obs_ids: list[str], data_rights: str | None
    ) -> Table:
        notes = []
        keep = np.ones(len(obs), dtype=bool)
        if obs_ids:
            wanted = set(obs_ids)
            found = [_text(o) for o in obs["obs_id"]]
            keep &= np.array([o in wanted for o in found], dtype=bool)
            missing = sorted(wanted - set(found))
            if missing:
                notes.append(f"configured obs_ids not returned by MAST: {missing}")
        if data_rights:
            rights = np.array([_text(r).upper() for r in obs["dataRights"]])
            not_public = keep & (rights != str(data_rights).upper())
            if not_public.any():
                notes.append(
                    f"dropped {int(not_public.sum())} rows with dataRights != {data_rights}"
                )
            keep &= ~not_public
        obs = obs[keep]
        if notes:
            self.record(sid, "filter observations", "ok", "; ".join(notes), required=True)
        if len(obs) == 0:
            self.fail(sid, "filter observations", "no observations left after filtering")
        return obs

    def _ranked(self, sid: str, scores: Table, sources: Table) -> list[dict[str, Any]]:
        positions = {
            _text(uid): (_float(ra), _float(dec))
            for uid, ra, dec in zip(
                sources["source_uid"], sources["ra"], sources["dec"], strict=True
            )
        }
        method_cols = [c for c in scores.colnames if c.startswith("score_")]
        rows: list[dict[str, Any]] = []
        seen: set[str] = set()
        for row in scores:
            uid = _text(row["source_uid"])
            if uid not in positions:
                self.fail(sid, "select targets", f"scored source_uid {uid!r} not in sources")
            if uid in seen:
                self.fail(sid, "select targets", f"duplicate source_uid {uid!r} in scores")
            seen.add(uid)
            ra, dec = positions[uid]
            rank_value = to_python(row["rank"])
            top = to_python(row["top_features"]) if "top_features" in scores.colnames else None
            if isinstance(top, list):
                top = ", ".join(str(t) for t in top)
            rows.append(
                {
                    "source_uid": uid,
                    "ra": ra,
                    "dec": dec,
                    "score": _float(row["score"]),
                    "rank": None if rank_value is None else int(rank_value),
                    "top_features": None if top is None else str(top),
                    "method_scores": {c[len("score_") :]: _float(row[c]) for c in method_cols},
                }
            )
        rows.sort(key=lambda r: (r["rank"] is None, r["rank"] or 0, r["source_uid"]))
        return rows

    def _targets(self, sid: str, ranked: list[dict[str, Any]]) -> Table:
        targets = Table(
            {
                "source_uid": [r["source_uid"] for r in ranked],
                "ra": np.array([np.nan if r["ra"] is None else r["ra"] for r in ranked], float),
                "dec": np.array([np.nan if r["dec"] is None else r["dec"] for r in ranked], float),
                "rank": np.array([r["rank"] or 0 for r in ranked], dtype=int),
            }
        )
        targets["ra"].unit = targets["dec"].unit = "deg"
        targets.meta.update(
            provenance=schema.Provenance.DERIVED.value,
            source=f"top {len(ranked)} by rank.score_anomalies, run {self.run_id}, sample {sid}",
        )
        return schema.validate(targets, schema.TARGET_COLUMNS, name="targets")

    def _skip_reason(self, cfg: Mapping[str, Any] | None, targets: Table) -> str | None:
        if cfg is None:
            return "not configured"
        if cfg.get("enabled", True) is False:
            return "disabled in config"
        if len(targets) == 0:
            return "no targets"
        return None

    def _relative_to_run(self, path_text: str) -> str:
        """Cutout path relative to the run dir when it lies inside it (for report links)."""
        path = Path(path_text)
        absolute = (path if path.is_absolute() else Path.cwd() / path).resolve()
        try:
            return absolute.relative_to(self.run_dir.resolve()).as_posix()
        except ValueError:
            return path.as_posix()

    def _matched_photometry(
        self, sid: str, sample: Mapping[str, Any], sources: Table
    ) -> tuple[Table, str | None]:
        """Join external matched-aperture magnitudes (optional stage, D-013).

        Returns the table to build features from and the aperture label to use (None: the
        pipeline catalogs' default aperture).
        """
        cfg = sample.get("matched_photometry")
        name = "photometry.join_matched_photometry"
        if not cfg:  # a per-sample option, not a stage every sample runs
            return sources, None
        label = str(cfg["label"])

        def join() -> Table:
            fetch_kwargs = {"max_bytes": int(cfg["max_bytes"])} if "max_bytes" in cfg else {}
            path = photometry.fetch_catalog(str(cfg["url"]), str(cfg["sha256"]), **fetch_kwargs)
            catalog = photometry.load_dja_catalog(
                path,
                [b.lower() for b in self.summaries[sid].bands],
                aperture=int(cfg.get("aperture", 1)),
            )
            return photometry.join_matched_photometry(
                sources,
                catalog,
                label,
                radius_arcsec=float(
                    cfg.get("radius_arcsec", photometry.DEFAULT_MATCH_RADIUS_ARCSEC)
                ),
            )

        joined = self.stage(sid, name, join, columns=schema.SOURCE_COLUMNS, required=False)
        if joined is not None:
            self.save(photometry.joined_columns(joined, label), sid, "photometry")
        if joined is None:
            self.summaries[sid].notes.append(
                "matched-aperture photometry failed: colours from pipeline catalogs"
            )
            return sources, None
        info = joined.meta.get("matched_photometry", {})
        self.summaries[sid].notes.append(
            f"colours from matched-aperture photometry '{label}' "
            f"({info.get('n_matched')} of {len(joined)} sources matched one-to-one; D-013)"
        )
        return joined, label

    def _quality(
        self,
        sid: str,
        summary: SampleSummary,
        sources: Table,
        feats: Table,
        ref_band: str,
        band_obs: Mapping[str, str],
        confirm_column: str | None = None,
    ) -> np.ndarray | None:
        """Quality gate (optional, D-011/D-014). Return the rows of ``feats`` to rank, or None.

        None (rank everything) when the gate is not configured, fails, or leaves fewer than
        ``min_ranked`` sources, so a broken gate can never abort or empty a run.
        """
        cfg = self.stages_cfg.get("quality")
        name = "quality.assess_sources"
        if cfg is None or cfg.get("enabled", True) is False:
            reason = "not configured" if cfg is None else "disabled in config"
            self.record(sid, name, "skipped", reason, required=False)
            return None

        def gate() -> Table:
            params = _stage_kwargs(
                cfg,
                (
                    "min_rel_weight",
                    "min_edge_arcsec",
                    "max_artifact_ci",
                    "min_detection_snr",
                    "require_multiband",
                ),
            )
            if params.get("require_multiband") and confirm_column:
                params["confirm_column"] = confirm_column
            map_kwargs = {"grid_arcsec": float(cfg.get("grid_arcsec", 1.0))}
            if cfg.get("step") is not None:
                map_kwargs["step"] = int(cfg["step"])
            uri = l3_image_uri(self.config.get("cloud"), band_obs[ref_band])
            weights = cutouts.sample_weight_map(uri, **map_kwargs)
            table = quality.assess_sources(sources, weights, ref_band=ref_band, **params)
            ref = getattr(weights, "reference_weight", None)
            if ref is not None and np.isfinite(ref):
                self.weight_refs[(sid, uri)] = float(ref)
            return table

        table = self.stage(
            sid, name, gate, columns=schema.QUALITY_COLUMNS, required=False, save_as="quality"
        )
        if table is None:
            summary.notes.append("quality gate failed: all sources ranked (ungated)")
            return None
        ok_by_uid = {
            _text(u): bool(v)
            for u, v in zip(table["source_uid"], np.asarray(table["quality_ok"]), strict=True)
        }
        uids = [_text(u) for u in feats["source_uid"]]
        missing = sum(u not in ok_by_uid for u in uids)
        keep = np.array([ok_by_uid.get(u, False) for u in uids], dtype=bool)
        min_ranked = max(2, int(cfg.get("min_ranked", 20)))
        reasons_of = {
            _text(u): set(filter(None, _text(r).split(",")))
            for u, r in zip(table["source_uid"], table["quality_reason"], strict=True)
        }
        problem = None
        if missing:
            problem = f"gate table lacks {missing} of {len(uids)} feature rows"
        elif keep.sum() < min_ranked:
            # The D-014 tests (low_snr, single_band) can be too strict for a small or single-band
            # sample: fall back to the D-011 image tests before giving up on gating entirely.
            d011 = set(quality.REASONS_D011)
            keep_d011 = np.array([not (reasons_of.get(u, set()) & d011) for u in uids], dtype=bool)
            if keep_d011.sum() >= min_ranked:
                summary.notes.append(
                    f"detection-confirmation tests (D-014) left {int(keep.sum())} of {len(uids)} "
                    f"sources (min_ranked {min_ranked}): ranked with the D-011 image tests only"
                )
                keep = keep_d011
            else:
                n_ok = int(keep.sum())
                problem = f"only {n_ok} of {len(uids)} sources passed (min_ranked {min_ranked})"
        if problem:
            self.record(sid, "quality gate (sanity check)", "failed", problem, required=False)
            for key in [k for k in self.weight_refs if k[0] == sid]:
                del self.weight_refs[key]
            summary.notes.append(f"quality gate unusable ({problem}): all sources ranked (ungated)")
            return None
        summary.quality = quality.summarize(table)
        return keep

    def _cutouts(
        self,
        sid: str,
        summary: SampleSummary,
        targets: Table,
        band_obs: Mapping[str, str],
        parent: str | None = None,
        render: bool = True,
    ) -> tuple[dict[str, dict[str, dict[str, Any]]], Table | None]:
        """Run cutouts per band; return ``({source_uid: {band: {path, quality...}}}, table)``.

        ``render=False`` leaves the contact sheet to the caller (D-019 screening).
        """
        cfg = self.stages_cfg.get("cutouts")
        reason = self._skip_reason(cfg, targets)
        if reason:
            self.record(sid, "cutouts.make_cutouts", "skipped", reason, required=False)
            return {}, None
        wanted = [str(b).upper() for b in cfg.get("bands") or []]
        bands = [b for b in wanted if b in band_obs]
        if not bands:
            bands = [summary.ref_band]
            summary.notes.append(
                f"cutout bands {wanted} not in this sample; used ref_band {summary.ref_band}"
            )
        kwargs: dict[str, Any] = {"out_dir": self.run_dir / sid / "cutouts"}
        if "size_arcsec" in cfg:
            kwargs["size_arcsec"] = float(cfg["size_arcsec"])
        quality_cfg = self.stages_cfg.get("quality") or {}
        if "min_rel_weight" in quality_cfg:  # the cutout low_weight flag uses the D-011 threshold
            kwargs["low_weight_frac"] = float(quality_cfg["min_rel_weight"])
        spike = cfg.get("spike")
        if spike:  # D-018: flag bright stars and star+galaxy blends on the cutouts
            kwargs["spike_radii_arcsec"] = tuple(float(r) for r in spike["radii_arcsec"])
            kwargs["spike_threshold"] = float(spike["threshold"])  # "screen" is used by the caller
            if "search_arcsec" in spike:
                kwargs["spike_search_arcsec"] = float(spike["search_arcsec"])
            if "host_annulus_arcsec" in spike:
                kwargs["host_annulus_arcsec"] = tuple(
                    float(r) for r in spike["host_annulus_arcsec"]
                )
        tables = []
        for band in bands:
            name = f"cutouts.make_cutouts[{band}]"
            try:
                uri = l3_image_uri(self.config.get("cloud"), band_obs[band])
            except (ConfigError, ValueError) as exc:
                self.record(sid, name, "failed", str(exc), required=False)
                continue
            table = self.stage(
                sid,
                name,
                partial(
                    cutouts.make_cutouts,
                    uri,
                    targets,
                    **kwargs,
                    **(
                        {"weight_ref": self.weight_refs[(parent or sid, uri)]}
                        if (parent or sid, uri) in self.weight_refs
                        else {}
                    ),
                ),
                columns=schema.CUTOUT_COLUMNS,
                required=False,
            )
            if table is not None:
                tables.append(table)
                summary.cutout_bands.append(band)
        if not tables:
            return {}, None
        try:  # still part of the optional stage: a failure here must not abort the run
            combined = vstack(tables, metadata_conflicts="silent")
            combined.meta.update(
                provenance=tables[0].meta["provenance"],
                source="; ".join(str(t.meta["source"]) for t in tables),
            )
            self.save(combined, sid, "cutouts")
            result: dict[str, dict[str, dict[str, Any]]] = {}
            for row in combined:
                result.setdefault(_text(row["source_uid"]), {})[_text(row["band"]).upper()] = {
                    "path": self._relative_to_run(_text(row["path"])),
                    "quality_flag": to_python(row["quality_flag"]),
                    "on_edge": to_python(row["on_edge"]),
                    "frac_nan": to_python(row["frac_nan"]),
                }
                for col in ("spike_s6", "host_ratio"):
                    if col in combined.colnames:
                        result[_text(row["source_uid"])][_text(row["band"]).upper()][col] = _float(
                            row[col]
                        )
            spiky = sorted(
                {
                    _text(r["source_uid"])
                    for r in combined
                    if "spikes" in _text(r["quality_flag"]).split(",")
                }
            )
            if spiky:
                summary.notes.append(
                    f"{len(spiky)} of {len(targets)} cutout targets show diffraction spikes "
                    f"(bright star or star+galaxy blend, D-018): "
                    + ", ".join("_".join(u.split("_")[-2:]) for u in spiky)
                )
        except Exception as exc:
            self.record(sid, "cutouts (combine bands)", "failed", _describe(exc), required=False)
            return {}, None
        if render:
            self._contact_sheet(sid, summary, combined, targets)
        return result, combined

    def _contact_sheet(
        self, sid: str, summary: SampleSummary, combined: Table, targets: Table
    ) -> None:
        """Render the top-k cutouts as one PNG for visual inspection (optional)."""
        name = "viz.contact_sheet"
        start = time.perf_counter()
        try:
            ranks = {_text(r["source_uid"]): int(r["rank"]) for r in targets}
            out = viz.contact_sheet(
                combined,
                self.run_dir / sid / "contact_sheet.png",
                ranks=ranks,
                title=f"{sid}: top {len(targets)} by anomaly score (unvetted)",
            )
        except Exception as exc:
            self.record(sid, name, "failed", _describe(exc), required=False)
            return
        summary.contact_sheet = self._relative_to_run(str(out))
        self.records.append(
            StageRecord(
                sid,
                name,
                False,
                "ok",
                n_rows=len(combined),
                seconds=time.perf_counter() - start,
                output=summary.contact_sheet,
            )
        )

    def _crossmatch(self, sid: str, targets: Table) -> dict[str, dict[str, Any]]:
        cfg = self.stages_cfg.get("crossmatch")
        reason = self._skip_reason(cfg, targets)
        if reason:
            self.record(sid, "crossmatch.crossmatch", "skipped", reason, required=False)
            return {}
        kwargs = _stage_kwargs(cfg, ("radius_arcsec", "services"))
        table = self.stage(
            sid,
            "crossmatch.crossmatch",
            partial(crossmatch.crossmatch, targets, **kwargs),
            columns=schema.XMATCH_COLUMNS,
            required=False,
            save_as="crossmatch",
        )
        if table is None:
            return {}
        meta = {
            k: to_python(v) for k, v in table.meta.items() if k in ("services", "radius_arcsec")
        }
        out = {}
        for row in table:
            cols = [c for c in schema.XMATCH_COLUMNS if c != "source_uid"]
            cols += [c for c in ("is_lens_related", "lens_types") if c in table.colnames]
            values = {c: to_python(row[c]) for c in cols}
            values.update(meta)
            out[_text(row["source_uid"])] = values
        return out


# -- report ------------------------------------------------------------------------------------


def _quality_line(counts: Mapping[str, int]) -> str:
    if not counts:
        return "not applied (all sources ranked)"
    excluded = {k: v for k, v in counts.items() if k not in ("n_sources", "n_ok")}
    detail = ", ".join(f"{k} {v}" for k, v in sorted(excluded.items())) or "none"
    return f"{counts['n_ok']} of {counts['n_sources']} passed; flags: {detail}"


def _cell(value: Any) -> str:
    text = "" if value is None else str(value)
    return text.replace("|", "\\|").replace("\n", " ").strip()


def _num(value: Any, fmt: str) -> str:
    return "" if value is None else format(value, fmt)


def _quality_cell(flags: Mapping[str, Any] | None) -> str:
    if not flags:
        return "not run"
    parts = []
    for band, q in sorted(flags.items()):
        text = f"{band}: {q.get('quality_flag')}"
        if q.get("on_edge"):
            text += ", edge"
        frac_nan = _float(q.get("frac_nan"))
        if frac_nan:
            text += f", NaN {frac_nan:.0%}"
        parts.append(text)
    return "; ".join(parts)


def _xmatch_cell(xm: Mapping[str, Any] | None) -> str:
    if not xm:
        return "not run"
    n_matches = xm.get("n_matches")
    if n_matches is None:
        return "result missing"
    if not n_matches:
        radius = xm.get("radius_arcsec")
        return "no counterpart" + (f" within {radius} arcsec" if radius is not None else "")
    sep = _float(xm.get("best_match_sep_arcsec"))
    text = f"{xm.get('best_match_service')}: {xm.get('best_match_id')}"
    detail = [str(xm.get("best_match_type") or "?")]
    if sep is not None:
        detail.append(f'{sep:.2f}"')
    if xm.get("is_star"):
        detail.append("star")
    if xm.get("is_known_object"):
        detail.append("known")
    return f"{text} ({', '.join(detail)}; {n_matches} matches)"


def render_report(
    context: Mapping[str, Any],
    config: Mapping[str, Any],
    runner: _Runner,
    *,
    status: str,
    error: str | None = None,
    db_path: Path | None = None,
) -> str:
    """Markdown report for one run (see module docstring)."""
    lines: list[str] = [f"# Anomaly-ranking run `{context['run_id']}`", ""]
    lines += [f"> **Read this first.** {DISCLAIMER}", ""]
    lines.append(f"**Status:** {status}" + (f" -- {error}" if error else ""))
    if status != "completed":
        lines += ["", "No candidates were stored: the run did not complete."]
    lines.append("")

    git = context.get("git") or {}
    cfg = context.get("config") or {}
    if git.get("commit"):
        dirty = {True: "yes", False: "no"}.get(git.get("dirty"), "unknown")
        code = (
            f"git `{git['commit'][:12]}` on `{git.get('branch')}` "
            f"(uncommitted code changes: {dirty})"
        )
    else:
        code = "not a git checkout (code version unknown beyond package version)"
    packages = ", ".join(f"{k} {v}" for k, v in (context.get("packages") or {}).items() if v)
    py = context.get("python") or {}
    plat = (context.get("platform") or {}).get("platform", "")
    lines += [
        "## Run context",
        "",
        "| Item | Value |",
        "|---|---|",
        f"| Run id | `{context['run_id']}` |",
        f"| Created (UTC) | {context.get('created_utc')} |",
        f"| Config | `{_cell(config.get('name'))}` from `{_cell(cfg.get('path'))}` "
        f"(sha256 `{(cfg.get('sha256') or '')[:16]}`; copy: `config.yaml`) |",
        f"| Code | {_cell(code)} |",
        f"| jwst-anomaly | {context.get('jwst_anomaly_version')} |",
        f"| Python | {py.get('version')} ({py.get('implementation')}) on {_cell(plat)} |",
        f"| Packages | {_cell(packages)} |",
        f"| Settings | {_cell(json.dumps(runner.settings()))} |",
        "",
    ]

    lines += [
        "## Stages",
        "",
        "| Sample | Stage | Required | Status | Rows | Time (s) | Notes / output |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in runner.records:
        note = r.message or (f"`{r.output}`" if r.output else "")
        n_rows = "" if r.n_rows is None else r.n_rows
        lines.append(
            f"| {_cell(r.sample)} | `{_cell(r.stage)}` | {'yes' if r.required else 'no'} | "
            f"{r.status} | {n_rows} | {r.seconds:.2f} | {_cell(note)} |"
        )
    lines.append("")

    lines += ["## Samples", ""]
    for s in runner.summaries.values():
        lines += [f"### {s.id} ({s.role})", ""]
        if s.description:
            lines += [s.description, ""]
        lines += [
            f"- Observations (observed): {s.n_observations}; catalog bands: "
            f"{', '.join(s.bands) or 'none'}; ref band {s.ref_band}",
            f"- Merged sources (derived): {s.n_sources}",
            f"- Features (derived): {len(s.feature_spec)}",
            f"- Quality gate (derived, D-011): {s.quality_text or _quality_line(s.quality)}",
            f"- Scored sources (model_prediction): {s.n_scored}; methods: "
            f"{', '.join(s.methods) or 'unknown'}",
            f"- Cutout bands: {', '.join(s.cutout_bands) or 'none'}",
        ]
        if s.topk:
            lines.append(f"- Top-k composition (derived, tracked metric): {_topk_line(s.topk)}")
        if s.tables:
            lines.append(f"- Tables: {', '.join(f'`{t}`' for t in s.tables)}")
        if s.contact_sheet:
            lines.append(f"- Contact sheet: [{s.contact_sheet}]({s.contact_sheet})")
        lines += [f"- Note: {_cell(n)}" for n in s.notes]
        lines.append("")
        cands = [c for c in runner.candidates if c["sample_id"] == s.id]
        if cands:
            methods = s.methods
            lines += [
                f"#### Top {len(cands)} candidates",
                "",
                "Score columns: `model_prediction`. Image quality: `derived`. "
                "Cross-match: `observed` (external catalogs). Status of all rows: `new`.",
                "",
                "| Rank | source_uid | RA (deg) | Dec (deg) | Score | "
                + "".join(f"{m} | " for m in methods)
                + "Top features | Image quality | Cross-match | Cutouts |",
                "|" + "---|" * (9 + len(methods)),
            ]
            for c in cands:
                links = ", ".join(
                    f"[{band}]({path})" for band, path in sorted((c.get("cutouts") or {}).items())
                )
                method_cells = "".join(
                    f"{_num((c.get('method_scores') or {}).get(m), '.3g')} | " for m in methods
                )
                lines.append(
                    f"| {c.get('rank')} | `{_cell(c['source_uid'])}` | {_num(c.get('ra'), '.6f')} "
                    f"| {_num(c.get('dec'), '.6f')} | {_num(c.get('score'), '.3g')} | "
                    f"{method_cells}{_cell(c.get('top_features'))} | "
                    f"{_cell(_quality_cell(c.get('flags')))} | "
                    f"{_cell(_xmatch_cell(c.get('xmatch')))} | {links or 'none'} |"
                )
            lines.append("")
        if s.feature_spec:
            lines += ["<details><summary>Feature definitions (derived)</summary>", ""]
            lines += [f"- `{k}`: {_cell(v)}" for k, v in s.feature_spec.items()]
            lines += ["", "</details>", ""]

    lines += ["## Limitations", ""]
    lines += [f"- {item}" for item in LIMITATIONS]
    db = f' --db "{Path(db_path).as_posix()}"' if db_path is not None else ""
    lines += [
        "",
        "## Next steps",
        "",
        f"- `jwst-anomaly candidates list{db} --run {context['run_id']}`",
        "- Vet a candidate (the `/vet-candidate` skill walks through the tests): "
        f"`jwst-anomaly candidates set-status{db} <source_uid> triaged --run {context['run_id']} "
        "--note ...`; record each test with `candidates add-vetting <source_uid> --test ... "
        "--outcome pass|fail|inconclusive --evidence ... --provenance ...`; conclude with "
        "`candidates set-status <source_uid> artifact|known_object|explained|unexplained "
        "--note ...`. Add `--sample <id>` when a source_uid occurs in several samples.",
        "",
    ]
    return "\n".join(lines)


# -- entry point -------------------------------------------------------------------------------


def run(
    config_path: str | Path,
    *,
    outputs_dir: str | Path | None = None,
    db_path: str | Path | None = None,
    samples: Sequence[str] | None = None,
) -> str:
    """Run the pipeline described by a YAML config and return the run id.

    ``outputs_dir`` defaults to ``paths.outputs_dir()``; ``db_path`` to
    ``<outputs_dir>/candidates.sqlite``; ``samples`` restricts the run to those sample ids.
    Raises :class:`ConfigError` before anything is recorded if the config is invalid, and
    :class:`PipelineError` (after recording a ``failed`` run and its report) if a required
    stage fails. Any other exception, including ``KeyboardInterrupt``, is also recorded as a
    ``failed`` run before it propagates.
    """
    config_path = Path(config_path)
    config = load_config(config_path)
    selected = list(config["samples"])
    if samples:
        wanted = set(samples)
        unknown = sorted(wanted - {str(s["id"]) for s in selected})
        if unknown:
            raise ConfigError(f"unknown sample ids {unknown}")
        selected = [s for s in selected if str(s["id"]) in wanted]

    out_root = (Path(outputs_dir) if outputs_dir is not None else paths.outputs_dir()).resolve()
    db = Path(db_path) if db_path is not None else default_db_path(out_root)
    context = provenance.capture_run_context(config_path)
    run_id = context["run_id"]
    run_dir = out_root / "runs" / run_id
    runner = _Runner(
        config, run_id, run_dir, manifest_path=paths.manifests_dir() / f"{config_path.stem}.ecsv"
    )
    run_dir.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(config_path, run_dir / "config.yaml")
    (run_dir / "run_context.json").write_text(json.dumps(context, indent=2), encoding="utf-8")
    report_path = run_dir / "report.md"
    log.info("run %s: config %s -> %s", run_id, config_path, run_dir)

    with CandidateStore(db) as store:
        store.add_run(run_id, context, config_name=config["name"])

        def finish(status: str, error: str | None) -> None:
            if status != "completed":
                runner.candidates.clear()  # nothing is stored for a failed run
            summary = {
                "status": status,
                "error": error,
                "settings": runner.settings(),
                "samples_selected": [str(s["id"]) for s in selected],
                "stages": [asdict(r) for r in runner.records],
                "samples": [asdict(s) for s in runner.summaries.values()],
                "n_candidates": len(runner.candidates),
                "report": report_path.name,
            }
            store.finish_run(run_id, status, summary=summary, error=error)  # the record of truth
            try:
                (run_dir / "run_record.json").write_text(
                    json.dumps(to_python(summary), indent=2), encoding="utf-8"
                )
                report_path.write_text(
                    render_report(context, config, runner, status=status, error=error, db_path=db),
                    encoding="utf-8",
                )
            except OSError as exc:
                log.error("run %s: could not write run record/report: %s", run_id, exc)

        def finish_failed(error: str) -> None:
            try:
                finish("failed", error)
            except Exception:  # never mask the original failure
                log.exception("run %s: could not record the failure", run_id)

        try:
            for sample in selected:
                runner.run_sample(sample)
            store.add_candidates(run_id, runner.candidates)
        except _StageFailure as exc:
            finish_failed(str(exc))
            raise PipelineError(str(exc), run_id=run_id, report_path=report_path) from exc
        except BaseException as exc:  # incl. KeyboardInterrupt: never leave a run 'running'
            interrupted = isinstance(exc, KeyboardInterrupt)
            finish_failed("interrupted" if interrupted else f"internal error: {_describe(exc)}")
            raise
        finish("completed", None)
    log.info(
        "run %s completed: %d candidates; report %s", run_id, len(runner.candidates), report_path
    )
    return run_id
