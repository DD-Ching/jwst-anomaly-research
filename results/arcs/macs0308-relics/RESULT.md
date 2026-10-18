# `arcs` on MACS J0308.9+2645: result (2026-10-18)

Plan: `PLAN.md` (commit 84c9f34, 19:12:28Z, before any `arcs` count; the "19:20Z" in its header was an estimate, the
commit time is the record). Script `scripts/arcs_macs0308.py` (reuses the SPT0615 BPZ conversion and astrometric
shift, unit-tested in `tests/test_arcs_spt0615.py`; wiring in `tests/test_arcs_macs0308.py`; driver
`arcs_rxj2129.decide` and outcome rule unchanged). Sheets `scripts/arcs_macs0308_sheet.py`. Labels: shapes
`observed`; photo-z `derived` (RELICS BPZ); predicted axes `model_prediction`; classes, counts and P values `derived`;
thresholds ASSUMPTIONs (D-024 defaults). Exploratory: RELICS BPZ, 95 % interval as z160–z840.

## Inputs (pinned)
- F200W catalogue `jw05293-o005_t001_nircam_clear-f200w_cat.ecsv` (MAST, PUBLIC, 1,578 s), sha256 in `summary.json`.
- RELICS BPZ `hlsp_relics_hst_acs-wfc3ir_macs0308p26_multi_v1_cat.txt` → 1,866 galaxy rows (stel < 0.8). Astrometry
  (`photoz_offset.json`): median JWST − BPZ offset +0.342″ (RA) and −0.104″ (Dec) over 868 pairs, above the 0.1″
  rule, so it was removed before the 0.3″ match (as planned).
- RELICS maps `macs0308-relics-glafic` (v2, primary) and `macs0308-relics-zitrin` (LTM-Gauss v1, second); sha256 in
  `lensmodel.MACS0308_RELICS_*`. Both pass `validate`: κ median relative difference 2.2 × 10⁻⁵ (GLAFIC) and
  1.1 × 10⁻⁵ (Zitrin) (`<model>/validate.json`).

## Screen (`summary.json`, `<model>/arcs.json`)
| Model | strong-shear background | aligned | anti | expected anti (shape noise) | P(≥ n_anti) | control (≥ 70 %) |
|---|---|---|---|---|---|---|
| GLAFIC v2 (primary) | 21 | 17 (81 %) | 1 | 0.309 | 0.266 | passes |
| Zitrin LTM v1 (second) | 18 | 15 (83 %) | 1 | 0.390 | 0.323 | passes |

**Frozen outcome: "does not replicate"** (primary passes its control; P 0.27 ≥ 0.05). The second model agrees.

## Strong-shear background rows, looked at (`background_glafic.png`, `background_zitrin.png`; 4″ F200W / F410M)
Deviation (secondary check only): program 5293 o005 has no F444W, so F410M is the second band.
- **2157**, the only anti in both models (m 20.6, S/N 505, e 0.53, z_phot 0.76 [0.67–0.85], r 19.5″): a bright,
  smooth elongated galaxy with a fainter neighbour, same axis in F200W and F410M. The predicted axis flips between
  z = 1 and z = 2 here (offset 82° at z ≤ 1, 8° at z = 2 and 4 in GLAFIC; 79° / 11° in Zitrin; g at z = 1: 2.1 and
  9.1), so the drawn z = 2 axis matches the observed one and the class rests entirely on the BPZ redshift being
  < 1. Zitrin puts it at g_min 1.19 (critical-curve flip region). Ordinary explanations: a critical curve for z ≈ 1
  sources passing near it, and the photo-z (a BPZ error to z ≳ 1 would make it aligned). Not a candidate.
- The `mixed` rows (290, 1140, 2043 in GLAFIC; 1140, 2288 in Zitrin) were looked at; 2043 sits on a bright, flat
  background (the stamp scale shows no structure) and 2288 is a stretched arc at g_min 5.7 / 33 (flip region).
- No row is sent to `/vet-candidate`.

## Limits
- One shallow F200W epoch (1,578 s); HST-only BPZ photometry; both models are HST-era (Acebron et al. 2018).
- 21 / 18 strong-shear background rows (expected anti 0.31 / 0.39): the test has power only against an excess of
  several anti rows, like S1063's 4 vs 0.52.
