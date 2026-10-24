# MARS MACS1149 map vs the Schuldt et al. 2024 new multiple images: gate H fails (D-034 addendum at)

Plan and image list frozen in 28d41c36; script and tests in 9a515ebf; both before the run. Numbers from `summary.json`
(per-image rows in `heldout_MACSJ1149.ecsv`) unless marked post hoc (`posthoc.json`, `posthoc_canucs.ecsv`).
Labels: Schuldt+2024 and MARS positions and listed z are observed (as the authors list them); residuals are model
predictions; thresholds are ASSUMPTIONs (`../mars_hff_heldout/PLAN.md`).

## Frozen result
| Field | Listed images | Frame offset O | Seen / held-out | Types E / N / M | Outside map | Gate H rms | Within 3″ | H | Training K rms |
|---|---|---|---|---|---|---|---|---|---|
| MACS1149 | 18 (Schuldt+2024 Table 2) | not applied (1 pair < 5) | 0 / 18 | 0 / 17 / 0 | 7 (301, 401) | 2.51″ (10 tested) | 70 % | **fail** | 0.742″ |

- **The HST-era MARS map does not reproduce the post-2023 MACS1149 images.** None of the 18 is in the MARS
  catalogue (no row within 0.5″), so every family is type N (source = mean of all its images, held-out ones
  included). SXa has no partner and is untested, as expected.
- Families 301 and 401 (7 images, ~50″ north of the BCG) fall outside the MARS map and are dropped and counted.
- Of the 10 tested images (families 201–204), 3 are beyond 3″: 202.1 (3.72″), 203.1 (3.47″), 203.2 (3.36″); the
  rest 1.30–2.41″. Not flags (frozen: no flag reading).

## Post hoc: the same families under the CANUCS v1 maps (training-set check, not held-out)
`scripts/mars_m1149_heldout_posthoc.py`, same rule 3 and grid step. The CANUCS v1 model was fitted to a gold list
that contains Schuldt+2024's systems (CANUCS readme), so this is not a held-out test. All 17 family images lie inside
its maps: rms 0.849″, 17 / 17 within 3″ (max 1.77″, 401.4); on the 10 images MARS tested, 0.716″, 10 / 10.

## Reading (narrow)
- The six Schuldt+2024 families are self-consistent under a model that fitted them (post hoc), which makes a wrong
  family assignment or z an unlikely cause of the MARS residuals. The miss is best read as MARS-specific: model
  incompleteness away from its HST-era constraints, coarse map pixels near critical curves, and a small map.
- Of three fields with a post-MARS image list, MARS passes gate H on one (MACS0416) and fails on two (Abell 370,
  MACS1149). MARS is therefore not a validated second model on MACS1149; no two-model extra-image stage there.
- Type N is a weaker test than type E (the held-out images anchor their own source), so a pass would have been
  modest evidence; the fail is not weakened by that.
- No image pixels were read; no contact sheet (catalogue positions only).

## Not done / limits
- 10 tested images in 4 families; 7 outside the map. One model's coarse maps; bilinear, mean source.
- The CANUCS v1 MACS1149 Lenstool input list is not public (readme: with v2), so the gold list is not stored here.
