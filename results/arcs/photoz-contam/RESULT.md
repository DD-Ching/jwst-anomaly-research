# `arcs`: photo-z contamination of the 3 single-component pooled anti rows (exploratory)

Plan: `PLAN.md` (committed before the run). Script: `scripts/arcs_photoz_contam.py`; rule functions unit-tested
before the run (`tests/test_arcs_photoz_contam.py`). Every number below is in `summary.json`. Inputs: RELICS HST
HLSP catalogues (DOI 10.17909/T9SP45), NED `objdir` (TAP), DJA v4.6 NIRSpec table; all read in the session
scratchpad and deleted, none stored.

## Frozen outcome: **mixed / undetermined** (derived)
All three rows are `undetermined`: no spectroscopic redshift exists, and the red-sequence test fails its own gate
in every cluster.

| Row | z_cl | Spec-z (NED 1.0″ / DJA 0.5″) | Stored photo-z (z160–z840) | RELICS BPZ zb (zbmin–zbmax, odds) | Second photo-z ≤ member cut | Label |
|---|---|---|---|---|---|---|
| WHL0137 2194 | 0.566 | none / none | DJA eazy 1.17 (1.14–1.27) | 1.079 (0.154–1.732, 0.165) | yes (cut 0.613) | undetermined |
| RXC J0600 581 | 0.43 | none / none | DJA eazy 0.70 (0.675–0.77) | 0.537 (0.457–0.665, 0.805) | yes (cut 0.473) | undetermined |
| MACS0308 2157 | 0.356 | none (one NED galaxy at 0.93″ without z) / none | RELICS BPZ 0.756 (0.671–0.853) | same catalogue | n/a | undetermined |

RELICS matches: 0.32″, 0.23″, 0.06″ (MACS0308 after the stored HST-to-JWST shift).

## Red-sequence gate (validation; derived)
F606W − F814W vs F814W, fitted to RELICS galaxies with F814W < 23 and |zb − z_cl| ≤ 0.05 (1 + z_cl):
| Cluster | n fit | σ_RS | Member recall | False-member rate f_bg (n) | Gate (recall ≥ 0.8, f_bg ≤ 0.3) |
|---|---|---|---|---|---|
| WHL0137 | 126 | 0.243 | 0.944 | 0.72 (50) | failed |
| RXC J0600 | 158 | 0.166 | 0.911 | 0.839 (62) | failed |
| MACS0308 | 102 | 0.093 | 0.882 | 0.618 (89) | failed |

Most BPZ-background galaxies brighter than F814W 23 sit within 2 σ_RS of the red sequence, so this colour cannot
tell members from background at these redshifts. The test is `undetermined` everywhere.

## Post hoc (not in the plan; never used for a label)
- Red-sequence offsets despite the failed gate: 2194 −1.86 σ, 581 +2.48 σ (redder than the sequence), 2157 +0.45 σ.
  With f_bg 0.62–0.84 these offsets carry little information either way.
- The two photo-z of 581 disagree: DJA eazy puts its 68 % interval at 0.675–0.77, and RELICS BPZ puts it at 0.537
  (0.457–0.665, whose lower end is below the member cut 0.473). Both best values lie behind the cluster, but HST
  photometry alone does not exclude a member. 2194's BPZ is unconstrained (odds 0.165, F814W 26.4).
- Contact sheet (`../pooled-census/census_anti.png`, looked at): 581 is a thin elongated disc with a bright core,
  consistent with a near edge-on disc whose large ellipticity (0.6875 stored, observed) dominates any shear. That is
  a reading by eye; it is not tested.

## Reading (narrowest statement)
- No public spectroscopic redshift exists for any of the 3 rows. Photometric tests do not decide membership: the
  F606W − F814W red sequence fails validation, and only 581 has a second photo-z that reaches the member cut.
- The open part of the pooled anti excess (D-024 addendum 9, post hoc P 0.024) is neither explained nor
  strengthened by this test. No candidate.

## Deviations
- The DJA table was streamed from the same URL as #781; its checksum was not compared this run. No row lies within
  0.5″ of any target.
- `posthoc_rs_offset_sigma` was added to the script after the first run (it reports the offset for gate-failed
  clusters); labels and the outcome are unchanged.

## Failed approach
- An HST F606W − F814W red sequence does not separate z_cl < z ≲ 1 background galaxies from members at z_cl
  0.36–0.57 (f_bg 0.62–0.84). Do not reuse it as a membership veto. A redder-baseline colour (e.g. F435W − F814W or
  JWST F090W − F200W) would need its own gate.

## Next
Membership needs a spectrum (none public: NED, DJA v4.6). Without one, the planned next test is a second-band shape
(F150W / F444W) of the 3 rows and the local shear error of the deciding models, not more redshift work.
