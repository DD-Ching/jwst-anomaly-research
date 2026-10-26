# Cosmic Gems arc (SPT-CL J0615−5746) vs the HST-era RELICS Lenstool v1 and GLAFIC v2 maps (2026-10-26)

Frozen plan `PLAN.md` (pushed in f52f24e4 before any map was evaluated at the arc; see Deviation). Numbers come from
`summary.json`, written by `scripts/spt0615_gems.py`. Clump positions: Adamo et al. 2024 Extended Data Table 3;
counter-image candidate and JWST-era predictions: Bradley et al. 2025 (observed / model_prediction as labelled in
the script). Crossings, magnifications and predicted images are model_prediction; classes derived; thresholds
ASSUMPTION.

## Frozen result (both maps usable: gates P and F from #537)
| z | Lenstool v1: crossing (μ A.1 / A.2) | nearest tangential curve | counter-image d to candidate | GLAFIC v2: crossing (μ A.1 / A.2) | nearest tangential curve | counter d |
|---|---|---|---|---|---|---|
| 10.0 | `miss` (−5.68 / −5.14) | 7.89″ | 2.526″ `offset` | `miss` (−4.03 / −3.74) | none within 15″ | 18.616″ `none` |
| 10.2 | `miss` (−5.61 / −5.09) | 7.95″ | 2.335″ `offset` | `miss` (−4.01 / −3.72) | none within 15″ | 18.898″ `none` |
| 10.4 | `miss` (−5.54 / −5.03) | 7.98″ | 2.141″ `consistent` | `miss` (−3.99 / −3.70) | none within 15″ | 19.169″ `none` |

- Rule 1: neither map has a tangential critical-curve crossing on the 1.447″ A.1 → A.2 segment at any z; rule 2
  agrees (A.1 and A.2 have the same, negative parity in both maps).
- Rule 4: **joint flag** (every usable map `miss` at every z). Pre-declared reading: an HST-era model limitation,
  not a lensing anomaly (the JWST-era Lenstool-A and glafic models put the curve through the arc even without the
  arc as a constraint, Bradley et al.).
- Lenstool v1's counter-image lies 2.141–2.526″ from the candidate, about as close as the JWST-era models'
  predictions (1.8″ and 2.2″ per Bradley et al.'s text; not stored); its distance to the JWST Lenstool-A prediction is 3.635–4.164″. GLAFIC v2's nearest
  counter-image is 18.616–19.169″ away (17.127–17.896″ from the Lenstool-A prediction).

## Cross-check against a known case
Salmon et al. 2018 (arXiv:1801.03103, abstract) give μ ~ 4–7 for SPT0615-JD from the RELICS models and predicted
counter-images. The maps here give |μ| 3.70–5.68 at the clumps (same order; the HLSP versions and the exact
evaluation points of that paper are not known), so the tool reproduces the published HST-era picture of a single,
moderately magnified image.

## Reading
- The HST-era maps describe the arc as one negative-parity image ~8″ (Lenstool v1) or > 15″ (GLAFIC v2) inside
  the z ≈ 10.2 tangential curve; the JWST imaging shows a mirror-symmetric pair across the curve (Adamo et al.). The
  arc's critical-curve location is therefore a residual of both HST-era maps; the ordinary explanation (stated, not
  tested here) is that no HST-era constraint lay near the arc. Photo-z within 10.0–10.4 does not change the class.
- Lenstool v1's third-image prediction already falls within ~2.3″ of the JWST candidate (`offset` at 10.0 and
  10.2, `consistent` at 10.4), in line with Bradley et al.'s statement that earlier models predicted a counter-image
  near it.
- Consequence for the repository: near the arc, neither SPT0615 RELICS map is reliable for magnification or
  parity (addendum bl's "usable for families 1–3" stays as is; it does not extend to the arc region).

## Limits
- One arc, two maps; the photo-z (no z_spec) and a single lens plane are ASSUMPTIONs. The GLAFIC v2 constraint
  list is undocumented (ASSUMPTION: the arc not used).
- Nearest-curve distances come from a 0.1″ grid in a ±15″ box; "none within 15″" is all that is stored for GLAFIC.
- Plot of both maps' z 10.2 critical curves with the clumps and the candidate looked at (not committed).

## Deviation
- The rule functions and their unit tests were written and passing before the run, but their commit failed on a
  lint error (line length) and was made after the run; only two docstring / string line breaks changed, and the
  run was repeated after the fix with the same summary.
- The plan says "the #537 frame offsets are applied": the configured frame offset is 0 for both maps (as in #537,
  where frame F passed; largest median offset −0.2236″ in `spt0615_images/summary.json`); nothing was applied.
