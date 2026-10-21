# MACS0717: D-041 cross-model readings with the pinned `macs0717-keeton` and `macs0717-sharon` (frozen plan)

Frozen 2026-10-21 before any cross-model run. Script: `scripts/m0717_crossmodel.py`.

**Why.** #420 pinned Keeton v4 (both gates pass) and Sharon v4cor (G2 fails as written; in-footprint use only) as
`MODELS` entries. The D-041 verdicts in `docs/fields/macs0717.md` Step 6 ("CATS-only extra image", "all three models
place an image here") came from scratch code (2026-10-08). This run re-reads **every** non-`observed` row of the CATS
v4.1 `images` table with the pinned entries, so the MACS0717 forced-photometry verdicts rest on pinned inputs. Same
method as `results/lens_consistency/macs1149_crossmodel/PLAN.md`, extended to two second models.

**Hypothesis.** No CATS extra image with a forced flag (`absent`, `confused`, `inconsistent_reference`) is predicted
by a pinned second model and left without an ordinary explanation (D-041: 51 flags, 0 surviving).
**Ordinary explanations** for any disagreement: CATS-only image (model-dependent), a catalogued image misplaced by
CATS (CATS rms 3.2″), model copies (D-041 copy rule), near-critical μ (D-036), a fitted source redshift (41 of 51
systems take z from `params.txt` `z_m_limit`, D-035). **Would change the plan:** a row with final outcome
`to_vet_candidate` (below), which goes to `/vet-candidate`.

**Inputs.**
- CATS rows: `lens_consistency.py --model macs0717-cats images` re-run as in Step 4 of the field note (F200W
  `_cat.ecsv`, `--step 0.25`, forced F277W `_i2d` by S3 byte range); tracked as
  `results/lens_consistency/macs0717-cats/images_predicted.ecsv` (`model_prediction`; forced classes `derived`).
  **Reproduction check (computable):** its forced-class counts must equal Step 4's (recovered 25, confused 38,
  absent 6, undetectable 13, ambiguous 0, no_reference 54, inconsistent_reference 5, off_image 3) and its image
  classes Step 4's (observed 40, demagnified 2, candidate 68, faint 42, no_flux_ref 34, outside 13). If not, the
  difference is reported as a deviation and the new table is used.
- CATS `arcs.txt` (`observed`); `macs0717-keeton` and `macs0717-sharon` x/y maps (sha256-pinned). No photometry
  beyond the CATS run is measured.

**Procedure (per second model M, per system with a non-`observed` row).**
1. Footprint: every catalogued image of the system and the CATS row inside M's maps (finite deflection); else the
   row is not tested by M.
2. M's images of the system's magnification-weighted least-squares source at the CATS `z_sys`
   (`m1149_refs.second_model_images`, grid step 0.25″), each assigned to the nearest catalogued image within 1.5″
   (`m1149_refs.split_second`); the rest are M's extras.
3. D-041 reading per CATS row (`m1149_refs.d041_verdict`, loose radius 3.0″, ASSUMPTION): `both models`,
   `catalogued image misplaced by CATS`, or `model-dependent (CATS only)`.

**Known-case gate (per model, computable).** Over the systems with every catalogued image inside M's maps, the
share of catalogued images that M reproduces (an M image within 1.5″) at the CATS `z_sys` must be ≥ 0.8
(ASSUMPTION). A model that fails is not used ("tool not validated" for that model). Because most CATS redshifts are
model-fitted, a pre-declared secondary gate counts only systems with a redshift in `arcs.txt`; it is reported, not
used.

**Frozen rules (functions `combine`, `outcome`, `final`).**
- `combine`: over the models that pass the gate and test the row: `both models` if any reads `both models`; else
  `catalogued image misplaced by CATS` if any reads it; else `model-dependent (CATS only)`; empty if no model tests
  the row.
- `outcome`: forced class not in {absent, confused, inconsistent_reference} → `no_flag`; empty combined →
  `untestable`; `both models` with forced class `absent` → `lensing_residual_input`; `both models` otherwise →
  `second_model_supported`; anything else → `explained_model_dependent`.
- `final`: a `lensing_residual_input` or `second_model_supported` row is looked up in the Step 6 record below by
  system and position (within 1.0″). If the record gives an ordinary explanation that does not rest on image
  existence in a second model (`independent`), the row is `explained_step6`; otherwise `to_vet_candidate`.

**Step 6 record (frozen before the run; from `docs/fields/macs0717.md`).** `independent` = the ordinary explanation
holds without any second model predicting or not predicting the image.

| System, RA, Dec | Step 6 verdict | independent |
|---|---|---|
| 5, 109.379428, 37.746088 | copy of 5.1 (3.1″) | yes |
| 6, 109.374860, 37.770130 | copy of 6.3 (3.2″) | yes |
| 7, 109.366048, 37.765120 | position error, near-critical μ (D-036) | yes |
| 7, 109.365638, 37.764574 | position error, near-critical μ (D-036) | yes |
| 34, 109.373742, 37.744904 | copy of 34.3 (2.8″) | yes |
| 68, 109.395425, 37.739381 | on a critical curve, μ spread > 2× (D-036) | yes |
| 19, 109.409131, 37.754055 | copy of 19.1 (2.3″) | yes |
| 19, 109.397479, 37.747960 | central image below sensitivity | yes |
| 20, 109.371506, 37.763382 | below sensitivity (predicted S/N 5.8) | yes |
| 67, 109.377758, 37.768356 | counterpart within the position uncertainty | yes |
| 70, 109.403161, 37.746014 | below sensitivity | yes |
| 75, 109.382840, 37.744173 | μ spread > 2× (D-036) | yes |
| 75, 109.386305, 37.742128 | μ spread > 2× (D-036) | yes |
| 3, 109.393233, 37.747986 | CATS-only extra image | no |
| 19, 109.392674, 37.752049 | CATS-only extra image | no |
| 45, 109.384136, 37.743892 | CATS-only extra image | no |
| 68, 109.391572, 37.744043 | CATS-only extra image | no |
| 68, 109.384180, 37.744954 | CATS-only extra image | no |
| 71, 109.398568, 37.743334 | CATS-only extra image | no |

The 28 `confused` copy rows and the 65 pair are covered by the D-041 copy rule (`independent`: yes); a flagged row
that matches none of the table and is `both models` is matched against the copy rule by its Step 6 system list
(confused copies: systems 1, 3, 7, 12, 17, 24, 32, 34, 45, 50, 56, 62, 67, 71, 72, 75, 76, 79, 80; 65 pair) and
otherwise is `to_vet_candidate`.

**Comparison (diagnostic, not used for the outcome).** The six Step 6 "CATS-only" rows: does the pinned reading also
give `model-dependent (CATS only)`? Differences are listed.
