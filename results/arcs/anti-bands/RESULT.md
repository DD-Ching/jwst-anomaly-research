# `arcs`: second-band shapes of the 3 single-component pooled anti rows (exploratory)

Plan: `PLAN.md` (committed before any second-band shape of these rows was read). Script:
`scripts/arcs_anti_bands.py`; rule functions unit-tested before the run (`tests/test_arcs_anti_bands.py`). Every
number below is in `summary.json`. Inputs: six pipeline `_cat.ecsv` from `stpubdata` S3 (sha256 in `summary.json`
`inputs_sha256`), read into the cache and deleted after the run. Shapes `observed`; classes and labels `derived`;
thresholds ASSUMPTIONs (`summary.json` `params`).

## Frozen outcome: **blend veto passed; the rows stay open (not a lensing test)** (rule 4: 3 `same_axis`, 0 `blend_flag`)

| Band gate (rule 1) | Control usable | Stay aligned | Δpa median / 95th pct | Passes |
|---|---|---|---|---|
| SW (F150W, 3 fields) | 30 | 29 | 1.7° / 11.6° | yes |
| LW (F444W; F410M for MACS0308) | 28 | 28 | 1.1° / 7.8° | yes |

| Row | F200W pa, e | F150W Δpa, e, S/N | LW Δpa, e, S/N | Class in SW / LW | Label |
|---|---|---|---|---|---|
| WHL0137 2194 | 70.8°, 0.795 | 0.2°, 0.821, 81.7 | 0.1°, 0.651, 99.5 (F444W) | anti / anti | same_axis |
| RXC J0600 581 | 30.6°, 0.688 | 0.1°, 0.691, 417.7 | 0.3°, 0.631, 618.7 (F444W) | anti / anti | same_axis |
| MACS0308 2157 | 176.4°, 0.531 | 4.4°, 0.487, 285.3 | 3.7°, 0.481, 173.3 (F410M) | anti / anti | same_axis |

All six band matches lie within 0.07″ of the F200W position; no match is shared.

## Contact sheet (`anti_bands.png`, looked at; reading by eye, post hoc, never used for a label)
- 2194: a thin straight streak with a faint knot near its centre, the same in all three bands; a fainter blob and a
  bright compact neighbour (each roughly 1–1.5″ away by eye, not stored) are separate, not blended into it.
- 581: a near edge-on disc with a central bulge in all bands; no second component.
- 2157: one compact elongated galaxy; a fainter neighbour (roughly 1.5″ away by eye, not stored) is separate.
- The red bar (F200W catalogue axis, placed through each cutout's WCS) follows the light in every panel.

## Reading (narrowest statement)
- In both passing bands the 3 rows keep their F200W axis to within 4.4° (control 95th percentile 11.6° / 7.8°), and
  each stays `anti` against its deciding model. A colour blend or colour-dependent structure does not set these axes.
- Band agreement cannot separate an intrinsic shape from a lensing stretch (addendum 4): both are achromatic. The
  rows stay open; the post hoc within-class excess (`../pooled-census/posthoc.json`, P 0.024) is neither explained
  nor strengthened. An intrinsically elongated galaxy (e.g. the edge-on disc 581) whose own axis dominates a weak
  shear remains the leading ordinary reading; it is not tested here. No candidate.

## Deviations
- The LW band of MACS0308 is F410M, as the plan states (no F444W catalogue for `jw05293-o005_t001`).
- `scripts/arcs_anti_bands_sheet.py` (contact sheet) was added after the run; it reads `summary.json` only for titles.

## Limitations
- 3 rows, chosen post hoc. Pipeline isophotal shapes (no PSF correction): the LW ellipticities are lower (PSF
  dilution), the axes are not.

## Next
The deciding models' local shear error at the 3 rows (model-to-model spread of the predicted axis and |g|), and
whether intrinsic shape noise at the observed e (0.531–0.795) and the deciding g_min (0.24–0.64, `g_min_200`) explains an anti
axis; membership only with a spectrum.
