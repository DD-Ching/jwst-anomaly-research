# Plan: vet MIRI's known NED `G_Lens` object (TASKS Then 5, 2026-11-07)

The first end-to-end run (`20261007T020124Z-4bfabaaa`) recovered one NED `G_Lens` object in the MIRI top 20; it was
never vetted. That run predates jwst 3.0.0, so the object is re-identified on the current MIRI config
(`configs/smacs_jwst3_miri_confirm.yaml`): every top-20 row whose cross-match lists NED `G_Lens` within the
classify radius, excluding the five already vetted (`docs/candidates/smacs_miri_f770w_top20_singles.md`).

**Hypothesis (ordinary, to test first).** The object ranks high because it is a lensed (magnified) background
galaxy of SMACS 0723 with PAH/dust-bright MIRI colours, or because of a MIRI processing or catalog effect (edge,
low weight, deblending, bright-neighbour light, single-band merge radius).

**Lensing-consistency question.** If it is a catalogued multiple image, does its position agree with the
Mahler+2022 model images of its system (within the model's rms)?

**What would change the plan.** An object that survives A–C of `/vet-candidate` goes to a needs-human candidate PR.
If no top-20 row carries `G_Lens` in the current run, record that the Oct-7 object no longer ranks and close the item.

## Frozen rules for the MIRI flux-ratio step (written 2026-11-07T07:22Z, before any flux was read)

Run `20261107T071516Z-7e9ac450` reproduces the MIRI top 20 of #754; its only `G_Lens` row is `f770w_118`, already
vetted as Mahler+2022 system 3 (ICLv2 label 3.3). The open lensing question is therefore whether system 3, and any
other ICLv2 system with ≥ 2 MIRI counterparts, has MIRI flux ratios consistent with the model (lensing is
achromatic, so D-032's NIRCam result should carry over to 7.7–10 µm).

- Images: ICLv2 `arcs.dat` (60), μ from ICLv2 `best.par` at each catalogued position and the system z
  (`z_m_limit` as in D-032), no frame shift.
- Match: nearest MIRI merged source within 0.5″ (ASSUMPTION; MIRI PSF FWHM 0.27″ at F770W); a source matching
  several images makes them untested.
- Flux: `<band>_isophotal_flux` for F770W and F1000W, S/N ≥ 10; implied source mag
  `m + 2.5 log10 |μ|`; |μ| ≤ 20 (D-032).
- Residual: leave-one-out against the median of the other usable images of the system; `flux_outlier`
  |residual| > 0.75 mag (D-032 ASSUMPTION). Untested: isophotal area > 2000 px (ASSUMPTION; host/ICL light).
- Changes the plan: a flagged pair whose ratio is the same in F770W and F1000W (achromatic) and not explained by
  host light in the cutout goes to `/vet-candidate`; a chromatic one is a measurement failure (D-032).
