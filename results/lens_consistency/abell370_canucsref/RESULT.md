# Abell 370: CANUCS DR1 reference extra-image test, gate K (result)

Plan: `PLAN.md` (frozen, commit 99c9738, before any run). Script: `scripts/a370_canucsref.py`. Run 2026-10-22, cloud,
1 min 50 s. Outputs: `summary.json`, `imageplane_abell370_abell370-canucs.ecsv`, `readings_abell370.ecsv`,
`canucs_a370_images.ecsv` (the model's 115 images, z fixed).

## Result (frozen rule)
- **Gate K fails** (derived): image-plane rms **2.34″** on 115 / 115 images (limit 1.0″); 95.7 % reproduced within
  3.0″ (limit 80 %, passes). It fails on rms alone. This reproduces D-044's `validate` rms 2.32″, as the plan expected.
- By rule 1 the model is not used: **no counter-image prediction, CANUCS flux or pixel was read**; no flag; nothing to
  vet. The Abell 370 line of the CANUCS-reference route stops here.

## Reported, not rules (post hoc)
- Five images carry the rms: 7.1 (18.1″), 5.3 (10.7″), 3a.2 (9.7″), 14a.3 (5.1″), 5.2 (4.0″); next is 14a.2 (2.3″).
  Without them, rms 0.60″ on 110 images (98.2 % within 3″). Leave-one-out on 29b and the family refits are in
  `summary.json` (`post_hoc`).
- Reading (hypothesis, ordinary): the model was fitted in the source plane (D-044), which does not penalise
  image-plane offsets of single images; large offsets of a few images of a source-plane fit are expected and are not
  a lensing residual. Not tested further here.
- Why the post hoc pass is not used: the five images were chosen after seeing the residuals, and the extra images
  that matter most for the test are predicted near the same caustics (36 families predict more images than
  listed). A route that keeps only fully reproduced families would need its own frozen plan and a known-case check.

## What would reopen it
An Abell 370 model fitted in the image plane with a public parameter file or maps, or a frozen family-restricted rule
validated first on MACS0416 / MACS1149. Checked after the run (2026-10-22, not yet pinned or validated):
- BUFFALO Niemiec et al. 2023 (MNRAS 524, 2883) model A: Lenstool image-plane fit, quoted rms 0.90″ on 98 gold
  images (paper, read by a research subagent). The MAST HLSP (DOI 10.17909/t9-w6tj-wp63) download script
  `hlsp_buffalo_hst_model_abell370_multi_niemiec-a_v1.0_lens-models-download.sh` lists x/y deflection maps (`main`,
  `par`, `full` variants), κ, γ, μ(z 1, 2, 4, 9) and `..._multi_niemiec-a-sl_v1.0_input.par` (listing read here;
  whether the image list with z is included is unverified). Best next model for gate K.
- Fallback: Lagattuta et al. 2019 (CDS `J/MNRAS/485/3738`): 138-image list (123 z_spec) and Gold_Shear dPIE
  parameters (image-plane rms 0.66″), but no member catalogue, so only an approximate rebuild.
