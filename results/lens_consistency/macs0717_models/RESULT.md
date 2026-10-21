# MACS0717: Sharon v4cor and Keeton v4 as map models — Keeton passes both gates; Sharon fails G2 as written (6.2 off the map)

Plan: `PLAN.md` (frozen before the run). Script: `scripts/m0717_models.py` (gates from `scripts/m1149_sharon.py`).
Numbers: `summary.json`, `imageplane_*.ecsv`. All model outputs are `model_prediction`; the image positions and
redshifts are `observed` (Sharon readme Table 7, transcribed into `sharon_gold_images.ecsv`, 26 images, 8 systems).

| Model | G1 κ median rel. diff (pixels) | G1 \|μ(z=2)\| median (p95; pixels) | G2 rms (solved) | G2 systems | Frozen verdict |
|---|---|---|---|---|---|
| `macs0717-keeton` | 0.00002 (11,578) | 0.00002 (0.0001; 10,308) | 0.853″ (26 / 26) | 8 / 8 | **pinned for use** |
| `macs0717-sharon` | 0.00025 (11,571) | 0.024 (0.119; 9,933) | 0.639″ (25 / 26) | 7 / 8 | **G2 fails as written** |
| `macs0717-cats` (comparison, no gate) | — | — | 2.77″ (26 / 26) | — | — |

**Sharon G2 failure.** Image 6.2 (109.362705, +37.752681) lies about 101″ west of the map centre, just outside the
200″ × 200″ maps (`summary.json` `post_hoc.systems_outside_maps`), so system 6's source is traced from 6.1 and 6.3
only and 6.2 cannot be solved. A map-footprint limit, as for MACS1149 system 9 (D-037 addendum 2026-10-20), not a
test of the model there. **Post hoc (deviation):** on the 7 systems with every image inside the maps, G2 passes
(rms 0.662″, every image within 3″). So `macs0717-sharon` is used only inside its map footprint (post hoc condition).

**Keeton** reproduces every image within 1.8″ (largest 1.1 1.79″, system 6 1.50–1.69″); its rms 0.853″ is close to
its readme's quoted 0.79″ on 29 images (it also fits system 15).

**CATS comparison.** CATS v4.1 places the same 26 images at rms 2.77″ (largest: system 3, 4.7–5.6″; 14.3 3.8″; 6.2
3.9″), consistent with its 3.2″ on its own 132-image list (D-035). Both independent models fit the gold set 3–4×
better, as the D-041 vetting assumed when it used them in scratch code.

**Limits.** The known cases are Sharon's own constraints (and largely Keeton's gold list), so G2 checks that the maps
load correctly, not that they predict independent data. Gold spectroscopic systems only.

**Next.** Re-run the D-041 cross-model readings of the MACS0717 forced rows with the pinned `macs0717-keeton` (and
`macs0717-sharon` inside its footprint), replacing the scratch readings in `docs/fields/macs0717.md` Step 6.
