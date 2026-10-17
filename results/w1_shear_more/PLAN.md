# W1 shear screen on MACS0417, MACS1423 and WHL0137 (plan, frozen before the run)

Cloud routine, 2026-10-17 (claim PR #305). TASKS "JWST focus" 1 (W1): more fields with public photo-z. MACS0717
still has none (DJA v7 index checked 2026-10-17). Three fields that already have a JWST F200W pipeline catalogue,
a photo-z catalogue and a map model that passed the `arcs` alignment control:

| Field | Catalogue (sha256) | Photo-z (sha256) | Model (z_l) |
|---|---|---|---|
| MACS0417 | `jw01208-o006_t003` F200W (`e55dffb9…`) | CANUCS DR1 `macs0417-clu` (`ca6b48aa…`) | `macs0417-canucs-maps` (0.443) |
| MACS1423 | `jw01208-o010_t005` F200W (`1a968574…`) | CANUCS DR1 `macs1423-clu` (`c257d8df…`) | `macs1423-canucs-maps` (0.545) |
| WHL0137 | `jw02282-o010_t001` F200W (`b0693065…`) | DJA v7.5 sunrise zout (`672d4bdf…`) | `whl0137-relics-lenstool` (0.566) |

All sha256 match `results/arcs/*/RESULT.md` and SOURCES.md. Maps and photometry streamed to the gitignored cache,
deleted after the run.

**Hypothesis (H).** A negative-mass (W1) lens of |M| ≥ 2 × 10¹² M☉ in a field gives a radial (E-mode positive in
the D-053 sign convention) aperture-mass peak.
**Ordinary explanations of a peak.** Diffraction spikes and streaks, PSF anisotropy and isophotal-moment noise
(B-mode floor), blends near bright members, residual cluster shear from model error (R g removal; these models are
preliminary or HST-era), edge apertures.
**Rule (D-053 defaults, unchanged; ASSUMPTIONs).** Schirmer Q_TANH, R_ap 10″, 1″ grid, sources 1–10″, rotation
null; a peak counts only if p_rot < 0.05 **and** S > the field's max |S_×|. Grid half-width 100″ (as MACS0416 and
MACS1149, the other CANUCS fields; ASSUMPTION for WHL0137 too). Spike veto: catalogue segments plus Gaia DR3 G < 17
within 4′ of the model reference point.
**Injections.** `inject_shear.py --fields <field>` seed 50, 200 per mass at the five D-053 masses, then a ten-field
combination with `--combine-only`.
**What changes the plan.** A peak passing the rule goes to `/vet-candidate` (streak check on the image first). R
outside (0, 1.5] or < 20 calibrating rows: the field is dropped (D-053). Null fields: ten-field limits replace the
seven-field ones.
**Known-case check.** R should fall near the 0.41–0.49 of the other seven fields (the cluster shear is seen); a
field with R far outside it is reported and not combined.

## Addendum A (written after the screen, before the vetting computation; 2026-10-17 11:18Z)
MACS1423 passed the D-053 rule (S_max 3.88 at (215.93975, 24.06457), p_rot 0.01, B-mode max 3.64). The 30″ F200W
sheet (`sheet_macs1423_smax.png`, looked at) shows a vertical stripe crossing the aperture centre and a bright
star (no Gaia G < 17 cross) about 8″ from it with spikes, one pointing near the centre. Frozen vetting test
(`scripts/w1_peak_vet.py`, computable):
- V1 streak: in the 30″ cutout, stripe columns are those whose column median deviates from the cutout median by
  > 5 × the MAD of the column medians; used sources within 1.0″ of a stripe column are removed.
- V2 star: the brightest catalogue source within 10″ of the peak is the star; used sources within 15″ of it whose
  catalogue axis is within 15° of the direction to the star are removed.
- S at the peak is recomputed without V1 ∪ V2 (and without each alone, reported). **Outcome rule:** if S without
  V1 ∪ V2 is ≤ the field's B-mode max 3.64 or its p_rot (same null draws) ≥ 0.05, the hit is "explained by
  streak/spike contamination (screen-level)"; otherwise it goes to `/vet-candidate`. All thresholds ASSUMPTIONs.
