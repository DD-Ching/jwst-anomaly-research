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
