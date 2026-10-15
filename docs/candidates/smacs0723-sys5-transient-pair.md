# SMACS 0723 system 5 transient pair (`g03415` / `g03423`): vetting record

## 2026-10-15 (cloud routine, session 01Syioec)

- **Source:** two point sources 3.2″ apart on either side of the Mahler+22 system 5 arc, found by the two-epoch
  grid (`scripts/transient_grid.py`, SMACS ERO o001 vs VENUS o057; DECISIONS D-027 addendum 2). Image A
  (110.825136, −73.452178), image B (110.823747, −73.452975), J2000. `observed`.
- **Verdict: known astrophysical source, a multiply imaged supernova reported by the VENUS team as SN Atalanta.**
  Not a candidate. No lens-model prediction is unmet so far; the third image is untested (below).

### Tests

| Test | Outcome | Evidence | Label |
|---|---|---|---|
| A. Edge, low weight, spikes | pass | `quality_flag` `ok` in every epoch except B in the 2022 F150W mosaic (`low_weight`; flux still measured); no bright star within the spike zone | observed |
| A. Persistence, moving object | pass | on sky in 19 (A) and 15 (B) clean exposures over 1.6 h in six filters at fixed positions (`persistence_check.py`, PR #237) | observed |
| A. Single-band artifact | pass | detected in F150W and F444W in two separate programs (VENUS 2026-06-05, DD 12774 2026-07-06) | observed |
| B. Catalog association | pass | same sky position in all four epochs' level-3 mosaics; stamps looked at (`smacs_transient_epochs.png`) | observed |
| C. Known transient (TNS) | none at the position | TNS cone search, 120″ around the pair (2026-10-15): only AT 2024jfz at (110.83416, −73.451752), 9.3″ from A, a different source | observed |
| C. Known transient (programs) | **match** | JWST DD 12774 and HST 18297 ("Taking aim at Precision Cosmology with SN Atalanta", PI Agrawal) target `SMACS0723-GLSN-1` at (110.82377, −73.45298), 0.03″ from image B; DD 12774 includes a NIRSpec prism slit | observed (MAST metadata, STScI program pages) |
| C. Literature | **match** | Fujimoto et al. (arXiv:2609.30440, SN Helios) list SN Atalanta at z = 1.48 (DD 12774, PI Agrawal) among VENUS's strongly lensed SNe | observed (literature) |
| C. Light curve | consistent with a supernova | absent in 2022-06-07 and 2024-05-11 (≤ 0.011 μJy in r = 0.15″, the arc's light); 0.22 / 0.28 μJy (F150W, A / B) on 2026-06-05; 0.077 / 0.065 μJy on 2026-07-06 (fading by 3–4× in 31 d observed, ≈ 12.5 d rest frame at z 1.48); F444W 0.16 / 0.17 → 0.076 / 0.087 μJy | observed (`results/transients/smacs_grid/pair_epochs.json`; small-aperture fluxes, no aperture correction) |
| Lens model | consistent | Mahler+22 ICLv2 maps A and B to one source point at z 1.40–1.45 with opposite parities (`pair_trace.json`); the catalogued arc redshift is 1.425, Atalanta's reported z 1.48 | model_prediction |

### Open (not run)
- The ICLv2 third image (μ ≈ 3, about 28″ south, at (110.82190, −73.46032)) is not detected in VENUS o057. Its
  time delay is not computed (`lensmodel` has no Fermat potential), and at μ ≈ 3 it would be ~2.4 mag fainter than
  B, near the o057 depth. Without the delay this is untested, not unmet. The VENUS / DD 12774 team's time-delay
  analysis will settle it.
- Why the catalog stage missed the pair (DECISIONS D-027 addendum 2) is still open; it is a pipeline question, not
  a property of the source.

### Regenerate
```
python scripts/smacs_transient_epochs.py --out results/transients/smacs_grid/pair_epochs.json
python scripts/smacs_transient_pair.py --out results/transients/smacs_grid/pair_trace.json
```

## Follow-ups `g00791` and `g02290` (2026-10-15, same session)

Same script: `python scripts/smacs_transient_epochs.py --out results/transients/smacs_grid/followup_epochs.json
--position g00791=110.772443,-73.453621 --position g02290=110.670772,-73.497949`. Neither position is on the GO-4043
mosaics; `g02290` is also off the DD 12774 mosaics. Stamps looked at.

| Source | Light curve (r = 0.15″, μJy; observed) | Counterparts | Verdict |
|---|---|---|---|
| `g00791` | F150W 0.154 (2022-06-07) → 0.006 (2026-06-05) → 0.009 (2026-07-06); F444W 0.081 → 0.012 → 0.006 | TNS: none within 30″. SIMBAD: none within 3″. NED: nearest a galaxy at 1.9″ (NDA2023 00711, z 1.883, flag S1L) | A compact source on the outskirts of a galaxy, seen only in the 2022 ERO (on sky in 50 clean exposures, PR #237). Most likely an unreported supernova in the ERO epoch (`hypothesis`; no spectrum, one epoch). Ordinary class; not a candidate. Its host's redshift and any lensed counter-image are not checked |
| `g02290` | F150W 0.451 → 1.052 (×2.3); F444W 0.193 → 0.345 (×1.8); 2022 and 2026-06-05 only | TNS: none within 30″. Nearest SIMBAD/NED entry, F150DA-005 (photo-z 3.01), is 2.3″ away and is a different source | A brightened compact galaxy nucleus. Most likely AGN variability, otherwise a nuclear transient (`hypothesis`). Inconclusive without a third epoch; not a candidate |
