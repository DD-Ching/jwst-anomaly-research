# SLICE JWST-era images vs RELICS Lenstool v1 maps (RXC J0232, RXC J2211, MACS J0553) — frozen plan

Frozen and pushed before any model map is read (2026-10-25; the maps were downloaded only to pin their sha256).
Rules mirror `results/lens_consistency/macs2129_clash_heldout/PLAN.md` (D-034 addendum bb), which reuses the
MACS0647 rules (addendum v); the differences are listed. Script `scripts/slice_relics_heldout.py`, rule and
parser functions tested in `tests/test_slice_relics_heldout.py` (next commit).

## Question and hypotheses
- **Question:** do the HST-era RELICS Lenstool v1 deflection maps of three SLICE clusters reproduce
  (a) the source galaxies SLICE first identified in JWST imaging, and (b) the spectroscopic families at their
  spectroscopic z?
- **Ordinary expectations:** an HST-era model built on 3–10 families misses some JWST-new images by more than its
  rms (incompleteness, substructure, a free z); a miss is a lensing-residual *input*, never a candidate.
- **What would change the plan:** a model failing gate K0 or check F is no usable reference and gives no reading
  (as on MACS J2129). A held-out miss by a usable model → a lensing-residual input to check against the SLICE
  Lenstool model (LAM atlas, closed on 2026-10-25) before any interpretation.

## Inputs (observed)
SLICE (Cerny et al. 2026) Table 8 blocks, parsed from the e-print `main.tex` (sha256 f7120d94…a9a7) into
`slice_images.ecsv`. **Clump 1 only:** of each source, the rows `S.1.k` (one point per image; the extra clumps
of a source are not independent constraints). Candidates (`c…`) are not in these blocks.

| Cluster (z_lens, SLICE Table 1) | RELICS map | Held-out families (paper text) | Spectroscopic families (reading B) |
|---|---|---|---|
| RXC J0232.2−4420 (0.2836) | `rxc0232m44_lenstool_v1` | none stated → no reading A | 5 (1.645), 7 (1.467) |
| RXC J2211.7−0349 (0.397) | `rxc2211m03_lenstool_v1` | 4–10 ("only visible in the JWST imaging"; 1–3 = Cerny+2018) | none: family 1 is 1.051 in Table 8 but 1.501 in the text, so it is not used |
| MACS J0553.4−3342 (0.412) | `macs0553m33_lenstool_v1` | 11–21 ("only fully confirmed using the JWST imaging"; 1–10 = Ebeling+2017) | 1 (1.14), 5 (3.32), 10 (3.76), 13 (3.35) |

Maps (x / y arcsec deflection at D_LS/D_S = 1, the RELICS Lenstool convention already used for MACS0417 and
SMACS 0723): sha256-pinned in `jwst_anomaly.lensmodel`; cosmology (70, 0.3) as for the other RELICS models.

## Rules (ASSUMPTIONs; per cluster, one model each)
1. **Redshift rule (fitted z):** as MACS0647 rule 1 (60-step log grid z_lens + 0.2 … 12; `m0417_images.best_z`).
2. **Gate K0** (the families not held out, fitted z): rms ≤ 2.0″ **and** matched fraction (≤ 3″) ≥ 0.8.
3. **Frame check F** (same set): |median| observed − nearest predicted ≤ 0.5″ in RA·cos δ and in Dec.
4. **Reading A** (held-out families, fitted z): an image whose nearest predicted image is > 3″ away gets
   `lensing_residual_input`, only if the model passes K0 and F.
5. **Reading B** (spectroscopic families at their z_spec, no fit, same 3″ match): a miss gets
   `redshift_residual_input`, only if the model passes K0 and F.
6. Families with more or fewer predicted than catalogued images are listed (not read).

## Not done here
No imaging, photometry or extra-image search. SMACS J2031 (no RELICS maps on MAST) and the other SLICE clusters
(no public HST-era maps found) are not tested. The SLICE Lenstool models stay the JWST-era reference to re-run
against when the LAM atlas reopens.
