# El Gordo arc difference imaging, PEARLS 1176 o241 → VENUS 6882 o051: plan (written before any run)

**Status of this unit.** This is an owner-approved, one-box exploratory pilot (2026-10-16), with a 4 h budget
from 22:14Z.
- The rule "full-chain recovery ≥ 0.50 before extending the search" (#277, failed on SMACS at 0.34) is waived for
  this box only. That failure record stands, and the waiver does not authorise a wider scan.
- Detection, data-quality and reality checks are unchanged. Known output-changing bugs must be fixed.
- Earlier El Gordo work used catalogues, lens models and arc shapes, not two-epoch differences of this box, so this
  is a first difference search here. It is still **exploratory, not a blind test**, because the tool settings come
  from earlier fields.

## Data (observed)
- **Epoch 1:** `jw01176-o241_t012_nircam_clear-{f150w,f444w}_i2d.fits`. PEARLS GTO 1176, 2022-07-29; F150W
  1890 s, F444W 2491 s.
- **Epoch 2:** `jw06882-o051_t051_nircam_clear-{f150w,f444w}_i2d.fits`. VENUS 6882, 2026-08-05; 1245 s per band.
- **Access:** S3 byte ranges, as `arc_diff.py` does. The pipeline versions are taken from the headers and stored in
  the summary.
- **Coverage:** these are the only two public NIRCam epochs at El Gordo (MAST query 2026-10-16). F150W and F444W
  are the only common bands; VENUS also has F210M and F300M.

## Box (chosen from catalogued images and coverage only, before looking at any difference)
- 24″ × 24″ centred on (15.7300, −49.2530).
- It holds 11 Caminha+2023 multiple images of 6 systems, the densest such box fully inside both epochs' F150W and
  F444W footprints.

## Observable and link to the project's hypotheses
- **Observable:** a compact, sky-fixed flux change (either sign) between 2022 and 2026 on or near lensed arcs.
- **Ordinary explanations:**
  - supernovae and other transients in lensed or cluster galaxies;
  - caustic microlensing of lensed stars near critical curves;
  - AGN or stellar variability.
- **Instrumental explanations:**
  - kernel residuals on bright cores (PSF-rotation rosettes);
  - diffraction spikes;
  - registration dipoles;
  - detector persistence;
  - mosaic edges;
  - pipeline-version differences between epochs.
- **Link (hypothesis only, `docs/exotic_lensing.md`):** a negative-mass or Ellis-type microlens on a caustic-
  crossing arc would make a compact image vanish or dim between caustic spikes (W3). Ordinary caustic
  microlensing produces both brightenings and fadings, so a single-epoch change is never evidence for it.

## Chain (existing tools, settings as SMACS #275 / review list #277; nothing tuned here)
1. `arc_diff.py --band f150w … --band f444w … --ra 15.7300 --dec -49.2530 --half-arcsec 12 --kernel-px 4
   --kernel-stamps --inject 100 --min-bands 2 --procs 4 --maps <scratchpad>` (per-band point-injection recovery
   at 27.5 / 28.0 / 28.5 AB).
2. `arc_diff_morph.py --all-hits`, with the default verdicts.
3. `arc_chain.py review` with position labels only; nothing is dropped for its position. `arc_chain.py` gains an
   `elgordo` image list: Caminha+2023 `elgordo-caminha23` images, shifted to the JWST frame as in
   `lens_consistency.py`.

## Controls (before the screen counts)
- **No known real transient** in this box is in the repository, so no real positive control exists for this pair.
  The SMACS ERO → VENUS pair (same instrument, the same VENUS program and bands) recovered SN Atalanta. That is a
  development control from another field, not a control here.
- **Applicable control (full-chain injection in this box):**
  - `arc_diff.py --inject-grid 80` with `--inject-grid-mags 26.0 27.5 28.0`, seed 0.
  - Gaussian point sources are added to epoch 2 of both bands on pixels above the epoch-1 90th percentile, 1.5″
    apart. These are galaxies, arcs and ICL; there are no blank-sky or bright-core-PSF positions.
  - Count = hit, unflagged and morph `pass`, by `arc_chain.py inject-review`.
  - **Gate:** whole-box efficiency ≥ 0.80 at 26.0 AB. If it fails, the box is a tool failure and the search stops.
  - 27.5 and 28.0 AB are reported with no threshold. The SMACS 0.34 is never applied here.
- **Data-quality gate:** both epochs have finite data over ≥ 95 % of the box in both bands, and `arc_diff.py`
  completes. Otherwise stop.

## Review and vetting (fixed now)
- **Review burden cap:** ≤ 15 review-list entries in total. If there are more, the list is saved and the run stops
  without vetting; no new rule is added to cut it.
- **Every entry gets the existing checks:**
  - the sheet (`arc_diff_sheet.py`);
  - `persistence_check.py` where it applies;
  - multi-epoch forced photometry with `hit_lightcurve.py` (every covering public mosaic; conservative errors);
  - a sky-fixed check: same position in both bands within 0.1″ and present in the epoch's two bands;
  - nearest catalogued object and known-transient name (DJA / VizieR at ≤ 1″ where reachable).
- **Classes (per entry; the rule is applied from these measurements):**
  - `artifact`: morph veto on the sheet, or a persistence / spike / edge / registration explanation found.
  - `ordinary transient or variable`: sky-fixed and real, with a host or a point source in which an SN, AGN or
    star explains it, or a multi-epoch light curve that does.
  - `open`: real and sky-fixed with no ordinary explanation yet. Only this class goes to `/vet-candidate`.
- **Worth testing a special model:**
  - **Only** an `open` entry on a lensed arc, where the lens model's counter-images and time delays can be checked.
  - Even then, the next step is a check, not a claim.
  - An SN or variable star is recorded but is never wormhole evidence.

## Outcomes
- **Zero entries, or only artifacts and ordinary sources:** the limited-sensitivity search found nothing in this
  box, at the injection efficiencies measured. No exclusion or abundance limit is derived.
- **Separately:**
  - MACS1423 (#282, merged) and #279 used the pre-fix `fit_models`. They are flagged for a re-check in the state
    files. They are not re-run here and not used for any conclusion.

## Addendum 1 (2026-10-16T22:35Z, after the review list, before the on-spike controls are measured)
- **Found in review:**
  - All 4 review entries lie within 2.6–3.9″ of the box's brightest point source, along its diffraction spikes
    (`spike/spike.json`).
  - PA_V3 differs by 0.84° between the epochs (headers), so a spike moves sideways by about 0.05″ at 3.5″.
  - The frozen dipole test of `scripts/elgordo_spike_check.py` (±1.5 px across the spike) passes only for h043.
- **On-spike control (rule fixed now):**
  - Take 8 positions on the same spike as the entry: the entry's angle from the source, distances 2.0–5.0″ in
    equal steps, skipping any within 0.3″ of the entry.
  - Measure them with `hit_lightcurve.py` like the entry (same aperture, conservative errors).
  - The entry is `explained by the spike shift` if ≥ 4 of the 8 controls show a change of the same sign as the
    entry, at ≥ 2σ, in the band where the entry's change is largest (F444W for h049). Otherwise it is
    `undetermined`, and an `undetermined` entry is treated as `open` for the vetting step.
  - Applied to every entry whose own forced-photometry change is ≥ 3σ in any band:
    - h049 (F444W −3.6σ);
    - h043 (F444W −4.1σ, on a bright galaxy; its rotation dipole is confirmed in F150W only).
    - h024 and h031 have no ≥ 3σ change in any band.
  - (Corrected before any control was measured: the first wording left out h043.)
