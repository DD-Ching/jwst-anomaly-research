# W1 calibrated streak veto: result (2026-10-17)

Plan: `PLAN.md` (committed before the runs; addendum A, a revised mask normalisation, committed after looking at
the first two masks and before any screen or injection run: **deviation** from the first frozen rule). Labels: pixels
and shapes `observed`; masks, S maps and limits `derived`; injections `simulated`; every threshold an ASSUMPTION.

## Masks (`masks/<field>_streak_mask.{fits.gz,json,png}`; `scripts/w1_streak_mask.py`)
STScI `acstools.findsat_mrt` (median Radon transform) on the streamed F200W mosaic (block 4, clip ±10σ, footprint
eroded 1.5″, rebin 2; status-2 trails only). Persistent trails per field: Abell 2744 0, MACS0416 2, MACS1149 0,
Abell 370 2, Abell S1063 3, SMACS 0723 3, El Gordo 0, MACS0417 4, WHL0137 1, RXJ2129 1, MACS0647 9, MACS1423 3.
Sheets looked at: the trails are vertical stripes and diagonal streaks or spikes. Missed by eye: a faint diagonal
streak in MACS1423, the fainter crossing streak at the MACS0647 S_max, and faint parallel diagonal streaks in the
S1063 S_max aperture (`sheet_abells1063_veto_smax.png`, from the superseded run; same S_max position). False trails: the MACS0647 mosaic gap and halo (PLAN addendum A).

## Known cases and adoption (`masks/adoption.json`, computed by `w1_streak_mask.evaluate`)
- K1: the MACS1423 mask covers the 2026-10-17 hit position: **yes**.
- K2: MACS1423 with the veto: S_max 3.41, p_rot 0.22, b_max 3.64, D-053 rule **not passed** (was 3.88, 0.01, passed).
- **Adopted** (K1 and K2 hold).

## Screen and injections with the veto (`inject_summary_<field>.json`; seed 50, 200 per mass)
The veto blanks the rows on trails inside `shear_screen` after the spike veto, PSF size and responsivity are
estimated from the whole table; vetoed rows are not lensed in injections, and painted images on trails are vetoed.
**Deviation (post hoc, after code review):** the first run dropped the rows from the table before the screen, which
also changed the catalogue spike veto (Abell S1063 then passed at p_rot 0.045; that run is superseded and not stored).



No recoveries at 2 × 10¹⁰ or 2 × 10¹¹ M☉ in any field. No field passes the D-053 rule with the veto. **MACS0647 loses
25.2 % of its used rows (> 20 %, PLAN):** combined, and named here; its trails include false ones along the mosaic gap
and halo (addendum A). Against the earlier summaries, recoveries move both ways (MACS1423 9 / 72 / 140 → 10 / 60 / 134;
MACS0647 0 / 1 / 1 → 0 / 5 / 5); the changes also include different rotation-null draws, so the veto's own cost is
not separated.

## Limits (`combined_12fields_veto.json`; `combined_11fields_veto_no_macs0647.json` without MACS0647)
- **Headline, twelve fields with the veto (MACS1423 back in), 70.5 arcmin²:** Σ(W1) < 4.4 × 10³ / 7.2 × 10² /
  3.5 × 10² deg⁻² at 2 × 10¹² / 8 × 10¹² / 2 × 10¹³ M☉ (border-corrected 4.9 × 10³ / 8.1 × 10² / 3.9 × 10²).
  Without MACS0647: 4.4 × 10³ / 7.3 × 10² / 3.5 × 10².
- Earlier eleven fields, no veto (`results/w1_shear_rxj_m0647/combined_11fields.json`): 4.5 × 10³ / 7.6 × 10² /
  3.8 × 10². This is the first W1 limit that runs a streak test through the injections.
- Still blind at ≤ 2 × 10¹¹ M☉.

## Limits of this result
- The mask finds bright and long trails. It misses faint parallel streaks (S1063 S_max aperture, by eye).
- The veto rule was revised once after the first two masks were looked at (addendum A).
- p_rot rests on 200 rotation draws.
- Inputs were staged by `scripts/w1_stage_inputs.py` (sha256 checked where recorded; the RXJ2129 photo-z tarball
  has no recorded hash, its member does).
- The mask FITS files are not tracked (`*.fits.gz` is gitignored): re-make them with
  `scripts/w1_streak_mask.py --fields <field> --out results/w1_streak_veto/masks`; each veto-run summary records the
  mask's sha256 and pixel count (`assumptions.streak_mask`). The committed runs used masks written with a gzip
  timestamp, so their file sha256 does not reproduce byte for byte: check a re-made mask by `mask_pixels`. Masks
  written now carry no timestamp, and summaries also record `data_sha256` (shape + pixel bytes).

## Second code-review (final diff, 2026-10-17)
Fixed without changing any number above: deterministic mask files and a pixel hash; `shear_screen` checks the
veto's dtype and length; `OMP_NUM_THREADS=1` for the trail-finder pool; no HEAD request before staging (the byte
count comes from the kept file); a stale comment; a test of the veto path (blanked rows, background, painted images).
Traced, no change:
- Aperture coverage and the injection footprint still count vetoed rows. Injected images that land on a trail are
  vetoed and count as not recovered, so the efficiency, and through it the limit, already pays for the vetoed area;
  only the `min_coverage` cut on grid centres is slightly generous next to trails.
- PSF width and responsivity R are estimated before the veto, over every row (frozen in PLAN, so the spike veto is
  unchanged). The injections use the same R, so recovery is self-consistent; an R measured without streak rows is
  untested.

## Open code-review findings (not fixed; no effect on the numbers above)
- The veto distance is measured from the nearest mask-pixel centre, so it is exact only to about half a mask
  pixel (≈ 0.18″) (documented in `veto_function`).
- `w1_stage_inputs.py` fetches the files without a recorded sha256 through a sentinel that compares equal to any
  hash. The fetch helpers should get an explicit no-verify path instead.
- `w1_streak_mask.py` streams each mosaic one strip at a time, with no concurrent range reads and no
  download/compute overlap (scripts/CLAUDE.md parallel topology).
