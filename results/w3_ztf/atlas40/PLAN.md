# W3 ZTF: ZTF alone vs ZTF + ATLAS, paired injections on the 40 calibration stars (feasibility trial)

Owner brief 2026-10-22 after #461:
- 4 h wall clock, ≤ 2 GB, ≤ 4 cores. Session start 08:56 UTC.
- ATLAS through the documented shared Astro Data Lab token. At most 2 of our tasks queued or running; only our own tasks
  are touched. No account, no e-mail.
- A feasibility trial on the already-used calibration sample, not an independent validation. No search group, no
  search.
- X4, ASAS-SN and the other paused units are untouched.

## 0. Prior state (checked 08:56–09:00)
- `main` at b85ca54 (#461 merged). Open PRs: #442 (dependabot) and the `needs-human` backlog; none on W3 / ZTF or ATLAS.
- #461 kept only 498 in-gap epochs of 4 stars (`../atlas/analysis.json`). The full light curves were deleted, so they
  cannot be recomputed and all 40 stars are requested again. The 4 #461 stars are **development data**:
  601209200002983, 601213300003718, 632215200005622 and 632211200001819.
- The official pages (output description, FAQ, API guide, schema), re-read 08:57, still say nothing about the value
  9999. The NOIRLab notebook still allows the Data Lab token.

## 1. Acquisition (committed before any request; `targets_script.py.txt` → `targets.json`, `acquire_script.py.txt`)
- **Targets:** the 40 stars of `../asassn/selection.json`, one request each, reduced mode (total flux), MJD 58150–61000.
  That covers both areas' ZTF carrier ranges (58198.4–60969.3).
- **Priority order** (committed counts only): most A0-insufficient W3 events, then most W3 events, then most fades,
  then key. No star is swapped for ATLAS coverage or results.
- **Limits:**
  - ≤ 2 tasks in flight, polled every ≥ 30 s; 429 honours `Retry-After`.
  - Never resubmitted.
  - A task unfinished 45 min after submission is recorded as `timeout` and left alone.
- **Acquisition deadline:** 170 min from the start of the scheduler. Unfinished tasks are then recorded for
  continuation, not as "no data".
- **Kept for recomputation:** compact per-star measurement tables (< 1 MB each; the columns needed by the quality cuts,
  time and flux), task ids, time stamps and source sha256. The token is in no product.
- The method (quality, time, baselines, errors, injection, rules, denominators, success and stop) is fixed in §2–§5,
  committed before any ATLAS flux of this round is read.

## 2. ATLAS measurements (fixed before any flux of this round is read; `prep_script.py.txt`)
- **Product:** the official reduced-mode total flux (uJy, duJy). Difference flux is not used.
- **Quality:** the FAQ filters, as #461:
  - duJy < 10000, err == 0;
  - 100 < x, y < 10460;
  - 1.6 < maj, min < 5;
  - −1 < apfit < −0.1;
  - mag5sig > 17, Sky > 17;
  - filter c or o only.
  A failing row is unknown: neither dim nor normal, and never zero flux. Low-S/N and negative fluxes are kept with
  their errors.
- **The value 9999:** no official page defines it. In #461 the one row with uJy = 9999 had m = 13.900 =
  −2.5 log10(9999) + 23.9 exactly (seen in the session, not stored). That looks like a 4-digit output limit. Rule,
  fixed now: rows with |uJy| ≥ 9999 are **unknown (possibly clipped)** and are excluded like failing rows. They are
  counted. A sensitivity run keeps them (secondary, not gating).
- **Time:** MJD (exposure start, UTC) + 15 s, converted to heliocentric MJD at the star, as #461 and as the ZTF
  series.
- **Baselines:** per star, filter (c, o) and ATLAS unit (`Obs` prefix, e.g. 01a): the median uJy of the passing rows.
  A (filter, unit) with fewer than 20 passing rows uses the star's per-filter median. Relative flux = uJy / baseline;
  σ = duJy / baseline. ATLAS c and o, and ZTF g and r, each keep their own baseline. Similar magnitudes are not taken
  as proof of a common calibration.
- **Errors (fixed now):** #461 found the scatter 1.7–4.4× the quoted errors.
  - Per star, filter and unit: k = the normal-consistent MAD of (rel − 1) / σ over the passing rows, using
    `jwst_anomaly.rank.robust_location_scale`, the repository's robust scale.
  - σ_cal = max(k, 1) × σ. The fallback for < 20 rows is as for the baselines.
  - Errors are never shrunk. A sensitivity run uses the quoted σ (secondary, not gating).
- **Kept:** `phot/<key>.csv.gz` per star (MJD, filter, uJy, duJy and every column the cuts use, plus Obs) and
  `calib.json` (baselines, k, counts). Each table must be < 1 MB.

## 3. Injection (fixed now)
- **Events and ZTF:** the replayed v23 events are unchanged: same samples, carriers, seeds (offsets 45000 / 46000) and
  ZTF noise.
  - A is ZTF-only bounds and must equal the stored v24 outcome for every event of both areas. If any differs, stop
    and find out why.
- **The same physical event on the ATLAS epochs:**
  - W3: `w3_dmag` with the same t0, t_E, u0, ρ and the cell's blend.
  - Fades: the same `control_model` closure.
  - ASSUMPTION: lensing is achromatic, with the same source fraction in ATLAS c / o as in ZTF. ATLAS's PSF (FWHM
    ≈ 3.7″) is wider than ZTF's; the 40 stars pass the 16″ clean cut, but extra blend light and chromatic blends are not
    modelled.
- **Flux injection,** reusing the W3 rule of `w3_microlensing.inject_w3` (OGLE):
  - F = 10^(−0.4 dm(t));
  - rel_inj = F + (rel − 1) × √max(F, 0.09);
  - σ_inj = σ_cal × √max(F, 0.09).
  The real residual is the noise, so nothing is added twice and **no new random number is drawn**. The rule is
  sky-limited when the source vanishes, an ASSUMPTION carried over from OGLE.
- **Thin interface** for `edge_bounds_alt`:
  - rel_inj > 0 → dm = −2.5 log10 rel_inj, σ_dm = 1.0857 σ_inj / rel_inj;
  - rel_inj ≤ 0 → σ_dm undefined, so the point is neither normal nor intermediate (bounds does not use dim points).
  - Unknown rows are absent.

## 4. Rules (no other change)
- **A:** bounds on ZTF r (v23 / v24; `vet` with EDGE_RULE 'bounds').
- **B:** A with the ATLAS points of the carrier as added edge points (`edge_bounds_alt`).
  - Run finding, span, EDGE_FRAC and the other three tests are A's.
  - No ZTF g, no #460 joint extension, no model fitting, no interpolation across gaps.
- **Secondary (reported, not gating):** B with the quoted σ, and B keeping the 9999 rows.

## 5. Denominators, success and stop (fixed now)
- **Primary:**
  - restored = W3 events on the 40 stars with A "edge_insufficient" and B "recovered".
  - **Success** if restored / 63 ≥ 0.50 **and** B fade passes / 79 fade injections on the 40 stars ≤ 0.10. Both areas
    are pooled. These are #459 / ASAS-SN plan §3, unchanged.
- **Also reported:**
  - all 154 W3 and 79 fade injections on the 40 stars (A and B outcome tables);
  - A-recovered W3 events lost in B;
  - fade insufficient → slow or → pass;
  - Wilson 95 % intervals;
  - per-star ATLAS status;
  - the breakdown without the 4 development stars.
- **Missing data:** a star without ATLAS data (failed, timed out or unfinished) keeps B = A and stays in every
  denominator.
- **Completion:**
  - **Complete** if every task has ended: data, or a definitive server error / no result.
  - **Incomplete** if any task is unfinished at the acquisition deadline. Then the primary numbers are reported as
    incomplete, with bounds:
    - lower bound: missing stars B = A;
    - upper bound: every A-insufficient event of a missing star is restored and no missing-star fade passes.
- **No second try:** no seed, threshold or rule change, and no second method after B is seen. A pass means only that an
  independent validation is worth doing. No search group is read.
- **Labels:**
  - simulated: every "dimming", "restored" or "pass" here is a simulated injection;
  - observed: ATLAS and ZTF epochs and their noise;
  - nothing is a real event.

## Deviation (09:17 UTC, before any ATLAS flux was read)
- **Why:** the shared queue runs our tasks one at a time. The first full-range task (MJD 58150–61000) had been running
  for 16 min by 09:17, so 40 full-range tasks cannot finish in the budget.
- **Change:** the 38 stars not yet submitted request only the hull of the windows that can affect B, from the ZTF
  replay alone (`hull_script.py.txt` → `hull_<npix>.json`):
  - per event with a bounds hit, [t_start − span, t_end + span] (`edge_bounds_alt` uses nothing farther away), and
    ± 60 d for the baseline.
- **Unchanged:** stars, order and method. The two submitted tasks keep the full range.
- **Note (added at review):**
  - For hull-range stars the baselines and k come from the hull, not the full range.
  - A normal point beyond one span can change the recorded U but not the class.
- **Restart:**
  - The scheduler was stopped by mistake (`pkill -f` also matched the calling shell) and restarted from its saved
    state, so nothing was resubmitted.
  - The acquisition deadline stays 11:49 UTC.

## Review fixes (10:00–10:10 UTC, after the code review, before any A / B outcome was read)
A smoke run at 09:48 checked only the calibration and the A = v24 check; no outcome was read.
- **Scheduler:**
  - a task past 45 min keeps its slot, so there are never more than 2 in flight;
  - a POST that raises stops the scheduler;
  - result GET / DELETE errors no longer crash it, and a failed download is retried up to 3 times.
  Restarted at 10:04 from its saved state; nothing was resubmitted.
- **k fallback:** a filter with < 20 usable rows uses the star's pooled k, not 1.0. A star with < 20 usable rows in
  all has no usable ATLAS (B = A).
- **Completion:** complete only if every target ends with data or "finished without result".
- **When incomplete:** success is also reported at both bounds.
  - best case: every missing-star A-insufficient W3 restored, with the observed fade passes;
  - worst case: none restored, and every missing-star A-insufficient fade passes.
- **Checks in `result_script`:**
  - the fixed denominators (154 W3, 63, 79);
  - every star with data has `atlas` true in the trial events;
  - `calib.json` is written from the trial calibration.
