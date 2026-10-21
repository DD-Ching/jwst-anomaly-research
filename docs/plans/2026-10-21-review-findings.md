# Parallel round 2026-10-21: /code-review findings, fixes and their effect on committed results

Companion to `docs/plans/2026-10-21-parallel-round.md`. Each row gives the fix, the test, whether the affected path
ran in the committed runs (with the evidence), and the effect on committed numbers. "The result file did not change"
is not used as proof on its own: the evidence is the run records (logs, JSON row statuses, timestamps, sha256) or a
targeted recomputation.

Run records outside git (session scratchpad, not committed): the line B download log (72 lines, URL / bytes / sha256),
the line C per-stage logs, the line A per-fit state and logs. Their content that matters is quoted below.

## Line B (`scripts/s1_chime.py`, `scripts/s1_chime_lib.py`, `tests/test_s1_chime.py`)
| Round | Finding | Fix (commit) | Test | Path ran in committed run? (evidence) | Effect on committed numbers (evidence) |
|---|---|---|---|---|---|
| 1 | `chi2_union` keeps the sub-bin δ with the best p even when its gains fail | best p among δ with detected gains (1531dd9) | failure case d = 0.2 / −0.3 | **Yes**, every pair and injection | **Verified by re-run** from stored profiles: no pass changed in 4,680 injections or 136 real pairs; frozen verdict identical; post hoc I3 +3 passes (≤ 0.2 pp). v1 outputs kept (`*_v1.*`) |
| 1 | `run --plots` ignores its directory | always writes into `--plots` (1531dd9) | — | Yes | Regenerated contact sheet byte-identical |
| 1 | downloads run ahead of reduction | bounded window of 4 (1531dd9) | fake downloader | Yes (old code) | Disk only (sample 374 MB); data integrity below |
| 2 | exception inside the fetch loop leaves files on disk | cancel, wait, delete, re-raise (3d676bb) | reducer that raises | **Yes, once**: fetch pass 1 stopped at 00:36 UTC after 17 transfers, on a reduction error (FRB20200314A has no calibration group; later code marks it ineligible) | None. The profiles behind every committed number come from the complete pass at 00:40 UTC (30 transfers, 30 distinct files); every repeated transfer of a file has the same sha256 (30 files, 70 transfers, 0 files with two hashes) and all 30 hashes equal the pins in `data/manifests/s1_chime_dynspec.ecsv`. Leftover raw files were deleted by hand (scratch data root now 20 kB) |
| 2 | partial file kept when curl fails all attempts | delete before returning (44d845b) | fake curl, 3 failures | Not shown by the log (failures are not logged) | None: a failed transfer is never reduced, in old or new code |
| 3 | head future invisible to cleanup during `.result()` | per-thread result registry (4912f2b) | SIGINT during wait (POSIX only, c3b1373) | No interrupt in the committed run | None (cleanup path only) |
| 3 | cleanup deletes files the run never touched | delete only touched paths (4912f2b) | reused + unrelated file survive | No (no `--keep-raw` files existed) | None |
| 3 | cleanup waits on retrying downloads | cancel event before retries / sleeps (4912f2b) | back-off cancelled < 3 s; portable child (c3b1373) | No | None (timing only) |
| — | Windows CI: `os.kill(SIGINT)` terminates the process; no `sleep` binary | skip on Windows; `sys.executable` child (c3b1373) | — | — | CI only |

## Line C (`scripts/w2_hst_blend*.py`, `tests/test_w2_hst_blend.py`)
| Round | Finding | Fix (commit) | Test | Path ran in committed run? (evidence) | Effect on committed numbers (evidence) |
|---|---|---|---|---|---|
| 1 | depth sums overlapping CAOM level-3 products | count each visit once (d9dab07) | `visit_exposure` | **Yes**, all 23 systems over-counted 1.9–7× | **Changed and re-run**: targets 10 → 7 (J0203+1612, J2257+2349, J0146-1133 out), 3 controls swapped; gate 9/10 = 0.9 unchanged; negative control S/N 13.0 unchanged (same stamp); control injection limit median 22.5 → 22.0. v1 kept (`*_v1`, `selection_diff.ecsv`) |
| 1 | `measure --only` drops other rows and rewrites the gate | merge rows, recompute on the full set (d9dab07) | merge test | No: the worker's command record has no `measure … --only` call (controls measured once with `--no-inject`, then three times with injections; targets twice, all full-set) | None |
| 1 | post hoc false-positive scaling mixes systems | same system for S/N, flux and error (d9dab07) | reference choice | Yes | None: one ePSF system, so reference 13.0 in v1 and v2 |
| 1 | `cmd_summary` crashes on an empty role | guard (d9dab07) | empty role | No (targets measured) | None |
| 1 | a failed download recorded as "no data" | distinct "download failed" (d9dab07) | failed vs no-data | No: every hapcut request in the v1, v2 and v3 download logs ends "[Done]"; no retry, error or traceback line; the one empty answer is DESJ0407-5006 with astroquery's explicit "Missing HAP files … Cutout not performed" | None (DESJ0407-5006 is a real no-data answer) |
| 1 | TAP COUNT(*) check and wide hapcut without retries | retries and back-off (d9dab07) | — | No (no error in the class and PSF logs) | None |
| 2 | `write_psf` / `fp_reference` crash on all-NaN S/N | guards, null in JSON (a1a4694) | all-NaN | No (finite; max 13.0) | None |
| 2 | `cmd_verdict` not idempotent | status always re-derived from the chain output (a1a4694) | pass and fail | Yes (run in d9dab07 and a1a4694) | None: all target statuses identical between the two runs |
| 2 | `cmd_posthoc` silently drops targets after `verdict` | controls-only on purpose, stated in provenance (a1a4694) | — | **Yes** (d9dab07 posthoc ran after verdict) | `posthoc.json` holds controls only, as the owner's rule requires; the v1 target diagnostics stay in `posthoc_v1.json`, withdrawn as readings |
| 2 | stale "download failed" never replaced; status file written only at the end | per-system status, overwrite on success (a1a4694) | — | No (no failure or interruption, logs above) | None |
| 2 | `cmd_psf --reuse` without an input check | inputs must be equal, else recompute (a1a4694) | mismatch refused | **Yes** (v2 reused 4 v1 rows) | **Checked**: all 5 reused rows have identical inputs; 3 labels normalised to "not detected"; no number changed |
| 2 | `--no-inject` leaves stale injection rows | drop them (a1a4694) | — | Only in the first v1 control pass, before any injection file existed | None |
| 2 | hard-coded "not interpreted" labels | read from `verdict.json` (a1a4694) | — | Yes | None under the failed verdict |
| 2 | `_CLASSIC` also matches individual exposures | explicit per-visit rule (a1a4694) | rule test | Yes | **Recomputed**: no depth changed in any of the 82 systems; `selection.ecsv` byte-identical (the 10 mixed visits pair a 674 s association with a separate 5 s exposure, and each has a single-visit mosaic) |
| 3 | partial cutout directory used by `run_system` / `psf_one` | require status "ok" (ed78c98) | partial folder refused | No (every request "[Done]", logs above) | None |
| 3 | NaN / inf into `allow_nan=False` JSON; reuse refusal bookkeeping; missing reuse column; NO_HAP re-queried | guards and bookkeeping (ed78c98) | six tests | No (committed JSON written without error; `reuse_refused_input_mismatch` is `[]`) | None |

## Line A (`scripts/x4_s1g.py`, `tests/test_x4_s1g.py`)
Fixes in 7dada57; run-record checks and small recomputations from stored inputs and best parameters in 898212e
(`results/x4/s1g_review_checks.json`). No fit was re-run.

| Round | Finding | Fix | Test | Path ran in committed run? (evidence) | Effect on committed numbers (evidence) |
|---|---|---|---|---|---|
| 3 | thread limits set after numpy import | set before import | — | **Exercised.** Every launch exported OMP/MKL/OPENBLAS = 1 in the shell before Python started (session record); lalsuite has no OpenMP (`ldd`: no libgomp; non-threaded FFTW; GSL with single-threaded CBLAS); s per generation 2.69–3.08 (validation) and 2.33–2.35 (R1) against the 2.64–2.77 benchmark; load 3.33 with 3 fits | **None (recomputed):** all 9 stored best fits re-scored with 1 thread and with default threads (OpenBLAS 4): max \|ΔFF\| 0.0 between modes and against the stored values, at pad 64 and 32 |
| 3 | box round trip not enforced for the validation draws | enforced at 1e-9 (post hoc tolerance; the plan named none) | committed truths pass | **Yes** | **None (recomputed):** `round_trip_max_diff` 0.0 for VA1 and VA2; objective and pad-64 FF change 0.0 after `from_cart(to_cart(x))` for VA1, VA2, V5, V9; validation set unchanged |
| 3 | feasibility outcome unreachable | plan order restored; flag exposed | three tests | Not entered: FF_conservative 0.947638 < 0.98790 | None: outcome 3 and its causes unchanged; `s1g_summary.json` only gains `feasibility_fail_verified: false` |
| 3 | refinement hook saves a rejected restart | accepted state only | yes | Not used: states are read only for fits stopped at the deadline (none); for all 9 fits the old hook's last state equals the final row | None |
| 3 | stale state files | run ID, cleared, ignored | yes | No: the crashed first launch (00:43:49, `TypeError` in the DE callback) never wrote state; the directory was empty at the 00:49 relaunch; the 9 state files present are the relaunch's (01:38–03:03) | None (0 rows use a partial state) |
| 3 | a raising fit kills the round | "error" rows | yes | No traceback or unfinished/error line in the relaunch log; 9 / 9 rows "finished" | None |
| 3 | `tool_final` without the completeness guard | guard | stored-summary diff | 6 / 6 ran | None ("pass" before and after) |
| 3 | R1 resume overwrites finished rows | keep finished rows | yes | R1 ran once, 02:31:22–03:03:22, no resume | None |
| 3 | test mutates the shared module | monkeypatch | itself | Tests ran in a separate process | None |
| 3 | specs resolved per job | resolved once in the parent | — | Spec sources unchanged during the run (`s1g_injections.json` since 00:42:54, afcdf9c) | None (recomputed: targets rebuilt from today's sources re-score every stored fit with \|ΔFF\| 0.0) |

Caveat: the shell launch lines and the empty state-directory listings come from the session transcript; the state
files and logs are in scratch, and the values used above are copied into `s1g_review_checks.json`.

## Round 4: bounded verification of the round-3 fixes (coordinator)
Scope: only the diffs of 7dada57, 4912f2b, c3b1373 and ed78c98. Result: 17 of 18 round-3 findings fixed as listed;
the refinement-state fix was partial, and the fixes introduced the regressions below. All are on paths no committed run
used (no fit, reduction or measurement was re-run after round 3), so none can change a committed number. Fixed by the
coordinator in one commit (code and tests only):

| Finding | Fix | Test | Effect on committed numbers |
|---|---|---|---|
| A: `refine_state` compares with the exact-lens start value, while `refine` may compare with the DE grid-lens value in round 0, so "accepted" can disagree with `refine` | the state is now named and documented as the best exact-lens point seen (`improves_exact`), which is what a partial row should re-score; `refine` unchanged | renamed test | None: states are read only for fits stopped at the deadline (none) |
| A: thread limits used `setdefault`, so an inherited value > 1 was kept | always set to 1 before numpy loads | — | None: every committed launch exported 1 (above) |
| A: a fresh `run` resumed from the committed `s1g_r1.json` | only the `r1` command resumes, and only if the stored budget equals the current one | fresh path ignores the file; resume refuses a budget mismatch | None: R1 ran once, before this code existed |
| A: resumed record kept the first run's budget | covered by the budget-equality check | same | None |
| A: budget recorded in two steps | one dict literal | — | None |
| C: `psf_one` deletes a complete wide cutout that has no `wide_*` record | not changed: `psf_one` deletes every wide cutout after use (cloud-disk rule), so a later run downloads it in any case | — | None |
| C: whole-file status writes can drop another stage's key | `set_status` re-reads the file, replaces one key and writes atomically | key-merge test | None (stages ran sequentially) |
| C: download-with-status sequence duplicated in `download` and `psf` | one shared `download_set` | existing download / partial-dir tests | None |
| C: `run_system` re-reads the status file per system | not changed (a cheap read; stages run sequentially) | — | None |

A last verification of the round-4 commit (188102b) raised two points: concurrent `psf_one` writes (not applicable:
`cmd_psf` calls `psf_one` sequentially in one loop) and two blank lines left by removed imports (removed, 5616671). No
further review round was opened (owner instruction: verify the fixes; do not widen the review).

## Review process
- Round 1: B and C first versions (10 findings). Round 2: their fix commits (10). Round 3: the final PR diff (10 in
  line A) and the B / C fix commits (10). Round 4: a bounded verification of the round-3 fixes only (above).
- The coordinator extended B's and C's code-fix window from 03:30 to 04:15 UTC (code only: no re-runs, no downloads).
  This changes the coordinator's allocation, not the owner's 6 h budget. Line A's computation stop (05:15) was not
  used past 03:03.
