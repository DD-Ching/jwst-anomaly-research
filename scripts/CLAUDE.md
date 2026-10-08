# scripts — search conventions (screen / vet / inject / limit)

Root CLAUDE.md applies. Per-search recipes and their failed-approach rules are directory-scoped skills in
`scripts/.claude/skills/` (`w3-survey`, `w12-lenscats`); load the one for the search you touch.

## Every exotic search
- Phases: prediction (`exotic_sim`) → screen → vet (cheapest ordinary test first) → injection-recovery → limit.
  Never write "discovery" or exotic language before every ordinary test has run; the strongest wording is
  "unexplained under tests X, Y, Z". Anything that survives every ordinary test goes to `/vet-candidate` and the
  owner is notified; nothing is announced outside the repo without the owner.
- A null is reported as a quantitative limit, valid only over systems where injection proved the test could find
  the signal. Inject through the *whole* chain (selection emulation, pre-screen, fit, every vetting test).
- Re-run injections after ANY pre-screen or vetting change before quoting a limit.
- Speed never removes a vetting test. A vetting test that removes injected signals must be calibrated (and its cost
  stated in the efficiency).
- Look at the contact sheet of flags and survivors yourself before concluding.
- Thresholds are ASSUMPTIONs, named in `Params`; outputs carry provenance labels.

## Parallel topology (owner decision 2026-10-08)
- Classify each stage: I/O-bound → 8–16 concurrent HTTP connections or range reads with back-off on 429/5xx —
  except services with published rate limits or etiquette (arXiv API: one request per ~3 s; CDS, MAST, Data Lab:
  batch queries instead of concurrency); a rate-limited query is retried, never recorded as a pass;
  CPU-bound → process pool sized to the cores with `OMP_NUM_THREADS=1`; never oversubscribe (load > cores stalls
  everything; two 4-process pools on 4 cores stalled both).
- Pipeline: download batch N+1 while processing batch N; write small tracked per-chunk tables under `results/`
  (D-059, D-062) so any session can resume or merge; delete raw batches after use; in cloud sessions never store a
  whole archive tar.
- Log per batch: items/s, CPU %, MB/s. If CPU < 70 % or the network idles, fix the pipeline before scaling.
- Pool workers must receive module-level settings through an initializer (spawned workers on Windows inherit none).
- Wait on a log line or an output file, never `pgrep -f` / `pkill -f` (they match the waiting shell itself).
