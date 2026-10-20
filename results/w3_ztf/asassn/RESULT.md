# W3 ZTF × ASAS-SN — light-curve service down (3 probes, 22:11–22:31 UTC); 40 stars selected and baselined, no flux read; trial not run

- **Pre-registration:** the plan was committed at 22:10:58 (5ac8dc2). Before it, about 22:05–22:10, connection tests
  were made by hand (not stored): the service message, the block-server list (`[]`), the catalogue list and one
  light-curve request, which failed at once because no block server existed. No flux was received.
- **Deviations:**
  - The PLAN time stamp was corrected after probe 1 (see PLAN).
  - The 60 s / 120 s timeouts apply to the raw requests only; the client's own calls use its defaults.
  - The ZTF-only baseline (section 3) was added after the plan (2e4ceeb); it uses no new data.
- **Labels:** service responses are `observed`; the selection uses the ASAS-SN catalogue (`observed`) and v22
  samples (`derived`); the ZTF baseline counts are `derived` from simulated events.
- **Scope:** a development trial on already-used calibration stars. No search group was read. ZTF psfcat
  restoration stays closed; 8124, 8126, 1624 and 1625 are untouched.

## 1. Entry (`probe.json`)
- **Access:** the official client (`skypatrol` 0.6.21) uses plain HTTP on port 9006. In this cloud session it works
  only through a CONNECT tunnel (`http_shim.py.txt`, an environment shim, not a method change).
- **Catalogue cones work:** 1.6–1.7 s per probe.
- **The light-curve service is unavailable.** Probes ran at 22:11:15, 22:21:26 and 22:31:27:
  - each time `get_block_servers` returned `[]`, so no light curve can be requested;
  - each time the service message said "The light curve service will be down for several days while we run
    mainanence on our cluster" (dated 13 October 2026; quoted as served).
- **As planned:** no light curve was read, and steps 3–5 (flux handling, paired injections) did not run. No other
  survey or path was opened; the ASAS-SN web forms need an e-mail submission.

## 2. Selection (`selection.json`, catalogue only; no flux)
- **Per area,** the 60 brightest calibration groups were examined:
  - **4991:** 58 matched an ASAS-SN `stellar_main` source within 3″, and 55 pass the contamination cut;
  - **2344:** 59 matched and 58 pass.
  - The contamination cut is `rp_00_1` ≥ 16″, i.e. neighbours inside 16″ add < 10 % of the G flux.
- **Selected:** the first 20 passing per area, 40 in total.
  - ZTF refmag_r is 14.51–14.74 in 4991 and 14.51–14.76 in 2344.
  - Gaia G is 14.54–14.82 and 14.34–14.88.
- Each star's ASAS-SN id, separation, Gaia id and magnitude, `rp_00_1` / `rp_01` / `rp_10` and `nstat` are stored.
  `nstat` is 0 for these stars, so the catalogue holds no variability statistic for them. All selected stars have
  `rp_00_1` = 99.9, read as the catalogue's "no neighbour found" value; it passes the cut.
- **Limit:** these are the brightest calibration stars of two areas, chosen for ASAS-SN's precision. Any later
  result applies to stars of about r 14.5–14.8 with clean 16″ apertures, not to the fainter bulk of the W3 sample
  (calibration median r 16.42 / 16.46; `../v22/<npix>/sample_info.json`).

## 3. ZTF-only baseline for the 40 stars (`ztf_baseline.json`, from `../v24/step12_events_*.json.gz`; no new data)
- **Events counted:** v23 events whose carrier is a selected star, flagged events only (events that were not
  flagged carry no carrier in the v24 file).
- **W3:**
  - 4991: 66 flagged events. v13 recovers all 66 flagged events; bounds recovers 27 and calls 39 insufficient.
  - 2344: 47 flagged events. v13 recovers 45 (2 slow); bounds recovers 21 and calls 24 insufficient.
- **Fades:**
  - v13 passes 7 of 32 (4991) and 5 of 41 (2344);
  - bounds passes 0, and all of those 12 are insufficient.
- **Undecided sides** (of bounds-insufficient events only; the fade row was restricted after review):

  | | A (no further covering exposure) | B (quality-flagged) | C1 (DR24 gap) |
  |---|---|---|---|
  | W3, 4991 / 2344 | 27 / 25 | 20 / 9 | 4 / 3 |
  | fade, 4991 / 2344 | 7 / 3 | 2 / 3 | 2 / 1 |

- **What a later ASAS-SN read must show:** whether ASAS-SN epochs fall in these A and B intervals and give a
  normal or significant intermediate point within 0.25 × span. These are 63 insufficient W3 events and 12
  insufficient fades (sums of the stored per-area counts).

## 4. Cost and time
- **Data:** under 1 MB of catalogue responses (not measured exactly; no light curve).
- **Time:** about 26 min from the start (about 22:05, not stored) to the last probe.

## 5. Answers available now
- **What ASAS-SN adds:** not measured. The light-curve service was down for the whole entry window.
- **W3 regained / fades removed:** not measured. There are 0 paired injections.
- **Which stars:** only the 40 selected bright, uncontaminated stars would be testable (section 2), with their
  ZTF-only targets fixed (section 3).
- **Next step:** re-run the same plan unchanged once `get_block_servers` is non-empty. Use the same 40 stars and the
  same ZTF baseline, so nothing has to be re-selected. Do not switch to another survey in its place without a new
  decision. Until then this route is **paused, not closed**: no evidence for or against ASAS-SN's value was
  obtained.
