# W3 ZTF × ASAS-SN — can ASAS-SN Sky Patrol supply the missing edge information? (one development trial)

Plan fixed 2026-10-20, 22:10 UTC (commit 5ac8dc2 at 22:10:58), before any ASAS-SN flux is read. Owner authorization after #414: 3 h, ≤ 2 GB, one
small measured trial. Uses the official client `skypatrol` (pyasassn 0.6.21) and the existing W3 code. ZTF psfcat
restoration stays closed (v24); 8124, 8126, 1624, 1625 untouched; no search group read.

## 0. Prior state (checked 22:05)
- No ASAS-SN access or evaluation exists in the repository (grep: only literature mentions). No open PR or branch
  claims this.
- **Access path (this cloud session only):** the client speaks plain HTTP to port 9006, while this sandbox egresses
  only through a CONNECT proxy. `http_shim.py.txt` forwards each request through a CONNECT tunnel. It is an
  environment shim, not part of the method.

## 1. Entry check (≤ 30 min from 22:05; `probe_script.py.txt` → `probe.json`)
- Record the service message, the light-curve block servers and one catalogue cone, with time stamps.
- The light-curve service counts as available only if `get_block_servers` is non-empty **and** one test light curve
  downloads within 120 s.
- **Retries:** up to 3 probes, at least 10 min apart, by 22:40 UTC. Each request has a 60 s connect / 120 s read
  timeout.
- **If unavailable:** save the evidence and stop the flux work. No second survey, no other ASAS-SN path (the
  web forms need an e-mail submission).

## 2. Source selection (catalogue only; `select_script.py.txt` → `selection.json`)
- **Candidates:** the 500 calibration groups of v22 4991 and of v22 2344, in ascending `refmag_r`; at most 60 are
  examined per area.
- **Match:** the nearest ASAS-SN `stellar_main` source within 3″ of the group position (one cone of 5″ per group).
- **Contamination (fixed now):** ASAS-SN photometry uses 16″ apertures (2 pixels of 8″). Keep a star only if the
  catalogue's `rp_00_1` ≥ 16″, i.e. neighbours within 16″ add < 10 % of the star's G flux (ATLAS Refcat2 definition,
  carried in `stellar_main`). A missing value fails.
- **Take the first 20 passing per area** (≤ 40 in total).
- **Recorded per star:** the ZTF key and refmag_r; the ASAS-SN id, separation, Gaia G, `rp_00_1`, `rp_01`, `rp_10`
  and `nstat`.
- Selection reads no light curve and no injection result.

## 3. Flux handling and trial (only if section 1 passes; fixed before any flux)
- **Data:** ASAS-SN light curves, both g and V, with all epochs kept: low-S/N and negative fluxes with their
  errors. Upper-limit magnitudes are never used as measurements. A missing epoch is never zero. Quality: drop only
  epochs with `quality` ≠ "G" (the client's flag).
- **Times:** HJD − 2400000.5, compared with ZTF HMJD.
- **Baseline:** per camera and band, the median flux over the whole light curve. The relative flux is f / median;
  there is no seasonal normalisation or detrending.
- **Joint model:** the same physical event, with dm(t) from `w3_dmag` or `control_model`, is applied to each
  survey's own times. Each instrument keeps its own baseline, and the blend fraction is the cell's blend.
- **Noise of added epochs:** each epoch's own flux error, from a separate fixed seed (5000 + cell, 5100 + control).
- **ZTF inputs:** unchanged, with the v23 seeds; ZTF and its noise draws are identical to v23.
- **Edge information:** a side is decided (bounds) if a normal ASAS-SN point (relative flux ≥ 1 − 0.24, i.e.
  dm < 0.3, and S/N ≥ 3) lies within EDGE_FRAC × span of the run edge. A significant intermediate ASAS-SN point
  (0.3 + 2σ ≤ dm < 1) at more than EDGE_FRAC × span is a slow edge. Everything else is as v13 / bounds.
- **Success (development):** in the subset of injections on the selected stars, bounds with ASAS-SN restores
  ≥ 50 % of the ZTF-only bounds "insufficient" W3 events, while fades passing stay ≤ 10 %.
- **Stop:** 01:00 UTC or 2 GB.
