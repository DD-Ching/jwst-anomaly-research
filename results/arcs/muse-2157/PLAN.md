# `arcs` anti row MACS0308 2157: redshift from the public MUSE cube (plan, written before any data access)

- **Question:** is the `robust` anti row MACS0308 2157 a cluster member (z ≈ 0.356) or a background source
  (stored photo-z 0.756)? Its anti classification (D-024 addenda) assumes a background source.
- **Data:** ESO MUSE cube `ADP.2020-01-04T23:56:39.902` (0104.A-0801(A)), streamed by HTTP byte ranges, aperture
  spaxels only; no full cube stored (Cloud disk decision).
- **Method:** circular aperture (ASSUMPTION r = 0.6″) at the stored 2157 position; local annulus sky; search for
  [O II] 3727, Ca H+K, Hβ, [O III], Hα at both redshifts.
- **Gate M (validation, before reading 2157):** the same extraction on a catalogued cluster member with a
  published MUSE redshift in this cube (or a bright star) must recover its redshift within 0.002.
- **Frozen reading:** `member` if ≥ 2 features agree with |z − 0.356| < 0.02; `background` if ≥ 2 features
  agree with z > 0.40; one feature only → `single_line`; none → `undetermined` (with the continuum S/N stated).
- **Would change the plan:** `member` makes the anti row an ordinary foreground galaxy (anti label withdrawn);
  `background` keeps it for the lens-model tests.

## Amendment before any cube data was read (2026-11-08T20:25Z)
- The ESO portal ignores HTTP Range (416 with unknown length) and streams the whole file, so the DATA extension
  (1.55 GB) is streamed once and only two 80 × 80-spaxel boxes and a white-light image are kept in memory; nothing
  is stored. Noise comes from robust scatter of random apertures in the box (STAT extension not read).
- No object in this field has a secure published spec-z in NED (2026-11-08 TAP query: only the cluster z 0.356 and
  photo-z rows). Gate M is therefore the central galaxy (NED WISEA J030855.93+264537.9) recovering the cluster
  redshift: ≥ 2 features (|S/N| ≥ 5) and |z − 0.356| < 0.01 (ASSUMPTION: allows ~±2200 km/s for the central galaxy's
  velocity offset and the cluster z precision).
- The reading is computed by `classify` in `scripts/arcs_muse_2157.py` (z scan 0–1.6, step 2e-4); `foreground`
  added for ≥ 2 features at z < 0.40 outside the member window.
