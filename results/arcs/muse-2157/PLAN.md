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
