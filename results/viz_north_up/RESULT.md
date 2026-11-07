# viz.north_up and shared contact-sheet scaling: validation (2026-11-07)

**Tool:** `viz.north_up(data, wcs)` resamples a cutout about its central pixel so north is up and east
left, from the local sky Jacobian of the cutout's own WCS (rotation, parity, skew and mosaic meridian
convergence included); NaN-aware bilinear by default. `contact_sheet(..., north_up=True, scale="panel" |
"band" | "all")` and `cutout_png(..., north_up=True)` use it. Defaults keep the stored orientation and
per-panel scale. The contact-sheet footnote now states the orientation shown, so default PNG text changes.

## Known cases first (synthetic, `tests/test_viz.py`; observed = none, simulated)
- Rotations 0°, 37°, −121° (east-left parity), 90°, 210° (east-right parity) and a mosaic-like cutout 12′ from the
  tangent point at dec −73°: a neighbour placed 0.27″ due north lands 9 px straight above the target, one placed
  due east 9 px to the left (within 1 px, the nearest-input-pixel placement), and the target stays central.
- Identity WCS reproduces the input (≤ 1e-6, TAN-plane curvature); 90° equals `np.rot90`; a one-pixel hot pixel
  survives a 210° flip (nearest neighbour would drop it, hence bilinear by default; failed approach below); an
  isolated NaN stays NaN (grown by about a pixel; any output pixel with < 0.99 valid weight is NaN, an
  ASSUMPTION); a galactic-frame WCS is turned to ICRS north.
- `scale="all"` gives one shared interval; `"band"` one per band; a cutout without a celestial WCS is drawn as
  stored and labeled "no WCS", never silently mis-oriented.

## Real data (observed pixels, derived positions)
`scripts/viz_north_up_check.py` on SMACS 2736 o001 F200W (S3 i2d and pipeline `cat.ecsv`; cutout PC rotation
35.39°, stored in `check.json`): 11 unordered catalogued pairs (target 20–23 mag, neighbour < 24 mag, 1.08–2.01″
apart; mutual pairs tested once; selection limits are ASSUMPTIONs listed in `check.json`). In the north-up cutouts
the brightest pixel within ±5 px of the neighbour position predicted from its east/north sky offset is not on the
window edge and lies ≤ 4 px from the prediction for 10 / 11 pairs (matched max 3.47 px; median over all 1.75 px).
The miss, c3046, has an elongated galaxy neighbour whose brightest pixel lies on the window edge, 4.71 px from the
catalogue centroid. Control: in the stored orientation the same window holds ≤ 0.199 × the north-up peak. The
edge and 4 px match rule was set after a first run with a ±3 px window (post hoc); that run counted each mutual
pair twice and accepted window-edge maxima (both fixed after `/code-review`). Contact sheets looked at: NIRCam
six-point spikes appear rotated by the same angle in every panel, and the neighbour layout matches the catalogue.

## Failed approach
- Showing pixels with ≥ 0.5 valid bilinear weight: an isolated interior NaN vanished after a 37° rotation
  (found by `/code-review`), so a masked bad pixel would look like sky.
- Nearest-neighbour (`order=0`) resampling: at 210° with a parity flip it dropped an isolated one-pixel
  value entirely (rotated grids skip input pixels), which could hide a hot pixel or cosmic ray during vetting.
- `mode="constant"`: the identity transform turned the last row and column NaN (matrix 1 − 1e-9 pushes edge
  samples past n − 1); `mode="grid-constant"` treats each pixel as ±0.5 px and keeps them.

## Limits
Display only: resampling smooths at the pixel level (bilinear), so measure on the stored cutout. The
0.99 coverage threshold is an ASSUMPTION. Multi-band panels are each rotated on their
own WCS; they are not put on one common pixel grid (that is `reproject`'s job; D-005 addendum).
