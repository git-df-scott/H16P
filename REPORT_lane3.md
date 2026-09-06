# REPORT lane 3 — exact certifier for normal-size limit cycles

Branch `fable/lane3-certify`.  Status: core certifier built and working; first
field certified end to end.  Independent `verify.py` and the full validation
table are in progress.

## What was built

`lane3/` holds a Poincare-Bendixson certifier that turns numerically located
limit cycles of a planar quadratic field into a machine-checkable certificate.
No interval ODE integrator is used anywhere.

- `exactgeom.py` — rational predicates: star-polygon certification, exact sign
  of the edge quadratic `X . n` on `[0,1]`, radial containment, box/polygon
  position, interval enclosure of a quadratic on a rational box.  Also an
  integer reformulation of the edge test (same decision, faster) used only
  inside the proposal loop.
- `fields.py` — exact rational `vec12` for Cherkas rows 1-8, the KKL control,
  Perko P3 (rotated by an exactly rational rotation) and a constructed control.
- `propose.py` — floating point only: polar return map, cycle location,
  multipliers, and the two-turn phase blend that proposes the polygons.
- `certify.py` — the driver and every exact check; writes
  `certificates/CERTIFICATE_<name>.json`.
- `verify.py` — independent re-checker (in progress).

## The construction

For each cycle, two closed polygons with rational vertices are built, one just
inside and one just outside, and the flow is proved to cross every edge towards
the cycle.  The proposal is the key idea: integrate one orbit for TWO returns to
a ray from the focus and blend the two turns at matched phase,

    C(tau) = (1 - tau) Z(T1 + tau (T2 - T1)) + tau Z(tau T1),  tau in [0,1],

which is closed exactly (both ends are `Z(T1)`) and drifts away from the cycle
at the rate at which the flow contracts towards it.  The flow is tangent to each
turn, so relative to `C` it has a normal component of one sign, of size about
`delta * |1 - multiplier| / period`.  Blending at matched PHASE rather than at
matched polar angle is essential: with angle matching the gap vector between the
turns is radial and its normal component collapses wherever the ray is oblique
to the orbit, which killed the first implementation.

Exact steps in the certificate: (T) one strict sign of the edge quadratic on
each closed edge, decided from its endpoint values and its vertex; (S) each
polygon simple, counter-clockwise, star-shaped about a rational centre;
(N) inner strictly inside outer; (E) no equilibrium in the closed annulus,
via resultants plus rational root isolation plus exact box placement;
(D) annuli pairwise disjoint.  Signs `(inner +1, outer -1)` give a positively
invariant annulus, `(inner -1, outer +1)` a negatively invariant one; either
yields a periodic orbit by Poincare-Bendixson.  A separate exact check that every
finite equilibrium has non-zero divergence rules out a centre and upgrades the
periodic orbits to limit cycles.

## Validation so far

| field | cycles found (r on the ray from the focus) | multipliers | certified | vertices (in+out) | wall |
|---|---|---|---|---|---|
| Cherkas row 1 | 0.2809, 1.0070, 3.0193 | 1.0035, 0.9901, 1.0834 | 3 limit cycles | 34k / 66k / 51k | 82 s |

Cherkas row 1 crossings of the paper's section `y = -1` are `x = 1.2809,
2.0070, 4.0193` against the paper's `1.26, 1.98, 3.95`; the paper's table is
quoted to two decimals and the third differs by 1.7%.  The certified count (3)
agrees with the paper's exact Dulac-function count.

## Open / next

- `verify.py`, independent re-check of the emitted certificates.
- Cherkas rows 2-8, KKL (3+1), Perko P3, the constructed one-cycle control.
- Cross-nest disjointness (bounding-box separation) for the (3,1) fields.
- Vertex counts scale like `1/(delta * |1 - multiplier|)`; the KKL innermost
  cycles with multipliers within 1e-3 of 1 are expected to be the hard case.

## Limits (documented, not worked around)

Cycles whose multiplier is within about 1e-6 of 1, and infinitesimal cycles such
as the Songling ones certified by Galias-Tucker at radii down to 1e-202, are out
of reach of this method: the required edge count grows like the inverse of the
transversality margin.  This certifier is for NORMAL-SIZE cycles.
