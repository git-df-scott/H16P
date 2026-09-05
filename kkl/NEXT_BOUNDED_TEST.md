# Next KKL test: a specific c>1 field, not a coefficient search

Status: **four initial remote controls completed; no U bracket found;
conditional origin work not run**. This bounded test follows the 202-call
checkpoint. The field is

\[
 \dot x=y+x^2+xy,\qquad
 \dot y=-10x^2+\frac{11}{5}xy+\frac{1001}{1000}y^2-\frac{196}{5}x.
\]

Here beta=0 and `K=32039/6250>1/64`. Exact rational interval checks give
one simple real remote equilibrium, negative trace, positive determinant,
and a negative focus discriminant. These are equilibrium certificates,
not periodic-orbit certificates.

The analytic reason for testing this particular field is the changed
multiplier-density geometry across c=1. Its quartic N has two positive
zeros in u=1+x, both greater than one, with a bounded negative interval
between them. Thus the necessary stability condition has a different
amplitude structure from the first explored c<1 component. This does not
prove a stationary return branch or any limit cycle exists. See
[the exact derivation](notes_other_strata.md).

The new test makes no persistence claim across either infinity boundary.
The compactified directions and their types must first be checked at the
new field. No interpolation from the previous large remote cycle is used.

## Bounded protocol

1. Record the exact finite and infinity gates. Keep the existing section,
   coordinate, CPU and total evaluation limits.
2. Test the remote downward section first, because an origin pair without
   a remote U cycle is insufficient. Allow at most eight initial remote
   return/derivative calls at this fixed field, with radii chosen within
   the prescribed range. Keep every nonreturn unresolved; never join a
   sign bracket across a change of itinerary or a failed return.
3. If an admissible remote U root is isolated numerically, use at most
   four initial origin controls, then check its stationary points if those
   controls supply a seed, using the quartic's
   amplitude condition as a guide. The quartic roots are x-amplitude
   restrictions, not section roots or return roots.
4. Cap this first test at 64 additional charged evaluations. Proceed to
   continuation only from an actually located return or stationary branch.
   Do not add a parameter grid if the fixed-field test supplies no seed.

A useful positive result is a remote U root plus an additional origin
stationary branch that can lead to the missing S/U pair. The actual
precursor gate remains three ordered origin S/U/S roots and the remote U
root at one common field. A five-cycle field requires the subsequent
negative-beta step and five rigorous return certificates.

If no seed appears within this test, report that result for this field and
these section intervals. Neither the c>1 stratum nor the full KKL box is
thereby excluded. The untouched budget remains available for a separately
justified, bounded continuation task.

## Result of the initial remote controls

Calls 203–206 all completed a numerical full return at the same field:

| Initial r | D(r) | R_r(r) |
|---|---:|---:|
| -8 | +0.1934757492 | 0.9382092882 |
| -512 | +259.8321636 | 0.3479611859 |
| -32768 | +29886.64503 | 0.0416319720 |
| -1048576 | +1037240.44514 | 0.0035049036 |

All four displacements point inward on this section. There is no U-root
sign bracket to refine, so the initial test stopped after four calls and
the conditional origin block did not run. The determinant and independent
speed/divergence derivatives agreed to less than 1.6e-12 in these controls.
The returns stayed left of x=-1 and below the coordinate guard.

These are ordinary numerical controls. They do not exclude roots in
unsampled intervals or prove absence of a remote cycle at this field.
No interval return certificate or new parameter continuation was made.
Raw results are in [c_gt_1_remote_controls.json](data/c_gt_1_remote_controls.json)
and the append-only ledger. Total charged evaluations are now 206,
leaving 3890 of the inherited allowance.
