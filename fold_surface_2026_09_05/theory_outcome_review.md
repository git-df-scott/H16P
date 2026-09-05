# Independent saved-outcome audit

The saved evidence supports a continued **numerical fold pair**, with an
inner stable and outer unstable cycle on its sampled pair side. It does
not currently support three origin cycles at any one sampled field. It
also does not prove global absence of a third cycle or complete the fold
component. No saved result examined here requires rejecting the existence
of the fold branch.

This bounded audit used no ODE evaluations. `theory_outcome_audit.py`
reads saved results, checks their internal sign evidence, and performs
exact rational equilibrium calculations. Its JSON records SHA256 hashes
of the snapshots. The ongoing quad continuation was read through its
sixth accepted point, r approximately 9.9257e13, c approximately
1.59161988024, K approximately 7.03091589407. Later points are outside this
snapshot and must not be attributed to this audit.

## Fold and pair evidence

| Saved sequence | Accepted folds / events | Saved pair brackets |
| --- | ---: | --- |
| increasing K | 18 / 18 | Two at 17 fields; one at K=6 |
| decreasing K | 4 / 4 | Two at all four fields |
| arclength | 17 / 27 | Two at all accepted fields |
| full angular | 1 / 2 | Two at its accepted field |
| two-half long double | 8 / 8 | Two at all eight fields |
| two-half binary128 snapshot | 6 / 6 | Two at all six fields |

Every saved bracket checked has ordered endpoints and opposite recorded
displacement signs. The raw adjacent profile sign changes reproduce the
saved bracket counts. No checked profile contains three sign brackets.
For assigned stability labels, the stored derivative/matching-residual
data agree with inner stable and outer unstable. The original full-angular
profile does not assign root stability labels. Root-refinement results
are approximate: a reported half-map `multiplier_at_match` is a cycle
multiplier only after matching, and a small absolute F alone can be weak
evidence at very large r because its sensitivity also tends to zero.

All accepted numerical fold curvatures checked are positive: L_zz for
full returns and G_z for half-map matching. Thus the earlier minimum-fold
assumption in the pair-side shift did not select the wrong side at these
stored folds. The owner has subsequently added the general curvature
sign and unresolved-refinement handling identified in the source review.

The K=6 one-root profile does **not** refute that fold. Its saved parameter
shift is +0.012, from c=1.53775968728 to c=1.54975968728, rather than the
much smaller curvature-scaled shift used subsequently. It is a coarse
excursion from the fold. Its sampled outer displacement remains negative,
and only the stable root is bracketed. Later arclength and half-map
profiles use smaller shifts and show the two-root side again. The missing
second bracket at this coarse field cannot be relabeled as a second
cycle, a third cycle, or a fold endpoint.

The early fixed-K corrector used a residual threshold proportional to K;
its largest accepted |L_z| is about 1.13e-7. Those points are approximate
fold locations, not exact rational double-cycle parameters. Arclength
acceptance has |L|,|L_z| below 2e-8 but its recorded independent derivative
discrepancy reaches 8.84e-6, near that chart's 1e-5 gate. Its two unresolved
returns and eight corrector limits record conditioning/chart failures;
they are not dynamical endpoints. The later angular and half-map repair
provides a better conditioned continuation beyond that range. The
long-double half-map residual maxima are approximately 4.25e-17 for F
and 4.29e-13 for G. Binary128 has substantially smaller saved residuals,
with no interval claim.

The nullcline and positive-horizontal sections use different r values.
For example the final arclength r approximately 214759 and the angular
r approximately 204584 describe different section coordinates near the
same parameter point. Their numerical radii must not be compared as
though the section were unchanged.

## Exact finite topology throughout the relevant rectangle

The following argument covers 24/25 <= c <= 8/5 and K>=0, containing
the requested c range and all audited fields. Write m=5(K+42)/(11c-5)
and u=1+x. Nonzero equilibria satisfy

    W(u) = (61/5-c)u^3 + (m-72/5+3c)u^2 + (11/5-3c)u + c = 0.

Here m>=50/3. For u>=0,

    W(u) >= (61/5-c)u^3 + (386/75)u^2 - (13/5)u + 24/25 > 0,

because the displayed quadratic has negative discriminant and positive
leading coefficient. Thus there is no other equilibrium on the origin
side x>-1. At x=-1 the first component is exactly 1, so this vertical
line is transverse everywhere and cannot be crossed by a closed orbit.

For u=-v<0, put h=61/5-c>0, b=m-72/5+3c>0 and a=3c-11/5>0. Then
W(-v)=-h v^3+b v^2+a v+c has exactly one positive root by its single
coefficient sign change and its endpoint signs. At that root,

    W'(u) = b v + 2a + 3c/v > 0.

The Jacobian determinant is x W'(u)/u>0. Hence there are exactly two
finite equilibria: the origin and one remote antisaddle with x<-1.
There is no finite saddle or finite equilibrium collision in this
rectangle. This finite topology is exact; it is not a root count for
limit cycles or a classification of all infinity graphics.

The origin has trace zero, determinant m>0, and the inherited first
Lyapunov sign is the sign of K. Thus every saved K>0 field has a repelling
weak focus at the origin. At K=0 the first coefficient vanishes, and at
the organizing value c=c0 the whole center annulus replaces isolated
cycles; ordinary regular-fold continuation does not pass through that
stratum unchanged.

## Exact remote stability test

For the remote point x=-(v+1),

    trace_remote = (v+1)[(1+2c)/v - (16-10c)/5].

For c<8/5 this trace is strictly decreasing in v. Its zero occurs at
v_H=5(1+2c)/(16-10c). Increasing K increases the unique remote v because
the coefficient b increases. Substitution of v_H into W(-v) gives

    K_H(c) = -441 J(c) / [125(16-10c)(1+2c)^2].

Consequently the remote antisaddle has positive trace for K<K_H(c),
zero trace for equality, and negative trace for K>K_H(c). This derivation
extends the valid test across the audited c>1.5 points; the old script's
`None` outside its inherited box does not mean the trace is unknown.
At c=8/5 the displayed trace is positive for every finite v.

Exact rational comparisons give K<K_H(c) at every audited fold and pair
field. The smallest checked gap K_H-K is about 0.002265, in the decreasing
sequence. Rational bisection additionally encloses the remote equilibrium
and its positive trace at the final accepted point of each sequence.
Selected trace values, rounded from those enclosures, are:

| Point | Remote trace |
| --- | ---: |
| K=1/8192 decreasing endpoint | 0.00031400649 |
| K=6 increasing endpoint | 7.06731246 |
| final long-double half point | 7.48330639 |
| sixth binary128 point | 7.53261397 |

These finite-point stability statements are exact for the saved decimal
parameters. They do not imply remote stability has been proved between
every continuation sample, and they do not exclude a remote limit cycle.

## What remains unproved

Positive curvature at samples, consistent pair stability, and no sampled
third root are useful numerical evidence. They do not exclude another
root between samples, outside the sampled radial range, beyond an angular
chart boundary, on an unsampled part of the parameter component, or on
another regular sheet incident to a degenerate boundary. A connected
sequence of predictors and correctors is also not a validated proof of
unique component tracking through every possible degeneracy.

No K1 certificate and no exact component exclusion is present in these
saved data. The finite center endpoint and the proposed neutral infinity
connection still require their stated local-to-global arguments. The
negative-K sheet near the center is not ruled out by the finite topology:
for sufficiently small negative K, m remains positive, origin stability
reverses, and the first-order Melnikov pair reverses stability. No exact
symmetry identifying that sheet with the positive-K sheet has been proved
in this review. It meets the latter only through the degenerate center
stratum, rather than at an ordinary regular fold point.
