# Astra second strike: the lobe region exists, but its first crossing must be late

2026-09-04. Continued directly from first-strike commit
`16f9406004c0820363d5d6e6cf1cb964b254f8e6`.

**PROVED:** the universal weighted-lobe region is a bounded, connected
analytic cell, parametrized by any three ordered interior primitive roots.
It contains an explicitly certified rational coefficient box. The exact
reconstruction nevertheless excludes that entire box for every `kappa>1`.
Every remaining five-zero target must satisfy

\[
\boxed{\kappa>\frac{21636}{19043},\qquad
\text{first primitive zero}>\frac5{11}.}
\]

No five-zero original-integral candidate was found. This strike reaches
stop condition **C**, with a new exact reconstruction criterion and actual
bounded construction attempts beyond proving lobe nonemptiness. It does
not close the Q4 route or establish a five-cycle quadratic system.

```yaml
FIVE Q4 ZEROS CERTIFIED: NO
FIVE-ZERO NUMERICAL CANDIDATE: NO
WEIGHTED-LOBE REGION NONEMPTY: YES
THREE-SIMPLE+DOUBLE POINT: NO
RECONSTRUCTION CAN SATURATE +2: UNKNOWN
Q4 STILL LIVE: YES
```

## The universal region is now a construction coordinate system

Use the inherited normalization

\[
t=\frac{\kappa-s}{\kappa-1},\quad
F={}_2F_1(1/6,5/6;1;t),\quad
M=1-6(1-t)F'/F,
\]
\[
q=A+Bt-1+(t-\eta)M,\qquad H(t)=\int_0^t uF(u)q(u)\,du.
\]

For three simple auxiliary roots `x1<x2<x3` and initial sign `sigma`, the
three strict weighted-lobe inequalities are exactly

\[
\sigma H(x_2)<0,\qquad \sigma H(x_3)>0,\qquad \sigma H(1)<0.
\]

**PROVED.** In this normalized chart `sigma=+1`. Every ordered triple
`0<y1<y2<y3<1` determines a unique `(A,B,eta)` whose primitive has these
three simple roots. The map is a global real-analytic diffeomorphism onto
the strict lobe region. Its boundedness, connectivity, contractibility,
boundary degeneration events, relation to the auxiliary discriminant and
excluded cusp neighborhoods, and a symmetric primitive-root spine are
proved in [Q4_LOBE_REGION.md](Q4_LOBE_REGION.md).

The exact endpoint functional is particularly useful:

\[
H(1)=\frac{18(9061A+6289B-2431\eta-7242)}{85085\pi}.
\]

The final lobe condition is thus an explicit rational half-space. The other
two conditions still matter. A topological spine does not imply that every
possible five-zero reconstruction chamber meets that spine.

**RIGOROUS COMPUTATION.** An explicit interior point is

\[
\boxed{(A,B,\eta)=
\frac1{10^{12}}(1243911778077,-86917392526,1460428426173).}
\]

Every coefficient triple within `10^-7` of it in the infinity norm lies in
the strict lobe region. Exact rational series and a proved tail majorant
certify primitive signs `+,-,+,-` at `1/8,3/8,5/8,7/8`. Therefore the box
has precisely three simple primitive roots, one in each successive open
interval between these four points. These are primitive witnesses, not
six alternating-sign witnesses for the original integral.

The certificate uses only integer and rational arithmetic, includes the
verifier hash, and is replayable by
[q4/q4_lobe_certificate.py](q4/q4_lobe_certificate.py). Its full result is
[q4/data/second_lobe_certificate.json](q4/data/second_lobe_certificate.json).

## Exact reconstruction identifies the missing oscillation

**PROVED.** Put `a=(kappa-1)/kappa`, `C=pi/sqrt(kappa-1)`. The original
integral follows from the uniquely specified center IVP

\[
(1-at)(1-t)Y''-\frac{1-a}{2}Y'+\frac{5a}{36}Y
=-\frac{H}{1152t^2(1-t)},
\]
\[
Y_0=\frac{3(1326A+864B-2431\eta-102)}{1361360}<0,
\qquad Y_1=-\frac32(1+a)Y_0-\frac\eta{192},
\]
\[
X'=\frac{Y}{(1-at)^{3/2}},\quad X(0)=0,\qquad
I(s(t))=-\frac{aC}{2}\sqrt{1-at}\,X(t).
\]

The negative forcing sign corrects the source's printed positive sign
after its change of variables. The correction is proved by exact operator
conversion and checked against the independent original area integral.
It leaves the inherited zero-count inequalities unchanged.

The elementary homogeneous solutions give a positive forward Green kernel.
Dividing by its positive loop-vanishing solution `y` reduces the lift to
two weighted integrations:

\[
Z=Y/y,\quad p=\sqrt{(1-t)/(1-at)},\quad
P=p y^2Z',\quad P'=-\Omega H,\quad \Omega>0.
\]

The original center data fix `Z(0)<0` and `P(0)=Y1-y'(0)Y0`. They are not
independent free shooting parameters. The complete necessary and sufficient
sequence of strict extremum inequalities, (S1) for four `P` roots, (S2) for
five `Y` roots, and (S3) for five original roots, is written in
[Q4_RECONSTRUCTION_GEOMETRY.md](Q4_RECONSTRUCTION_GEOMETRY.md).

That file also derives compact-domain limits as `kappa` tends to `1+` or
infinity. Neither limit is silently extended uniformly to the loop.
Positivity of the Green kernel alone is not a global variation-diminishing
theorem for the original lift with its coefficient-dependent center data.

**PROVED.** If `P(0)<=0`, a lobe-region point cannot produce five original
zeros. Uniform center-data estimates exclude all original five-zero
targets for `1<kappa<=21636/19043`. A stronger bound on the first possible
positive `Z` maximum excludes every lobe-region point whose first primitive
root is at most `5/11`, for every lift parameter. The latter proof reduces
its only fractional-power comparison to the exact integer inequality
`11*3^11 < 2^7*5^6`.

**DISPROVED:** the constructed rational lobe box can saturate the lift to
five original zeros. Its first primitive root is below `3/8<5/11`; the
exclusion is uniform in `kappa`, not inferred from failed samples. The
proof also bounds original sign changes by three on this box. No global
multiplicity-three bound is claimed.

## Bounded construction attempts and the surviving fold target

**NUMERICAL.** Three primitive-anchor triples
`(r,1/2,1-r)`, with `r=1/8,1/4,3/8`, were transported at `kappa=2,4`.
All six trials passed the primitive-lobe screen. None exhibited an original
sign change at the 65 fixed samples on `[.01,.99]`. The omitted endpoint
intervals are not covered by these samples. The data and bounded replay are
[second_spine_diagnostic.json](q4/data/second_spine_diagnostic.json) and
[q4_lobe_anchors.py](q4/q4_lobe_anchors.py).

A separate one-dimensional shooting of `P(0)` for the certified rational
point selected

\[
a=\frac{189417314263391}{400000000000000},\qquad
\kappa=\frac{400000000000000}{210582685736609}.
\]

It exhibited four `P` crossings near `.144465,.382127,.583860,.810063`.
However all four corresponding `Z` extrema were near `-0.00457918`.
It reached the numerical (S1) sign pattern and decisively failed (S2).
This intermediate success is not a five-zero candidate. The frozen failure
is retained in [second_green_shoot.json](q4/data/second_green_shoot.json)
and [q4_green_shoot.py](q4/q4_green_shoot.py); the all-parameter analytic
exclusion independently rules out this coefficient box.

The scalar reconstruction also agreed with independent original area
quadrature at `kappa=2` and `4`. Those agreement checks are numerical
regressions, not interval certificates for original signs.

**PROVED conditional reduction; no explicit boundary point found.** After
(S1) and (S2), an original ordinary double zero is exactly an extremum
height `X(vj)=0`, where `Y(vj)=0` and `Y'(vj)!=0`. Replacing one of the
four interior strict inequalities in (S3) by equality, keeping the others
and `X(1)>0` strict, gives three simple original zeros plus that ordinary
double. A transverse change of the extremum height splits the double.
This replaces five arbitrary interpolation equations by one extremum-height
shooting equation with explicit strict sign guards. No such admissible
point was constructed, and no quadratic realization gate is triggered.

**CONJECTURAL / UNRESOLVED:** the remaining nested inequalities admit a
five-zero point or the guarded ordinary fold. Neither their satisfiability
nor their incompatibility is proved.

## Verification, preserved canon, and the single next strike

The first-strike structural script and all seven inherited unit tests pass.
The new exact certificate, reconstruction identities and exclusions, frozen
anchor and shooting diagnostics, and four independent second-strike tests
also pass. The latter compare the series and anchor solve with direct
hypergeometric quadrature and compare the transported PF reconstruction
with original area quadrature. Full commands, outputs, dependency versions,
and resource limits are recorded in
[q4/data/second_verification.txt](q4/data/second_verification.txt).

Each replay is single-threaded, run at lowered priority, and capped at ten
CPU seconds. No old sweep, global coefficient sampling, large symbolic
elimination, or literature restart was performed. The first-strike proofs,
original integral evaluator, corrected beta strip, and cubic P2 filter are
preserved. The withdrawn upper cutoff on `kappa` stays withdrawn.

**Single strongest third strike:** use the exact primitive-root chart with
`5/11<y1<y2<y3<1` and `kappa>21636/19043` to force the first positive Green
maximum in (S2), subject to all of (S1). Target a guarded zero extremum in
(S3), or six alternating original signs, only after that first-maximum
obstruction is overcome. Derive the delayed-root/loop asymptotics before
bounded continuation. This directly attacks the surviving reconstruction
inequality and avoids revisiting the now-excluded symmetric examples.
