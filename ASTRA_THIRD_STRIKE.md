# Astra third strike: threshold crossed; no five-zero counterexample found

2026-09-04. Continued from main commit
`616519009cdec6617c2969254a44ea05f86bdb42`. All first- and second-strike
artifacts are preserved unchanged.

**No counterexample was found.** The first primitive root crosses `5/11`
along an exact analytic path, and an explicit rational coefficient box is
certified on the late-root side. Eight targeted kappa shots reached the
four-crossing intermediate pattern numerically; all failed the first
positive Green maximum. Four primitive-boundary shots and four reverse
Green-tangency lines also produced no five-zero lead.

The joint late-root/large-kappa limit has two exact cancellations. A
matched asymptotic proof nevertheless blocks its required first `P`
crossing for fixed strict anchor ratios and fixed positive scaling ratio.
This does not close the full delayed-root region or the Q4 architecture.

This strike ends under **stop G**, after the bounded construction,
boundary-shooting, reverse-tangency, and endpoint-asymptotic tasks recorded
below. Stops A-F have not been achieved. The original five-zero target
remains open; threshold crossing is not substituted for that target.

```text
FIVE Q4 ZEROS CERTIFIED: NO
FIVE-ZERO NUMERICAL CANDIDATE: NO
PRIMITIVE ROOT r > 5/11: YES
KAPPA > 21636/19043 USED: YES
POSITIVE GREEN MAXIMUM FORCED: NO
RECONSTRUCTION SATURATES +2: UNKNOWN
THREE-SIMPLE+DOUBLE POINT: NO
Q4 STILL LIVE: YES
FIVE LIMIT CYCLES REALIZED: NO
```

## Exact threshold crossing and a rational late-root certificate

**PROVED.** The root in the threshold is the **first zero of the universal
primitive**

\[
H(t)=\int_0^t uF(u)[A+Bu-1+(u-\eta)M(u)]\,du.
\]

The Green maximum is the first maximum of `Z=Y/y`, where `Y=G/C` and
`y` is the positive homogeneous solution vanishing at the loop. The two
inherited inequalities

\[
\boxed{r>5/11,\qquad \kappa>21636/19043}
\]

are strict necessary conditions for five original zeros. Equality and the
lower sides are excluded. Neither condition is sufficient.

Use the exact analytic coefficient map from the second strike:

\[
\gamma(r)=T\left(r,\frac{1+r}{2},\frac{3+r}{4}\right),\qquad0<r<1.
\]

Its three prescribed primitive roots are simple, and the entire path lies
in the strict weighted-lobe region. Its first root is identically `r`.
Thus `r_-=2/5` and `r_+=3/4` give a rigorous threshold crossing, with the
crossing exactly at `r=5/11`. This is an analytic construction, not a
numerical optimization of root positions.

**RIGOROUS COMPUTATION.** The following rational coefficients give an
explicit late-root point:

\[
\boxed{
(A,B,\eta)=\frac1{10^{18}}
(1210581187245108808,-125731163118386543,1212211767298108636).}
\]

Exact rational series with a proved tail bound certify `H` signs `+,-,+,-`
at `23/32,13/16,29/32,31/32`. The primitive therefore has exactly three
simple roots, with the first in `(23/32,13/16)`, strictly beyond `5/11`.
Every coefficient triple within `10^-8` of this point in the infinity norm
preserves those signs. The minimum base sign margin exceeds `10^-5`.
The certificate also evaluates `H(5/11)>0` directly.

Proof, exact enclosures, path formulas, and endpoint asymptotics are in
[Q4_THRESHOLD_PATH.md](Q4_THRESHOLD_PATH.md). The frozen certificate is
[third_threshold_certificate.json](q4/data/third_threshold_certificate.json).
These are primitive witnesses; they are not six alternating original-
integral witnesses.

## Construction attempts after crossing the threshold

The original reconstruction and its source-sign correction remain exactly
those of the second strike. Write

\[
P'= -\Omega H,\quad \Omega>0,\qquad
Z'=\frac{P}{p y^2},\quad Z(0)=Y_0<0,
\]
\[
X'=\frac{Y}{(1-at)^{3/2}},\quad X(0)=0,
\qquad I(s(t))=-\frac{aC}{2}\sqrt{1-at}\,X(t).
\]

Five original zeros require, sequentially, four simple `P` crossings (S1),
five simple `Y=yZ` crossings (S2), and five simple `X` crossings (S3).
At the first `P` zero before the first primitive root,

\[
Z''=-\frac{\Omega H}{p y^2}<0.
\]

Thus an actual crossing there produces an ordinary first maximum. The
unresolved issue is its positive height and the subsequent alternating
extremum heights, not merely the existence of a critical point.

**NUMERICAL.** The initial controlled probes used four path locations
`r=.5,.75,.9,.99`, each at `kappa=1.137,2,4,16`. This starts immediately
above the inherited kappa threshold. The probes bracketed scalar shooting
conditions; they were not a coefficient-space sweep.

Next, solving `P(tau2)=0` and moving a specified fraction into the adjacent
four-crossing band gave the following eight targeted shots. The second
path has anchors `(r,1-(1-r)^2,1-(1-r)^3)` and tests unequal distances
from the loop.

| Primitive-root path | r | Tuned kappa | First Z maximum, numerical |
|---|---:|---:|---:|
| affine gaps | .5 | 2.176993690 | -0.003934275289 |
| affine gaps | .75 | 2.585137225 | -0.003420052400 |
| affine gaps | .9 | 3.273485702 | -0.002994048345 |
| affine gaps | .99 | 5.563932424 | -0.002532118321 |
| affine gaps | .9999 | 8.107104320 | -0.002362670910 |
| power gaps | .6 | 2.411599781 | -0.003599582185 |
| power gaps | .75 | 2.975077238 | -0.003139178640 |
| power gaps | .9 | 4.231007907 | -0.002719393737 |

All eight exhibited the numerical (S1) pattern and failed the first (S2)
height. These floating bands are not interval certificates. The last affine
band is especially narrow, so it is retained as a diagnostic, not promoted
to a rigorous crossing count. Original-integral sample values showed no
five-sign lead; finite samples do not cover the omitted endpoint slivers.

The reproducible bounded evaluator is
[q4_green_max_3.py](q4/q4_green_max_3.py), with frozen initial, tuned, and
shape data in `q4/data/third_*_shoot.json`.

**NUMERICAL boundary tests.** Four confluent primitive contacts
`H=H'=H''=0` at `.6,.9,.99,.9999` were constructed by exact-form small
linear systems and shot in kappa to `P=0`. Every corresponding `Z` value
was below `-0.002`. This is a primitive triple contact with an auxiliary
double root, not the previously excluded auxiliary triple-root cusp.
No original three-simple-plus-double point was found.

An independent reverse construction imposed `Y(t*)=Y'(t*)=0` directly.
For fixed `(a,t*)`, these are two linear coefficient conditions, giving
one coefficient line. Four selected lines were tested. Two failed the
necessary lobe coefficient bounds with large numerical margins. The two
later lines survived those bounds but supplied no three-root primitive
point. An additional ratio-determinant analysis checked their potential
thin lobe intervals. All these statements are numerical diagnostics,
not a proof that every reverse-tangency line misses the lobe region.

## The joint endpoint escape has a proved, scoped obstruction

**PROVED.** The affine-gap path tends to

\[
(A_*,B_*,\eta_*)=(94/77,-17/77,1),\qquad Y_{0,*}=-3/1232,
\]

with positive interior primitive

\[
H_* =\frac{6t(1-t)^2}{77}[5F-36(1-t)F']>0.
\]

At the compact-interval limit `a=1`, positive-series beta integrations give

\[
\int_0^1\Omega_1 H_*\,dt=\frac1{14784},\qquad
\int_0^1(1-t)^{-2/3}\Omega_1H_*\,dt=\frac{25}{14784}.
\]

Consequently both `P_*(1)=0` and `Z_*(1)=0`, while `P_*(t)>0` and
`Z_*(t)<0` for every interior `t`. This exact balance exposed a genuine
endpoint question; it is not a finite admissible five-zero parameter.

Let `epsilon=1-r` and consider fixed strict endpoint anchor ratios and
fixed `lambda=1/(kappa epsilon)>0`. The rescaled primitive tends to

\[
h(c)=e-Vc+Qc\log c+Dc^2,\quad c=(1-t)/\epsilon,
\]

where `D>0`, `-D<e<0`, `Q<0`, and `0<V<D`. Exact Green asymptotics,
with the matching constant fixed by the original center data, yield a
strictly positive limiting first `P` minimum. Therefore sufficiently late
members of this regime cannot satisfy the first negative extremum required
by (S1). The result is locally uniform for lambda in compact subsets of
`(0,infinity)` and for fixed strict anchor ratios.

The matching is supported by an explicit integrable kernel remainder,
not by an assumed integration constant. Full proof and its independent
audit are in [Q4_GREEN_MAX_3.md](Q4_GREEN_MAX_3.md) and the linked notes.
The theorem does not cover all non-asymptotic late-root points, arbitrary
degenerating root separations, or unrestricted coupling as lambda tends
to zero or infinity. The overall positive-maximum architecture remains
unclosed.

## Verification and completion scope

**NUMERICAL hostile checks.** Two newly shot parameters were transported
back to the original four-integral basis and evaluated independently by
40-digit area quadrature. The discrepancies from the scalar reconstruction
were `1.09e-18` and `1.16e-16`. Both passed the corrected beta-strip and
cubic P2 filters. These comparisons check normalization and signs; they
are not interval proof of any original root count.

The existing eleven unit tests and four new tests pass, together with the
exact endpoint identities and all bounded new replays. The new tests
check the strict threshold arithmetic, closed moments against the positive
series, prescribed-root path residuals, and the exact rational late-root
certificate. Source hashes, commands, resource limits, and outputs are
recorded in [third_verification.txt](q4/data/third_verification.txt).

Every computation was a small linear solve, a prescribed one-dimensional
shooting or tangency calculation, or an exact identity/certificate replay.
Replays use one numerical thread, lowered priority, and ten-second CPU
ceilings. No old production sweep, original cusp experiment, old certified
box exclusion, global coefficient optimization, or literature restart ran.

The completion audit is [third_completion_audit.md](q4/third_completion_audit.md).
Candidate freezing, rigorous original sign certification, and quadratic
realization are conditional gates that were not triggered. No candidate or
realization document is supplied to imply otherwise.

**Single strongest fourth strike:** solve the reverse Green-tangency
problem `Y(t*)=Y'(t*)=0` in the remaining late-root regime, using its
one-dimensional coefficient line and the primitive ratio-determinant
boundaries to locate a certified intersection with the strict lobe region.
Use the newly proved endpoint obstruction to exclude the fixed-ratio,
fixed-lambda limit, and resolve the remaining coupled limits before further
shooting. An intersection at the first maximum, followed by a transverse
height change, would address the missing positive maximum directly;
it still must pass the remaining (S2)-(S3) signs to yield five originals.
