# FASTRA — Astra afternoon block, 2026-09-05

**CE: NO. D1: OPEN. D2: OPEN.**

Base reviewed: `claude/conjecture-progress-report-ixsmgv`, commit `e3594aa`.
This work is on `astra/fastra-afternoon-2026-09-05`. Fable was not contacted or awaited.
No neutral-hemicycle, Yu–Han, Q4, negative-sheet, far-positive-sheet or random sweep was run.

## 1. Hostile review of Proposition A

**Literal proposition: GAP. Corrected implication required by D2: VERIFIED.**

The affine-divergence calculation and its polynomial factorization are correct,
but the stated converse is false, even without the two listed chart degeneracies.
There is also an omitted denominator case, `H=a^2-b(l+1)=0`.

The corrected statement is:

> For a real system in the zero-trace Shi chart with first focal value zero,
> if a non-origin equilibrium has zero divergence, then the origin is a center.

This holds including `b=0` and `l=-1`. In particular a genuine order-two weak
focus cannot coexist with a neutral non-origin saddle; hence it cannot be
surrounded by a homoclinic loop through such a saddle. No converse is needed.

### Independent elimination

Write

\[
P=-y+lx^2+mxy+y^2,\qquad Q=x(1+ax+by),\qquad
T=\operatorname{div}X=(b+2l)x+my.
\]

The undivided first-focal relation is
\[
(l+1)m=a(b+2l).
\]

For `l != -1`, substitute `m=a(b+2l)/(l+1)`. Define
\[
H=a^2-b(l+1),\qquad C_3=a^2(b+2l+1)-(b+1)(l+1)^2.
\]
At `(0,1)`, `T=m`; its vanishing puts the system either on `a=m=0`
(reversible) or on `b+2l=m=0` (Hamiltonian).

On `1+ax+by=0`, for `b != 0`,
\[
T=-\frac{(b+2l)(Hx+a)}{b(l+1)}.
\]
When `b+2l != 0` and `H != 0`, solving **the line and zero-trace equation
simultaneously** gives
\[
x_0=-a/H,\qquad y_0=(l+1)/H,\qquad
P(x_0,y_0)=-C_3/H^2.
\]
The uncancelled numerator using denominator `b^2(l+1)H^2` is exactly
\[
-b^2(l+1)C_3.
\]
That polynomial is the numerator, not the value of `P` itself.

If `b+2l=0`, divergence is identically zero. If `H=0` while `b != 0`,
`l != -1` and `b+2l != 0`, then `a != 0` and the displayed affine trace is
constant and nonzero: **there is no zero-trace candidate on the line**.

Vanishing of a factor of the second focal value does not, by itself, prove a
center. The independent audit supplies the missing sufficient center argument:
for `C3=0`, let `d=(l+1,a)^T` and
\[
R=I-2dd^T/(d^Td).
\]
The polynomial identity `X(Rz)+R X(z)=0` holds on `C3=0`.
The exact residual factors by `C3` in both components, as recorded in the JSON
certificate. This is reversibility about the line perpendicular to `d`.
With the nonsingular imaginary linear eigenvalues at the origin, reversibility
forces a center. The reflection denominator is positive for real parameters
when `l != -1`.

### The two chart degeneracies

**`b=0`, `l != -1`.** If `a=0`, then `m=0` and the origin is a reversible
center. If `l=0`, `m=0` and divergence is identically zero. In the remaining
case `a*l != 0`, the nonvertical-axis equilibria have `x=-1/a`.
Zero trace forces `y=(l+1)/a^2`, and
\[
P\left(-1/a,(l+1)/a^2\right)
=-\frac{a^2(2l+1)-(l+1)^2}{a^4}.
\]
This is the `b=0` specialization of `C3`; the same reflection proves the center.

**`l=-1`.** Do not use the divided formula for `m`. The first-focal relation is
`a(b-2)=0`, so there are two subcases.

- `b=2`: `T=my`. If `m=0`, the whole field is Hamiltonian, for every `a`.
  If `m != 0`, zero trace requires `y=0`, and then `P=-x^2`, so only the
  origin has zero trace.
- `a=0`, `b != 2`: if `m=0`, the origin is reversible. For `m != 0`, `(0,1)`
  is not neutral. If also `b=0`, `Q=x` and there is no other equilibrium.
  For `b != 0`, a neutral equilibrium on `y=-1/b` must have
  `x=m/[b(b-2)]`. Its remaining equilibrium equation is
  \[
  N=(b+1)(b-2)^2-(b-1)m^2=0,
  \qquad P=N/[b^2(b-2)^2].
  \]
  On `N=0`, the same exact reflection construction with `d=(b-2,m)^T`
  makes the field reversible. This proves a center and completes the degeneracy.

### Exact counterexample to the converse

Take `(a,b,l,m)=(1,1,0,1)`:
\[
x'=-y+xy+y^2,\qquad y'=x+x^2+xy.
\]
The origin is a reversible center under `R(x,y)=(-y,-x)`.
The only other equilibria are `(-1,0)` and `(0,1)`, with traces `-1` and `+1`.
Here `b=1`, `l=0`, `C3=0`, but `H=0`. Thus even the existential statement
“center implies some non-origin equilibrium has zero divergence” fails.

Exact coefficient vector, in order `[P:1,x,y,x²,xy,y²; Q:same]`:

```text
[0,0,-1,0,1,1, 0,1,0,1,1,0]
```

Reproduction: `python astra_afternoon_2026_09_05/audit_proposition_a.py`.
This script imports no Fable code. It checks elimination identities and explicit
reversibility, rather than only re-running the inherited factorization.

## 2. Statement (C) and the classification audit

**(C), in its intended form “no homoclinic loop/graphic surrounds a genuine
third-order weak focus”: VERIFIED from the published classification.**

Llibre–Schlomiuk, *Canadian Journal of Mathematics* 56 (2004), 310–343,
[DOI and article](https://doi.org/10.4153/CJM-2004-015-2),
[full PDF](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/FCA46BEF322C4B8A02C2C47270177DF5/S0008414X00031758a.pdf/the-geometry-of-quadratic-differential-systems-with-a-weak-focus-of-third-order.pdf):

- **Theorem 16(III), p. 337:** the only graphic portraits are **W13, W15, W18**;
  each has one graphic surrounding the **strong** focus. Parts (I) and (II)
  cover the other fifteen portraits and contain no graphics.
- **Figure 2, p. 330**, and its enlargement **Figure 3, p. 331**, were inspected
  visually. W13, W15, W18 are in Figure 3; they do not place the weak focus
  inside the graphic. Their graphic vertices are at infinity.
- **Theorem 12, p. 334:** the diagram has algebraic singularity bifurcations
  and the connection set `G3`. It explicitly qualifies the uniqueness and shape
  of `G3` as numerically observed. Do not cite the entire diagram as having
  an unqualified exact component enumeration.
- **Section 7.1(viii), p. 328:** the no-limit-cycle result around a third-order
  weak focus is attributed to Li (1986). It is not, by itself, the stronger
  no-enclosing-graphic result needed here.

The implication (C) uses Theorem 16's location of all graphics, not a claim
that numerical positioning of `G3` is exact. No portrait contradicting (C) is
present. The literal version in section 9 of the inherited D2 note, which
singles out only `(b,l)=(2,-1)`, is too narrow: center intersections must be
allowed wherever the third focal value vanishes, not just at that one point.

**Artés–Llibre–Schlomiuk (2006): full portrait/component audit NOT COMPLETED.**
The publisher's full PDF returned HTTP 403. The accessible publisher preview is
only p. 3127. The author's [supplementary page](https://mat.uab.cat/~artes/articles/qvfwf2o.html)
and publication listing supply the abstract and supplementary materials, but
explicitly do not supply the article. Their advertised 373 parameter subsets
and 126 portraits (95 with order-two focus) are not a verified list of the loop
components relevant to D2. I cannot give a checked theorem number or individual
second-order loop-portrait number from an unread catalogue.

The source is [Artés–Llibre–Schlomiuk, IJBC 16 (2006), 3127–3194](https://doi.org/10.1142/S0218127406016720).
The 2004 classification fully settles the third-order question independently of
this access limitation. The second-order catalogue still needs to be read and
mapped to the Shi chart before it can close the remaining component check.

**D2 therefore remains OPEN.** Corrected A proves `sigma != 0` on genuine
order-two loop fields; continuity makes `sign(sigma*eta2)` constant on each
connected component. Neither that fact, (C), nor 104 previously sampled negative
values proves that all components were sampled. The missing deliverable is an
exhaustive component list with a justified sign on each, including parameter
chart boundaries. No such list or global sign proof is claimed here. Likewise,
the generic square-factor observation near some center crossings is not a
substitute for handling nongeneric crossings and all components.

## 3. D1 — positive centerward sheet only

**24 exact K values from `1e-10` through `1/512`; no sign-pattern disagreement.**

The grid is

```text
1e-10, 2e-10, 5e-10, 1e-9, 2e-9, 5e-9,
1e-8, 2e-8, 5e-8, 1e-7, 2e-7, 5e-7,
1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5,
0.0001, 0.0002, 0.0004, 0.0008, 0.0012, 0.001953125.
```

Each inherited fold was corrected in binary128. As in the previous branch,
a small curvature-scaled displacement in `c` opens the pair; K is fixed and
beta is zero. The actual field is rationalized **before** any full-return call:

\[
x'=y+x^2+xy,\qquad
 y'=-\frac{5(K+42)}{11c-5}x-10x^2+\frac{11}{5}xy+cy^2.
\]

All full coefficient vectors are in the JSONL and CSV files. These are pair-side
fields, not a claim that a rounded numerical fold is an exact semistable cycle.

The ray is the positive horizontal ray at the origin; radii are measured in
the original `(x,y)` coordinates. The internal conditioning `y=6Y` does not
change the section. Fold radii range from approximately 6.757944 to 6.759395;
the chosen inherited pair-opening rule produces roots near 5.53 and 8.24.
Different sections or unfolding sizes explain why these numbers differ from
the handoff's approximate 6.8–12 description.

Two roots are bracketed on all 24 fields. The sampled local profile reaches
`|log(R/r)| = 8.43923e-17`. Root-bracket endpoints are re-evaluated at tighter
tolerance. Just outside the outer root (0.02 farther in log radius), **all
24 displacements are positive**, ranging approximately from `5.99e-15` to
`1.17e-7`. These signs agree at local tolerances `2e-25` and `2e-28`.

The initially inherited radius endpoint `e^40` still returns for every field,
and its sign is positive in all 24 cases. It is a grid endpoint, not a proven
return-domain boundary. I therefore extended the **same** 24 fields toward
loss of the full angular return. `domain_edge_map.jsonl` records the last
successful/first unresolved log-radius brackets and tighter evaluations one
log-radius unit inside the successful endpoint. These brackets range roughly
from log radius 113 at `K=1/512` to 629 at `K=1e-10`.

A numerical angular-chart failure is not proof of a mathematical domain edge:
it could require a different integration chart or more resources. All endpoint
claims are therefore about the **numerically resolved returning side**. The
initial two smallest-K tight checks hit their 250,000-evaluation cap; their
records preserve those failures and the targeted increased-budget rechecks. Both
rechecks returned successfully. The final resolved-edge comparison is **24
agreements, zero disagreements, zero unresolved sign comparisons**; unresolved
trajectory behavior beyond those endpoints remains explicitly unclassified.

The sign-map summary and exact endpoint comparison results are generated in
`astra_afternoon_2026_09_05/summary.json` and `domain_edge_map.csv`.

**This does not close D1.** Equal signs outside the known outer root and near
the domain edge permit an additional even number of roots. For example,
`D(r)=(r-a)(r-b)` is positive at both ends of an interval enclosing `a<b` but
has two interior simple zeros. A new pair need not change the endpoint signs.
No exhaustive root exclusion, derivative sign certificate, interval enclosure,
or proof that this is the only possible finite-radius birth region is claimed.

### Exact rational field for cross-verification

For `K=1/10^10`, beta zero, the smallest-K pair-side field used here is

```text
[0, 0, 1, 1, 1, 0,
 0,
 -2100000000005000000000000000000000000000000000/56548269692408210225266324857453840606122539,
 0, -10, 11/5,
 9686206335673473656842393168859440055102049/10000000000000000000000000000000000000000000]
```

Its two positive-horizontal roots are near `5.5304134744` and `8.2378533486`.
Use the JSONL brackets, not these shortened decimal values, to reproduce signs.

## 4. Precision handoff

Published self-contained files:

- `astra_afternoon_2026_09_05/full_return128.cpp`
- `astra_afternoon_2026_09_05/full_return128.py`
- `astra_afternoon_2026_09_05/README.md`
- `example_request.json` and `example_result.json` in that directory.

One call accepts a rational coefficient vector, ray angle and log-radius grid;
it returns decimal-string log/radial displacements and explicit statuses.
GCC/libquadmath is compiled automatically into a temporary cache. No conversion
through double is used for coefficients or input coordinates. Optional exact
rational equilibrium translation and automatic anisotropic scaling are provided.

The interface rejects floating-point inputs and nonfinite intermediate/error
values; failures cannot masquerade as zero displacement. It integrates a full
turn, rather than inferring a return from matching two half-maps. Angular
monotonicity is required; failures are explicitly unresolved. It is numerical,
not interval arithmetic.

Six analytic/failure controls passed, covering near-integrable displacement,
radii beyond `e^36`, large coefficients, arbitrary-ray anisotropic centers,
a nonlinear exact center, and a failed angular chart. The actual 24 nonlinear
pair fields were checked at two tolerances as described above. The handoff is
ready for targeted calls; no additional field sweep was launched.
