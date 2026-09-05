# KKL pilot: exact equilibrium gates and a bounded continuation segment

2026-09-04. Independent analytic lane. No ODE solve, return evaluation,
parameter search, or cycle certificate was performed here. This note proves
algebraic gates for the authorized pilot; the existence and persistence of
the origin and remote cycles remain numerical or conditional until their
return maps are validated.

Use the exact family

```text
x' = y+x^2+xy,
y' = -10x^2+(11/5)xy+c y^2+alpha x+beta y.
```

The pilot starts at `(c,alpha,beta)=(7/10,-80,0)` and retains the inherited
ceilings: 4096 return/derivative evaluations overall, with a first pilot
limited to 64 total evaluations or 16 accepted continuation steps,
whichever occurs first. This algebraic lane spends no such evaluations.
Failed returns and derivative calls still belong in the main pilot's
ledger. These limits cannot establish coverage of every root sheet.

## 1. Exact starting-point gates

Write `A=-alpha>0`. At the origin

```text
J0 = [[0,1],[-A,beta]],
K = A((11/5)c-1)-42.
```

At the starting point `K=6/5`, with margin
`K-1/64=379/320>0`; `det J0=80` and `tr J0=0`.
The inherited first-order weak-focus calculation gives, in the specified
normalization, `l1=K/(8 sqrt(A)^3)>0`. This is a repelling order-one weak
focus on the Hopf plane. This gate alone proves no nonzero cycle.

Every nonorigin finite equilibrium has

```text
y=-x^2/(1+x),
T(x)=(c-61/5)x^3+(alpha-111/5-beta)x^2
     +(2 alpha-10-beta)x+alpha=0.
```

The excluded denominator is harmless: on `x=-1` the first component is
identically `1`, so this line contains no equilibrium and cannot be crossed
by a periodic orbit. At the starting point, put `s=-x`; then

```text
R(s)=115s^3-1022s^2+1700s-800,
R(34/5)=-1688/5 < 0,
R(69/10)=10223/200 > 0.
```

The uniqueness proof in section 2 therefore gives one and only one
nonorigin real equilibrium, with the exact isolation

```text
34/5 < s* < 69/10,
-69/10 < x* < -34/5,
1156/145 < y* < 4761/590.
```

It is a hyperbolic stable focus by the uniform inequalities below. The
interval contains no multiple equilibrium. These are exact rational
isolation and type certificates, not a validated surrounding limit cycle.

For the quadratic homogeneous terms, the finite-slope infinity equation is

```text
Q2(1,z)-z P2(1,z) = -10+(6/5)z+(c-1)z^2,
Delta_infinity = (1000c-964)/25.
```

At the starting point its discriminant is `-264/25`. There are no real
finite-slope infinity singularities. The vertical direction persists; in
the chart `(u=x/y,v=1/y)` its desingularized radial/angular eigenvalues
are `(-c,1-c)`, hence `(-7/10,3/10)` at the start. It is a saddle.
Antipodal time reversal does not change that classification.

## 2. A uniformly admissible segment on the Hopf plane

**PROVED.** The entire exact segment

```text
7/10 <= c <= 3/4,
alpha(c)=-216/(11c-5),
beta=0
```

has `K=6/5`, exactly one nonorigin real equilibrium, a hyperbolic stable
remote focus, and no change in the infinity singularity portrait.
It lies in the prescribed coefficient box. Its endpoint is
`(c,alpha)=(3/4,-864/13)`.

Here is a proof covering all parameters on that segment, rather than
sampling its endpoints. With `s=-x` and `s!=1`, the nonorigin equilibrium
equation on `beta=0` becomes

```text
A = F(s,c),
F(s,c)=10s+(11/5)s^2/(s-1)-c s^3/(s-1)^2.
```

For `s<0`, put `z=s/(s-1) in (0,1)`; then
`F=s[10+(11/5)z-c z^2]<0`, which cannot equal `A>0`.
For `0<s<1`, the second and third terms are negative, so `F<10`.
On our segment `864/13<=A<=80`, excluding these roots as well.
The only remaining possibility is `s>1`.

For `w=s-1>0`, differentiation gives

```text
F_s = 61/5-c+(3c-11/5)/w^2+2c/w^3,
w^3 F_s = (61/5-c)w^3+(3c-11/5)w+2c.
```

Uniformly in the segment,
`61/5-c>=229/20`, `3c-11/5>=-1/10`, and `2c>=7/5`.
For `0<w<=1`, the latter polynomial is strictly greater than
`-1/10+7/5=13/10`. For `w>=1`, its first two terms are at least
`(229/20-1/10)w>0`. Thus `F_s>0` for every `s>1`.
Since `F(1+,c)=-infinity` and `F(+infinity,c)=+infinity`, exactly one
remote root exists and it is simple.

Moreover,

```text
F(5,c) <= 1865/32 < 864/13 <= A,
F(7,c) >= 19397/240 > 80 >= A.
```

Consequently `5<s*<7`, uniformly. In particular the root remains on the
same side of `x=-1`, and `25/4<y*<49/6`.

At a root, the remote trace and determinant are

```text
tr J* = s[(1+2c)s/(s-1)-21/5],
det J* = s(s-1) F_s(s,c).
```

The determinant identity follows by differentiating the second component
along the first nullcline, not by treating `alpha` as varying in the
Jacobian. Since `5<s<7`,

```text
7/6 < s/(s-1) < 5/4,
-49/5 < tr J* < -43/8 < 0.
```

An alternative derivative expression is

```text
F_s = 10+(11/5) s(s-2)/(s-1)^2
      -c s^2(s-3)/(s-1)^3.
```

For `s>5`, the second term is positive, while
`0<s^2(s-3)/(s-1)^3<1`. Hence `F_s>10-c>=37/4`, and

```text
det J* > 185,
(tr J*)^2-4 det J* < 2401/25-740 = -16099/25 < 0.
```

These strict rational bounds prove stable-focus type everywhere on the
segment, with substantial margins.

For infinity,

```text
-264/25 <= Delta_infinity <= -214/25 < 0,
7/10 <= c <= 3/4 < 241/250 < 1.
```

The finite-slope directions therefore stay nonreal, while the vertical
saddle eigenvalues remain bounded away from zero. The margins to the
two named portrait boundaries are at least
`241/250-3/4=107/500` and `1-3/4=1/4`.
This proves the singularity portrait remains unchanged; it does not prove
that a numerically selected orbit retains its whole return itinerary.

## 3. Rational boundary checks for later continuation

The proven segment can be used directly without repeated root searches.
For continuation outside it, use rational isolation rather than a
floating-point equilibrium label:

1. At exact rational parameters, isolate all real roots of `T` by a
   Sturm sequence (or an equivalent exact real-root algorithm).
   Isolate the required remote root in a rational interval strictly left
   of `-1`; certify uniqueness and absence of any other nonorigin real
   root if that remains an experiment gate.
2. Substitute `y=-x^2/(1+x)` into `tr J`, `det J`, and
   `(tr J)^2-4 det J`. Rational interval evaluation on the isolating
   interval must prove respectively negative, positive, negative.
   Refine the root interval if interval dependence prevents a sign
   decision. An unresolved sign is not a type certificate.
3. Check `K>=1/64` exactly. Check the infinity boundaries exactly before
   moving to another portrait. A step crossing `c=241/250` or `c=1`
   requires a new compactified itinerary analysis; it is not an ordinary
   failed finite return.

For `beta=0`, useful exact remote-boundary parameterizations are available
without a resultant. With `s>1`, `alpha=-F(s,c)`, the trace-zero boundary is

```text
c = 8/5-21/(10s).
```

Put

```text
N(s,c)=(61-5c)(s-1)^3+(15c-11)(s-1)+10c,
E(s,c)=s[(10c-16)s+21]^2-20 N(s,c).
```

Then

```text
det J* = s N/[5(s-1)^2],
(tr J*)^2-4 det J* = s E/[25(s-1)^2].
```

Thus `N=0` is the equilibrium multiplicity/determinant boundary, and
`E=0` is the focus-node boundary. These equations specify exact algebraic
gates for a rational interval continuation cell. They are not claims that
the curves have been traced. At nonzero beta, use the original Jacobian
and cubic; the beta-zero formulas must not be reused unchanged.

## 4. Proposed path and its limits

The constant-K tangent is

```text
d alpha/dc = 2376/(11c-5)^2 > 0.
```

It avoids spending the pilot on an accidental crossing of `K=0` and stays
far from the infinity boundaries. As `c` increases along the segment,
the unique remote equilibrium also varies smoothly. In fact `s*`
decreases: implicit differentiation gives

```text
ds*/dc = [ -2376/(11c-5)^2 + s*^3/(s*-1)^2 ] / F_s < 0,
```

because the negative term has magnitude at least `38016/169>224`,
while the positive term is at most `343/36<10`.
This predicts equilibrium motion only, not a return root or its multiplier.

An optional rational mesh consistent with at most 16 accepted steps is

```text
c_j=7/10+j/320,
alpha_j=-69120/(864+11j),
beta_j=0,
j=0,...,16.
```

Adaptive smaller steps or early stopping remain appropriate; this mesh
does not override the 64-evaluation pilot ceiling or require using all
16 steps. Record both the seeded origin cycle and remote cycle at the
same parameter values, with the corrected first-return derivative.
Do not infer remote-cycle persistence from the stable-focus gate.

The path is coherent because all the exact geometry gates are preserved,
not because an analysis here predicts a fold. A missing-pair fold, a
three-nonzero-origin-cycle precursor, and coexistence with the remote
cycle remain unproved. No result in this note excludes disconnected
root sheets, certifies any cycle count, or upgrades a bounded negative
pilot to a theorem over the full coefficient box.
