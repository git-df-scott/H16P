# KKL: other infinity strata and exact remote-focus restrictions

2026-09-04. New bounded analytic note following the published
[202-call checkpoint](../STRIKE5_PRECURSOR.md). No ODE job, return
evaluation, coefficient sweep, continuation or new cycle certificate was
performed here. Earlier checkpoint notes and data are unchanged.

The exact family remains

```text
x'=(1+x)y+x^2,
y'=-10x^2+(11/5)xy+c y^2-alpha_abs x,
beta=0,
1/2<=c<=3/2, 10<=alpha_abs<=200,
K=alpha_abs(11c/5-1)-42>=1/64.
```

Below write `m=alpha_abs=-alpha`, `d=16-10c`, `e=1+2c`, and
`J(c)=305+634c-11c^2-1000c^3`. Both `d` and `e` are positive in this box.
The main new restriction is

```text
remote trace<0  if and only if  K>K_H(c),
K_H(c)=-441 J(c)/[125 d e^2].
```

This is an equilibrium statement. Negative trace alone is not a focus
certificate and does not imply a surrounding remote unstable cycle.

## 1. There is no additional finite-equilibrium branch to seed

The first nullcline is `y=-x^2/(1+x)`. The line `x=-1` has `x'=1`,
contains no equilibrium, and cannot be crossed by a periodic orbit.
Put `s=-x`. Every nonorigin equilibrium satisfies

```text
m=F(s,c),
F(s,c)=10s+(11/5)s^2/(s-1)-c s^3/(s-1)^2.
```

The K gate gives `m>420/23>18`. For `s<0`, write
`z=s/(s-1) in (0,1)`; then `F=s[10+(11/5)z-cz^2]<0`.
For `0<s<1`, `F<10`. Thus only `s>1` remains.

With `w=s-1>0`,

```text
w^3 F_s=(61/5-c)w^3+(3c-11/5)w+2c.
```

For `0<w<=1`, this is greater than `-7/10+1=3/10`.
For `w>=1`, the first two terms are at least
`(107/10-7/10)w=10w>0`. Consequently `F_s>0` everywhere on
`s>1`, throughout the specified c interval. Its endpoint limits are
`F(1+)=-infinity` and `F(+infinity)=+infinity`.

Therefore every admissible field has exactly two finite real equilibria:
the origin and one simple remote equilibrium with `x*<-1`. The same
proof remains valid at K=0 in the stated c interval. A different return
sheet cannot be explained by an overlooked finite equilibrium in this
family. This makes no claim about the number of cycles.

## 2. The remote stability boundary is tied exactly to J

At the remote equilibrium,

```text
trace = s[(1+2c)s/(s-1)-21/5],
determinant = s(s-1)F_s > 0.
```

Its trace vanishes exactly at

```text
s_H=21/d,
m_H(c)=21(1000c^2+1021c+481)/[25d e^2].
```

The derivative of the trace with respect to s is negative:

```text
trace_s=-(d/5+e/(s-1)^2)<0.
```

Since `F_s>0`, increasing m or K decreases the trace. Thus negative
trace is equivalent to `m>m_H`. Direct polynomial simplification gives

```text
m-m_H=[125d e^2 K+441J]/[25d e^2(11c-5)],
K_H=m_H(11c/5-1)-42=-441J/[125d e^2].
```

At K=0 the remote trace is negative exactly when `J>0`; at `J=0`
it is zero. On the `J<0` side a positive stability threshold must be
cleared. In particular,

```text
c=1: J=-72, m_H=973/25, K_H=588/125.
```

The inherited `K=1/64` path is therefore incompatible with a stable
remote equilibrium at c=1, independently of the numerical return-radius
cap. Moving to a negative-J shape while keeping K arbitrarily small
cannot preserve that remote stability requirement.

This boundary concerns the trace. A candidate still needs a strictly
negative focus discriminant, and the remote cycle must be found at the
same parameters. The necessary experimental region is narrowed to
`K>=1/64` and `K>K_H(c)`, with those remaining gates separate.

## 3. A remote Hopf bifurcation on positive K has the wrong cycle type

The following coefficient calculation is analytic and was independently
checked by the audit agent without ODEs or symbolic execution.

On `x<-1`, put `u=1+x<0`, `rho=(-u)^(-c-1)>0` and use a local
coordinate satisfying `dz/dx=rho`. With the notation W and B from
[the Liénard derivation](notes_lienard.md), set

```text
H=xW/u,
f(z)=-B(x),
g(z)=rho H(x).
```

The transformed equation is `z''+f(z)z'+g(z)=0`. At remote Hopf,
`f=g=0`, `g'=D=determinant>0`. Derivatives in the following coefficient
are with respect to z. In coordinates
`U=z-z*, V=-z'/sqrt(D)`, with time scaled by `sqrt(D)`,

```text
l1=(f'g''-f''g')/[16D^(3/2)].
```

To check the numerator directly in original coordinates, at the
equilibrium write `w=s-1`. Then

```text
H'=s w F_s,
H''=-2(2s-1)F_s-s w F_ss,
B'=d/5+e/w^2,
B''=2e/w^3,
f'g''-f''g'
 =[B''H'-B'H''+2(c+1)B'H'/u]/rho^2.
```

At Hopf, `s=21/d`, `w=5e/d`, and `B'=21/(5w)`. The numerator reduces
to

```text
B's[2(1-c)F_s+wF_ss]/rho^2
 = B's * 2J/[25e^2 rho^2].
```

The polynomial identity behind the last step is

```text
5e^2(1-c)(61-5c)-c(c+1)d^2=J(c).
```

Therefore

```text
l1_remote = 441J/[5000e^3 rho^2 D^(3/2)].
```

Changing the amplitude normalization changes the positive prefactor;
the sign is exactly the sign of J.

For `J<0`, the remote Hopf occurs at positive `K_H` and is
supercritical. Its small cycle is **stable**, on the side with positive
remote trace (`K<K_H`). It does not create the required unstable remote
cycle on the stable-focus side (`K>K_H`). At `J>0`, the remote Hopf is
instead at negative K, outside the precursor's K margin. At `J=0`,
the cubic vanishes and this ordinary Hopf argument is inapplicable;
the separate center-organizer analysis handles that degeneration.

Thus an ordinary remote Hopf cannot be used as a new small-U-cycle seed
inside the positive-K precursor conditions. This does **not** exclude a
remote unstable cycle of finite amplitude on another sheet.

## 4. Infinity topology before and after the two boundaries

For finite slope `z=y/x`, the infinity directions satisfy

```text
p(z)=-10+(6/5)z+(c-1)z^2=0,
Delta=(1000c-964)/25.
```

At a simple finite direction, the radial/angular eigenvalues are
`(-(1+z),p'(z))`. The vertical direction has eigenvalues `(-c,1-c)`.
The types below are projective directions; antipodal sphere charts may
reverse time, without changing saddle versus node type.

| c range | Finite directions | Vertical direction |
|---|---|---|
| `c<241/250` | None real | Saddle |
| `c=241/250` | Double direction `z=50/3`, angular eigenvalue zero | Saddle |
| `241/250<c<1` | Two positive directions: one saddle and one node | Saddle |
| `c=1` | One finite saddle `z=25/3` | Nonhyperbolic, angular eigenvalue zero |
| `1<c<=3/2` | One positive and one negative direction, both saddles | Node |

For the last row the negative root is less than `-1`, as follows from
`p(-1)=c-61/5<0`. This determines the radial eigenvalue sign.

The two new strata are not excluded by the finite-equilibrium index or
stable-focus conditions alone. Conversely a numerical return from the
first stratum cannot be continued across these boundaries while silently
retaining its old itinerary. The selected remote root reaching `-2^20`
near c=0.9301046 was an experimental limit, not proof that the true cycle
had disappeared or reached either infinity boundary.

## 5. One exact finite-amplitude checkpoint, with no existence claim

The following rational point illustrates a genuinely different amplitude
geometry just beyond c=1; it is not a numerical cycle seed:

```text
c=1001/1000, alpha=-196/5, beta=0,
K=32039/6250.
```

It was selected algebraically near the c=1 multiplier-polynomial
transition, without sampling coefficients. Exact rational evaluation gives

```text
F(88/25,c)-196/5=-102412/1771875 < 0,
F(353/100,c)-196/5=32990493/581900000 > 0.
```

Hence the remote root satisfies `88/25<s*<353/100`. Its lower endpoint
exceeds `s_H=2100/599`, so its trace is negative. Also `trace>-5` and
`determinant>79`, giving `trace^2-4determinant<-291`. Thus this is an
exact stable-focus gate beyond c=1, not evidence of a surrounding cycle.

For the multiplier quartic N in [notes_lienard.md](notes_lienard.md),
use `u=1+x`. Its coefficients, in ascending powers of u, are

```text
(3006504501, 592888296, -920742694, -122238304, 6708201)/100000000.
```

They give

```text
N(1)=32039/1250 > 0,
N(4)=-15459777419/100000000 < 0,
N(64)=76770488465461/100000000 > 0.
```

Descartes' rule supplies at most two positive roots counting multiplicity;
these signs supply two distinct roots. They are therefore simple, both
greater than one, with N negative only between them on u>0. The negative
multiplier-density band on x>0 is bounded, unlike the negative leading
term for c<1. These exact fractions were evaluated directly; no orbit
or parameter scan was involved.

This difference is a reason to consider a separate finite-amplitude
stationary-return branch, not a proof that one exists. If this checkpoint
is ever used computationally, first establish a remote U full return
within the unchanged section cap and audit the new c>1 itinerary.
Only then can origin stationary-return work be relevant to the common
precursor. N controls a necessary multiplier balance on periodic orbits;
its sign samples are **not** displacement or `D_r` signs, and do not
give root brackets for a stationary-return solver.

### Exact compactified preflight at this checkpoint

The two finite infinity directions and their radial/angular eigenvalues
are explicitly

```text
z_plus  = -600+100sqrt(37),
eigenvalues = (599-100sqrt(37), +sqrt(37)/5): saddle;
z_minus = -600-100sqrt(37),
eigenvalues = (599+100sqrt(37), -sqrt(37)/5): saddle.
```

The vertical direction is a node, with eigenvalues
`(-1001/1000,-1/1000)`. The signs follow directly from `6<sqrt(37)<7`
and, for the first radial sign, `599^2<10000*37`. The small vertical
angular eigenvalue must be retained in any near-infinity numerical
diagnostics. None of these directions is an original finite equilibrium.
This is a new c>1 starting field; no orbit is presumed to persist across
the c=1 degeneration from the published path.

Use the original downward `y=0` remote section. Its permitted component
satisfies `r<alpha/10=-98/25`, because
`Q(r,0)=r(-10r-196/5)<0` there. The existing cap remains
`-2^20<=r<-98/25`. For any attempted return, record the intermediate
upward crossing, complete downward crossing, winding about the exact
isolated remote equilibrium, coordinate ranges, and which compactified
charts and saddle/node neighborhoods were approached. A trajectory that
escapes or lacks a verified numerical full-return itinerary is unresolved;
it is not a zero or proof of nonexistence. Crossing `x=-1` also invalidates
the claimed remote periodic itinerary.

### Suggested charged control protocol, not executed

Use at most **eight** charged remote return/derivative calls, including
failures and retries. First try the four fixed section values

```text
r=-8, -512, -32768, -2^20
```

at this single exact parameter vector. These are section controls, not a
coefficient sweep or a proof of section coverage. Record D and its
corrected event derivative on every valid full return. In increasing
section order, the numerical sign orientation for a simple U root is
`D(left)<0<D(right)`. If such a bracket exists, spend the remaining at
most four calls on safeguarded refinement and derivative checking. If
there is no such bracket, a return is invalid, or the allowance is spent,
report only the tested controls and unresolved intervals. Do not expand
the section cap or infer that no remote cycle exists.

Only if this stage isolates a numerical remote U root with the correct
itinerary, a separate conditional block of at most **four** origin
return/derivative calls may use `r=1,4,16,64`, recording `D,D_r,D_rr`.
These locations sample finite amplitude relative to the exact N band;
they are not analytic stationary-root brackets. The two blocks together
have an absolute suggested ceiling of **12** charged calls, within the
unchanged overall 4096 allowance. Failure to seed a stationary branch in
these controls leaves that question open. There are no computed cycles
or displacement signs at this checkpoint in this note.

## 6. Narrow recommendation and remaining task

Keep the search within the inherited family, box and evaluation budget.
Before assigning expensive returns in another infinity stratum, enforce
the exact remote-stability threshold and focus discriminant. Do not try
to recover the required remote U orbit by its ordinary positive-K Hopf:
the local cycle has the wrong stability and lies on the wrong focus side.

The remaining construction task is still a **separate finite-amplitude
origin S/U pair** coexisting with a stable origin cycle and the remote U
cycle. Neither the existing positive displacement maximum nor the
algebraic checkpoint above supplies it. A coherent next step must seed
an actual additional stationary-return branch, check full-return
itineraries, and maintain the same-parameter remote ledger. No assertion
of existence beyond c=1, no absence theorem over the box, and no new
cycle count is made here.
