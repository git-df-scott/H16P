# Third strike: independent Green audit and bounded reverse tangency

2026-09-04. This note contains exact identities, an independent endpoint
matching audit, and explicitly numerical diagnostics. No positive first
Green maximum or five-original-zero candidate was found here. No global
exclusion of the delayed lobe region is asserted.

## 1. The first-root threshold remains a strict necessary condition

The inherited threshold is about the first primitive root of
`H(t)=integral_0^t uF(u)q(u)du`, not a root of `q` or of the original
Abelian integral. Five distinct original zeros force five simple zeros by
the multiplicity-five bound and force the strict shooting chain (S1)--(S3).
They therefore require

\[
 t^H_1>5/11,\qquad \kappa>21636/19043.
\]

Neither inequality is sufficient for a positive Green maximum, nor do the
two inequalities together establish any original root. The excluded
complements include equality. In the dimensionless parameter
`a=1-1/kappa`, the second condition is `a>2593/21636`.

The exact first-root estimate can be checked without decimal evaluation.
Writing `m=-Y(0)>0`, the lobe bound is `m/eta<601/136136`, and
`P0=C_a*m-eta/192` with `C_a<9/4`. If the first primitive root is at most `5/11`, then
the first possible Green maximum satisfies

\[
 Z(p_1)<-m+\frac89P_0<m-\frac\eta{216}<0.
\]

Here the integral bound follows from
`y(t)>(1-t)^(5/6)` and `p(t)>=sqrt(1-t)`. Its only nontrivial rational
comparison reduces to
`11*3^11=1948617<2000000=2^7*5^6`, and
`1/216-601/136136=395/1837836>0`.

## 2. An exact moment form for the first Green maximum

Use the inherited positive homogeneous factor and integrating factor:

\[
 y(t)=\frac{O(x(s(t)))}{O(x(k))},\quad
 p(t)=\sqrt{\frac{1-t}{1-at}},\quad
 Z=Y/y,\quad P=py^2Z',\quad P'=-\Omega H,
\]

where `s=k-(k-1)t`, `x(s)=asinh(sqrt(s-1))`, and

\[
 E(x)=\frac{5\cosh(x/3)-\cosh(5x/3)}4,\qquad
 O(x)=\frac3{10}\sinh(5x/3)+\frac32\sinh(x/3).
\]

The Wronskian identity gives the following elementary primitive exactly:

\[
 \mathcal R(t):=\int_0^t\frac{du}{p(u)y(u)^2}
 =\frac{2O(x(k))^2}{\sqrt{k(k-1)}}
 \left[\frac{E(x(s(t)))}{O(x(s(t)))}
             -\frac{E(x(k))}{O(x(k))}\right].
\]

Integration by parts then yields

\[
 Z(t)=Y_0+P(t)\mathcal R(t)
       +\int_0^t\mathcal R(u)\Omega(u)H(u)\,du.
\]

At the first zero `p1` of `P`, which must precede the first `H` root
in an (S1) configuration,

\[
 P_0=\int_0^{p_1}\Omega H,\qquad
 Z(p_1)=Y_0+\int_0^{p_1}\mathcal R\Omega H.
\]

The measure `Omega*H*dt` is positive on this first interval. Thus positive
first maximum is exactly a weighted-moment inequality; its mean of
`Rcal` must exceed `(-Y0)/P0`. This identity motivates the reverse
construction below but does not by itself bound that mean for all lobes.

## 3. NUMERICAL: four fixed reverse-tangency lines

For a prescribed pair `(a,t*)`, the two conditions `Y(t*)=Y'(t*)=0`
are affine in `(A,B,eta)`. Four basis IVPs and one two-by-two solve give
the line `A=A0+A1*eta`, `B=B0+B1*eta`. The only pairs examined were

\[
 (.75,.8),\quad(.75,.95),\quad(.9,.999),\quad(.99,.999).
\]

The first two lines fail elementary necessary lobe inequalities throughout
the numerical eta strip. These are floating-point diagnoses, not interval
certificates of the whole lines. In particular, for the first line `B<-1`
throughout `1<eta<54/31`; for the second, `q'(0)<0` would require an eta
larger than the allowed strip.

For the two later lines, the remaining numerical eta intervals are
approximately `(1.02026738,1.27360457)` and
`(1.00040215,1.47941770)`. A bounded 257-value check on each line found
only zero or one sampled primitive crossing, never three.

To avoid relying solely on eta sampling, write `H=U+eta*V`. Stationary
heights of `eta(t)=-U(t)/V(t)` away from poles satisfy

\[
 U(t)q_{\rm dir}(t)-q_{\rm base}(t)V(t)=0,
\]

because `H'=tFq` and `tF>0`. The fixed 1601-point determinant mesh covers
`0.01<=t<=1-1e-10`, with loop clustering. Sign-change brackets are refined
with Brent's method. The detected values were:

| a | t* | Stationary t | Stationary eta | Endpoint H(1)=0 eta |
|---:|---:|---:|---:|---:|
| .9 | .999 | .7748605960 | 1.0126713259 | 1.1413797115 |
| .99 | .999 | .9647058127 | .9910953746 | 1.0206129074 |

The detected stationary heights lie below their necessary eta intervals.
The ratio has a detected pole at `.9927647683` on the first line and
`.9959388455` on the second. No extra thin three-root interval was
detected by this one-dimensional investigation. A finite mesh can miss
stationary roots, including even-multiplicity roots, and does not rule
out all line intersections or any other coefficient family.

The bounded replay is `q4_reverse_tangency_third.py`; the frozen numerical
record is `data/third_reverse_tangency.json`. It imposes a ten-second CPU
ceiling, lowered priority, and one numerical thread. The saved replay used about 0.111 CPU
seconds. Its IVPs, ratio roots, and reported signs are not rigorous
interval computations.

## 4. Independent audit of the exact late-root path

The path `T(r,(1+r)/2,(3+r)/4)` uses the proved global primitive-anchor
map. Its first root is exactly `r`; threshold crossing and path membership
therefore follow from the anchor theorem. The independently certified
rational box in `Q4_THRESHOLD_PATH.md` has first primitive root greater
than `23/32`, but that fact has no implication for the sign of its first
Green maximum.

The endpoint coefficient expansion in that document was checked directly
from

\[
 tFq=(tA+t^2B-\eta)F-(t-\eta)K,
 \qquad K={}_2F_1(-1/6,1/6;1;t).
\]

Its affine remainder estimate is coefficientwise. After the displayed
`U,v,E` rescaling, the three exact anchor matrices converge to the stated
invertible matrix; the largest error coefficients are `O(epsilon*L^2)`.
Thus convergence of the scaled solution follows from matrix inversion;
boundedness of that solution is not an unproved premise. The leading
constants, including `(A*,B*,eta*)=(94/77,-17/77,1)`, check.

## 5. Independent closure of the joint endpoint matching constant

This verifies the formerly outstanding remainder issue in (G5) of
`notes_green_third.md`. Fix `lambda=1/(k*epsilon)` in a compact subset
of `(0,infinity)`, and fix positive inner coordinate `c`. It does not
cover `lambda` tending to zero or infinity or degenerating anchor ratios.

Write `O_s=O(asinh(sqrt(s-1)))`,
`C_O=3/(5*2^(1/3))`, `b=5/2^(4/3)`, and

\[
 K_k(u)=\frac{k^{3/2}}{O_k}
          \frac{O_{1+(k-1)u}\sqrt u}{[1+(k-1)u]^{3/2}},
\]
\[
 T_k(u)=u^{-1/6}+bk^{-2/3}(u^{-5/6}-u^{-1/6}),
 \qquad E_k=K_k-T_k.
\]

The two exact positive-series beta moments are `1` and `25`. They give
the integral of `T_k*R` as `1+24*b*k^(-2/3)`, with
`R(u)=[5F(1-u)-36uF'(1-u)]/(1-u)`. Those terms cancel the exact
center-data expansion of `14784*P0star`; the remaining center error is
`O(k^(-1))=o(epsilon^(5/6)*L)`.

Here is a uniform bound sufficient to fix the remaining constant. The
elementary hyperbolic formula gives

\[
 n_k:=\frac{C_Ok^{5/6}}{O_k}=1-bk^{-2/3}+O(k^{-1}),
\]
\[
 \frac{O_s}{C_Os^{3/2}}
    =s^{-2/3}+bs^{-4/3}+O(s^{-5/3}).
\]

For `s=ku+1-u`, the mean-value theorem compares the first two powers
with their values at `ku`. On `epsilon*c<=u<=1`, `ku` is bounded
below by a fixed positive number. Multiplying by `k^(2/3)*sqrt(u)`
therefore gives a comparison remainder `O(k^(-1)*u^(-7/6))`.
The additional normalization errors are
`O(k^(-1)*u^(-1/6))` and `O(k^(-4/3)*u^(-5/6))`; both are bounded by
the same majorant because `u<=1` and `k>=1`. Hence

\[
 \boxed{|E_k(u)|\le C k^{-1}u^{-7/6}},
 \qquad \epsilon c\le u\le1.
\]

After `u=epsilon*v`, this becomes the integrable bound
`epsilon^(1/6)*|E_k(epsilon*v)|<=C_lambda*v^(-7/6)`.
Also `R(epsilon*v)/L` is uniformly bounded for `v>=c` and converges
pointwise to `5/(2*pi)`. Dominated convergence consequently fixes the
star contribution to the unique finite-part primitive

\[
 \mathcal B_\lambda(c)=\frac65c^{5/6}
     +6b\lambda^{2/3}c^{1/6}
     -\int_c^\infty
       [\omega_\lambda(v)v^2-v^{-1/6}
                         -b\lambda^{2/3}v^{-5/6}]\,dv.
\]

The explicit hyperbolic primitive in the PF notes has this derivative
and this large-c expansion with no constant term; it is therefore the
same primitive. Its positivity and normalization check directly.

For the coefficient perturbation, split the integral at any fixed
`0<u0<1`. On the inner piece the exact coefficientwise endpoint expansion
and the bounded scaled anchor coefficients give

\[
 2\pi\Delta H(1-u)
 =\epsilon^2L[e_\epsilon-V_\epsilon v
                     +Q_\epsilon v\log v]
 +O(\epsilon L^2u^2[1+|\log u|]),\quad v=u/\epsilon.
\]

The cancellation of the potentially larger logarithmic term uses exactly
the definition of the scaled coefficient `V`; it must not be omitted.
The exact kernel satisfies `Omega<=C*u^(-13/6)` on this inner interval.
The remainder integral is thus `O(epsilon*L^2)`, since
`integral_0^u0 u^(-1/6)*(1+|log u|)du` is finite. Dividing by
`epsilon^(5/6)*L` leaves `O(epsilon^(1/6)*L)`, which tends to zero.
On the outer piece `u>=u0`, coefficient smoothness and the center
`H/t^2` cancellation give the same `O(epsilon*L^2)` estimate. The
rescaled main integrand is bounded by
`C*v^(-13/6)*(1+v+v*|log v|)`, which is integrable for `v>=c`.

These bounds prove the full limit, including its integration constant:

\[
 \frac{2304\pi P(1-\epsilon c)}{\epsilon^{5/6}L}
 \longrightarrow
 D\mathcal B_\lambda(c)
 -\int_c^\infty\omega_\lambda(v)
          [e-Vv+Qv\log v],dv.
\]

The scaled matrix argument extends to any fixed strict anchor ratios:
the three functions `1,c,c*log(c)` form an ECT family, with Wronskian
`1/c>0`, so their three-point evaluation matrix is invertible. On compact
subsets of the strict ratio domain its inverse remains bounded.
For fixed strict anchor ratios with largest inner root one, the elementary
quadratic Rolle argument gives `-D<e<0`, `Q<0`, and `0<V<D`.
Therefore `e-Vv+Qv*log(v)<0` for `v>=1`, and the displayed limit at
`c=1` is strictly positive. In particular, the first local minimum of
`P` is positive for sufficiently small epsilon in this controlled joint
regime. The first required negative minimum in (S1) fails.

This is a rigorous asymptotic exclusion for each fixed finite positive
lambda, locally uniform on its compact subsets. It does not exclude
finite delayed configurations, `lambda` escaping those compact sets,
or anchor ratios that themselves coalesce or separate at new rates.

## 6. Final hostile consistency review

`ASTRA_THIRD_STRIKE.md`, `Q4_THRESHOLD_PATH.md`, and `Q4_GREEN_MAX_3.md`
were checked against the preceding derivations. The first-root and kappa
thresholds are strict necessary conditions; the late rational box is a
primitive certificate; the joint endpoint obstruction is restricted to
the proved scaling regime. The final documents preserve the distinction
between a positive first maximum of `Z` and a negative first minimum of
`P`. Their (S1)--(S3) equivalences include the required center and loop
endpoint signs. The original ordinary-fold claim remains conditional on
all the other strict signs and does not assert that such a point was found.

The reverse-tangency and shooting statements remain numerical, with no
rigorous whole-line or global exclusion. No stronger multiplicity bound,
five-zero candidate, five-cycle realization, or reconstruction-saturation
claim is inferred. Stop G records exhaustion of these bounded tasks rather
than closure of the remaining architecture. A companion-period notation
collision, a missing inner-interval qualifier on the Omega bound, and a
prefactor-five ambiguity in beta-sum explanatory prose were reported to
the document owners; the mathematical constants themselves check.
