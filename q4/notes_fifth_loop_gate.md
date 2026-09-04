# Strike 5: exact original-loop functional for the final interior-zero gate

2026-09-04. This auxiliary identity evaluates the original endpoint functional without
an ODE reconstruction. It does not determine its sign on the remaining
fibre. Use the canonical elliptic Hamiltonian and original coefficient
basis from `Q4_THEORY.md` and `basis_mp` in `q4_integrals.py`. The
orientation is the positive area orientation used by that evaluator. Put

\[
k=\kappa>1,\quad d=k-1,\quad a=d/k,\quad c=k^{-1/2},
\quad h_{\rm loop}=-2c/3.
\]

The saddle is `(x,y)=(0,c)`. Parametrize a ray from that saddle by
`y=c+m x`. Exact substitution in the Hamiltonian gives

\[
\mathcal H(x,c+mx)-h_{\rm loop}
=x^2\left[c(km^2-d)+\frac{x}{3}D(m)\right],
\qquad D(m)=km^3-3dm+2d.
\]

For `-sqrt(a)<=m<=sqrt(a)`, the nontrivial boundary point has

\[
L(m)=\frac{3c(d-km^2)}{D(m)}\ge0.
\]

The denominator is strictly positive on this closed interval. Indeed
`D'=3k(m^2-a)<=0` there and
`D(sqrt(a))=2d(1-sqrt(a))>0`. The bounded loop domain is precisely
`0<x<L(m)` between these two tangent rays. The change of variables
`(x,m)->(x,c+mx)` has positive area Jacobian `x`. The parametrized region
is bounded, connected, and contains the center `(1,1)`: its ray slope is
`1-c` strictly between the two tangent slopes, and the center Hamiltonian
value `-2/3` is strictly below `h_loop=-2c/3`. Thus this is the bounded
loop interior of the original period annulus, rather than another branch
of the nodal cubic.

Its area moments therefore satisfy

\[
I_{00}=\frac12\int L^2dm,\qquad
I_{10}=\frac13\int L^3dm,
\]
\[
I_{01}=\int\left(\frac c2L^2+\frac m3L^3\right)dm,
\quad I_{-1,0}=\int Ldm,
\quad I_{-1,1}=\int\left(cL+\frac m2L^2\right)dm,
\]

with all integrals over `[-sqrt(a),sqrt(a)]`. In particular the apparent
singular basis terms cancel at the loop:

\[
2I_{-1,0}+3kh_{\rm loop}I_{-1,1}
=-\sqrt{k}\int mL^2dm.
\]

For the original four-term coefficient vector `mu`, define the affine
quartic

\[
N_\mu(m)=\left(-\frac{\mu_1}{3}+\frac{\mu_3}{2}-k\mu_4m\right)D(m)
 +(\mu_2+\mu_3m)(d-km^2).
\]

Then the exact original-loop functional is

\[
\boxed{I_{\rm loop}(k,\mu)=
9k^{-3/2}\int_{-\sqrt a}^{\sqrt a}
\frac{(d-km^2)^2N_\mu(m)}{D(m)^3}\,dm.} \tag{T5}
\]

The denominator is positive and the weight is nonnegative, vanishing only
at the two tangent endpoints. This is an elementary rational integral
with an affine numerator in the original coefficients. Its expanded
numerator is

\[
\begin{aligned}
N_\mu(m)={}&-k^2\mu_4m^4
+k(-\mu_1/3-\mu_3/2)m^3
+k(3d\mu_4-\mu_2)m^2\\
&+d(\mu_1-\mu_3/2-2k\mu_4)m
+d(\mu_2-2\mu_1/3+\mu_3).
\end{aligned}
\]

The bounded exact replay `check_fifth_loop_gate.py` verifies the saddle-ray Hamiltonian
factorization, positive-area Jacobian, reconstruction from all four
original area basis terms, and the polynomial expansion. It imposes one
numerical thread, reduced priority, and a ten-second CPU ceiling. It uses
no numerical evaluation or integration.

For the remaining `+,-,+` branch, the last condition in (T3) is
`X(1)<0`, equivalently `I_loop>0`, since the multiplier between the two
is strictly negative. Thus a proof that this loop functional is
nonpositive would exclude the final gate; positivity would satisfy that
gate and would not prove or disprove four interior roots. No sign claim
for (T5) on the residual fibre is made here.


This is only the final functional in the original interior-zero gate. No
statement about cycles born at the loop or endpoint cyclicity is made.
