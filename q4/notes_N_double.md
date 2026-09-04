# Theorem N: a global kernel comparison replaces the double-limit split

2026-09-04. Astra mathematical work on Part 3 of the fourth-strike
handoff. This is a new analytic proof, not the separate hostile-review
lane. All first-three-strike and Claude files are preserved. The proof
uses the inherited exact homogeneous equation, positive homogeneous
factor, and two corner beta moments without re-proving them.

The proposed separate double-limit estimate is unnecessary: a strict
global comparison removes the lift parameter from the entire positive
first primitive lobe. Together with the universal anchor monotonicity
proved by the other mathematical agents, it proves Theorem N.

## 1. Objects and the corrected affine constant

Let `0<a<1`, `k=1/(1-a)`, `delta=1-t`, and

\[
 y_a(t)>0,\qquad y_a(0)=1,\qquad
 p_a(t)=\sqrt{\frac{1-t}{1-at}},
\]
\[
 \mathcal R_a(t)=\int_0^t\frac{du}{p_a(u)y_a(u)^2},
 \qquad y_{2,a}=y_a\mathcal R_a.
\]

The exact homogeneous equation is

\[
 \mathscr A_a y:=(1-at)(1-t)y''-\frac{1-a}{2}y'
                            +\frac{5a}{36}y=0.
\]

Both `y_a` and `y_{2,a}` solve it. Their normalizations imply
`y_{2,a}(0)=0`, `y_{2,a}'(0)=1`.
The kernel in the new necessary condition is

\[
 W_a(t):=\mathcal R_a(t)\Omega_a(t)
 =\frac{y_{2,a}(t)}{1152t^2(1-at)^{3/2}(1-t)^{3/2}}>0.
\]

Thus `Phi_a(t)=Y0+integral_0^t W_a H`.
For reference, because
`H=(A-1)K0+B*K1-eta*K2+K3`, the constant coefficient of this affine
functional is

\[
 c_0(a,t)=-\frac{306}{1361360}
                   +\int_0^t W_a(u)[K_3(u)-K_0(u)]\,du.
\]

The fourth-strike handoff omitted `-K0` from this constant. That omission
is not used here; the original handoff is left untouched.

## 2. A scalar comparison equation

Define

\[
 v_a(t)=\frac{y_{2,a}(t)}{(1-at)^{3/2}}.
\]

Direct differentiation gives the exact conjugation

\[
 \mathscr A_a[(1-at)^{3/2}v]
             =(1-at)^{3/2}\mathscr L_a v,
\]
\[
 \boxed{\mathscr L_a v=(1-at)(1-t)v''
          -\frac{1+5a-6at}{2}v'+\frac{8a}{9}v.}
\]

Therefore `L_a v_a=0`, `v_a(0)=0`, `v_a'(0)=1`, and

\[
 W_a(t)=\frac{v_a(t)}{1152t^2(1-t)^{3/2}}.
\]

At the limiting value `a=1`, the explicit solution with the same initial
data is

\[
 v_1(t)=\frac32\left[(1-t)^{-4/3}-(1-t)^{-2/3}\right].
\]

Its residual for every finite lift has an especially simple sign:

\[
 \boxed{\mathscr L_a v_1
   =\frac{(1-a)[22-7(1-t)^{2/3}]}{6(1-t)^{7/3}}>0,\qquad0<t<1.}
\]

The short symbolic check `check_N_kernel.py` verifies only this new
conjugation, the two initial values, and the exact residual. It uses one
numerical thread, lowered priority, and a ten-second CPU ceiling. No
shooting or parameter scan is involved.

## 3. Positive causal Green kernel: no unproved maximum principle

The residual sign is sufficient because `L_a` has a known positive
homogeneous solution

\[
 h_a(t)=\frac{y_a(t)}{(1-at)^{3/2}}>0.
\]

To make the comparison precise, write
`L_a=A_a D^2+B_a D+C_a`, where `A_a=(1-at)(1-t)>0`, and let
`w=v_1-v_a`, `z=w/h_a`. Then `w(0)=w'(0)=0`, and

\[
 z''+\left(2h_a'/h_a+B_a/A_a\right)z'
          =\frac{\mathscr L_a v_1}{A_a h_a}>0.
\]

On any compact interval `[0,T]` with `T<1`, the integrating factor

\[
 \rho(t)=h_a(t)^2
           \exp\left(\int_0^t B_a(u)/A_a(u)\,du\right)>0
\]

gives `(rho*z')'>0`, with `z(0)=z'(0)=0`. Hence `z'(t)>0`, `z(t)>0`
for every `0<t<=T`. As `T<1` was arbitrary,

\[
 \boxed{0<v_a(t)<v_1(t),\qquad
        0<W_a(t)<W_1(t)\quad(0<t<1,;0<a<1).} \tag{K}
\]

This argument is the positive causal Green comparison after conjugation
by a positive homogeneous solution. It does not assume a maximum
principle for an operator with a positive zeroth-order term.

## 4. Singular endpoints and the corner baseline

The inherited corner point is

\[
 (A_*,B_*,\eta_*)=(94/77,-17/77,1),\qquad Y_{0,*}=-3/1232,
\]
\[
 H_*(t)=\frac{6t(1-t)^2}{77}
              [5F(t)-36(1-t)F'(t)]>0.
\]

The limiting kernel is

\[
 W_1(t)=\frac{3}{2304t^2}
          [(1-t)^{-17/6}-(1-t)^{-13/6}].
\]

The two inherited beta moments imply

\[
 \boxed{\int_0^1 W_1(t)H_*(t)\,dt
 =\frac32\frac{25-1}{14784}=\frac3{1232}=-Y_{0,*}.} \tag{B}
\]

There is no hidden divergent integral in this step. At the center,
`H_*(t)=O(t^2)` and `W_1(t)=O(1/t)`, so the product is `O(t)`.
At the loop, with `delta=1-t`,
`H_*=O(delta^2*(1+|log delta|))` and
`W_1=O(delta^(-17/6))`. Their product is bounded by
`C*delta^(-5/6)*(1+|log delta|)`, which is integrable.

For every fixed finite `k`, the second homogeneous solution has a finite
positive loop limit. Hence `W_a=O(delta^(-3/2))`; its product with
`H_*` is `O(delta^(1/2)*(1+|log delta|))`, also integrable. Alternatively,
(K) directly supplies the integrable majorant `W_1*H_*`.

Consequently the baseline corner sign holds for every finite lift:

\[
 \boxed{\Phi_{a,*}(1)=Y_{0,*}+\int_0^1 W_aH_*<0,
                   \qquad0<a<1.}
\]

Continuous dependence of the regular initial-value problem gives
`v_a(t)->v_1(t)` at each fixed `t<1` as `a->1`. Dominated convergence
using (K) and (B) therefore proves `Phi_{a,*}(1)->0` from below.
No conjectured decay rate or unquantified asymptotic sign is needed.

## 5. How the universal anchor comparison closes the double limit

The other mathematical agents proved the following universal facts for
every strict lobe point with first primitive root `tau1`:

\[
 Y_0<Y_{0,*},\qquad0<H(t)<H_*(t)\quad(0<t<\tau_1).
\]

Their proof uses the strictly convex curve
`x=K1/K0`, `m=K2/K0`. For a variation
`deltaH=alpha*K0+beta*K1+gamma*K2` with two distinct interior roots and
`gamma>0`, convexity gives
`alpha>-gamma/6`, `beta>-(1105/462)*gamma`, whence

\[
 \delta Y_0=\frac9{3080}
  \left[\alpha+\frac{144}{221}\beta+\frac{11}{6}\gamma\right]
 >\frac9{3080}\frac{25}{231}\gamma>0.
\]

The normalized primitive cardinal functions and the alternating signs of
`H'(tau_j)` make every anchor derivative have `gamma>0`. Both `Y0`
and the primitive at each earlier argument therefore increase as any
anchor increases. Moving all anchors toward one with their strict
endpoint-distance ratios fixed yields the displayed comparison with the
known corner limit. These signs and endpoint constants were independently
checked in this subtask; the complete universal proof belongs to the
other agents' fourth-strike notes.

Using (K), the universal comparison, and the positive integrable corner
majorant gives directly

\[
 \begin{aligned}
 \Phi_a(\tau_1)
 &=Y_0+\int_0^{\tau_1}W_aH\\
 &<Y_{0,*}+\int_0^{\tau_1}W_1H_*\\
 &<Y_{0,*}+\int_0^1W_1H_*=0.
 \end{aligned}
\]

Thus the double limit, including arbitrarily degenerating anchor ratios
and arbitrary coupling between the lift and the anchors, is covered by
a global pointwise inequality. No separate compact/corner partition is
needed for Theorem N. This statement is about the strict lobe region;
an extension to a global original zero count must separately handle
coefficient points outside that region. No such extension is assumed
in this note.
