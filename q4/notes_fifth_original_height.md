# Strike 5: an exact functional for original extremal heights

2026-09-04. **PROVED:** the proposed original-height functional has a
strictly positive kernel and equals the original primitive at every zero
of `Y`. It gives additional exact necessary weighted-height tests.
**OPEN:** no uniform sign obstruction at the second `Y` zero is established
here. No search, quadrature, or new endpoint asymptotic is used.

## 1. Exact identity and a positive derivative

Fix a finite lift `0<a<1`. Use the inherited positive homogeneous solution
`y`, Green integral `R=Rcal`, and variables `Z,P,Phi,X`. To avoid confusion
with the fibre constant `ell(V)`, name the second integrated solution `j`:

\[
 u(t)=\int_0^t\frac{y(s)}{(1-as)^{3/2}}\,ds,\qquad
 j(t)=\int_0^t\frac{y(s)R(s)}{(1-as)^{3/2}}\,ds.
\]

Thus `u'>0`, `j'=R u'`, and `j>0` for positive time. Define

\[
 U(t)=u(t)-\frac{j(t)}{R(t)},\qquad
 \mathcal F(t)=X(t)-\frac{j(t)}{R(t)}Z(t).
\]

These formulas have removable center values `U(0)=mathcal F(0)=0`.
For positive time,

\[
 U=\frac1{R(t)}\int_0^t[R(t)-R(s)]u'(s)\,ds>0,
 \qquad U<u,
\]

and exact differentiation gives

\[
 \boxed{U'=\frac{jR'}{R^2}>0,
 \qquad\mathcal F'=U'\Phi.} \tag{O1}
\]

Indeed `X'=u'Z`, `Z'=R'P`, and `Phi=Z-RP` make all other derivative
terms cancel. In particular

\[
 \boxed{\mathcal F(t)=\int_0^t U'(s)\Phi(s)\,ds.} \tag{O2}
\]

It is a new positive-weight anchored primitive of `Phi`, not a new
independent reconstruction equation.

## 2. The positive affine forcing kernel

Substitute `Phi=Y0+int_0^s R Omega H` in (O2) and interchange integrals
on a compact interval. The center singularity is removable in their
products because `H=O(t^2)` and `R=O(t)`. The exact result is

\[
 \boxed{\mathcal F(t)=U(t)Y_0+
 \int_0^t K_t(s)\Omega(s)H(s)\,ds,\qquad
 K_t(s)=R(s)[U(t)-U(s)].} \tag{O3}
\]

Equivalently,

\[
 K_t(s)=U(t)R(s)-R(s)u(s)+j(s).
\]

Since `U` and `R` are strictly increasing from zero,
`K_t(0)=K_t(t)=0` and

\[
 \boxed{K_t(s)>0\quad(0<s<t).} \tag{O4}
\]

The proposed one-maximum description also checks directly:

\[
 \partial_sK_t(s)=R'(s)[U(t)-u(s)].
\]

Because `0<U(t)<u(t)` and `u` increases strictly from zero, the derivative
changes from positive to negative exactly once. This is a verification
of the kernel shape, not an assumption required for positivity.

For reference, the equivalent variation-of-parameters identity is

\[
 X=u\Phi+jP-
 \int_0^t[R(s)u(s)-j(s)]\Omega(s)H(s)\,ds.
\]

## 3. Connection with the first original height

At any zero of `Y`, the positive factor `y` gives `Z=0`, and therefore

\[
 \boxed{\mathcal F(t)=X(t)\quad\text{whenever }Y(t)=0.} \tag{O5}
\]

Consider the dangerous orientation `H:+,-,+`, with primitive roots
`r<s` and `Y0<0`. If the first two shooting gates hold, write the three
`P` roots as `p1<p2<p3` and the four `Y` roots as `v1<v2<v3<v4`.
The first required positive original extremal height is exactly

\[
 X(v_2)=\mathcal F(v_2)>0.
\]

Thus (O3) is an affine-in-coefficients representation of that height at
each fixed test time. The root `v2` itself still depends on the
coefficients. Since `v2<p2<s`, its forcing integral contains only the
first positive `H` lobe and possibly part of the negative middle lobe.
If `v2>r`, the latter contribution is strictly negative. One valid upper
bound is

\[
 \frac{\mathcal F(v_2)}{U(v_2)}
 <Y_0+\int_0^r R\Omega H=\Phi(r).
\]

If `v2<=r`, (O2) and the strictly increasing first-lobe `Phi` give instead
`mathcal F(v2)/U(v2)<Phi(v2)<=Phi(r)`. Thus the same strict bound holds
throughout. The remaining two-root branch already permits `Phi(r)>0`,
so this bound does not establish the desired nonpositive original height.

There is also an exact interlacing consequence. From
`Phi(pj)=Z(pj)` and `Phi(vj)=-P(vj)R(vj)`, the three simple roots of `Phi`
have the positions

\[
 \xi_1\in(v_1,p_1),\qquad
 \xi_2\in(v_2,p_2),\qquad
 \xi_3\in(v_3,p_3).
\]

Its signs are `-,+,-,+`. Equation (O1) then makes `xi1,xi2,xi3` the
successive extrema of `mathcal F`. Four original crossings necessarily
imply the additional strict height tests

\[
 \boxed{\mathcal F(\xi_2)>0,\qquad
        \mathcal F(\xi_3)<0.} \tag{O6}
\]

These follow respectively from `X(v2)>0`, `X(v3)<0`, (O5), and the
strict monotonicity between the displayed extrema. They are necessary
filters, not sufficient original-zero certificates.

## 4. Why positivity of the kernel alone does not close the fibre

Use the exact fibre endpoints from
[notes_fifth_two_anchor.md](notes_fifth_two_anchor.md): the loop-zero
baseline `B` and the center-zero point `C=B+lambda_c V`. For every fixed
`0<t<=r`, Theorem N's endpoint extension gives `Phi_B(t)<0`, while
`Y_C(0)=0` and `H_C>0` on the first lobe give `Phi_C(t)>0`. Therefore
(O2) proves

\[
 \mathcal F_B(t)<0<\mathcal F_C(t).
\]

At that fixed time `mathcal F_lambda(t)` is affine in `lambda`. It has
a unique zero strictly between zero and `lambda_c`, and becomes positive
before the center value changes sign. Consequently there is no valid
extension asserting negativity of this new functional everywhere on the
first lobe of the entire dangerous fibre. Such positivity at a fixed time
does not make that time a `Y` zero and does not verify any original-zero
gate.

The remaining meaningful target is a **joint** constraint using
`Z(v2)=0`, the first-P-root conditions, and (O3), or a sign proof for one
of (O6). The present derivation verifies an exact original-height
functional and isolates those tests; it does not supply that missing
sign proof.
