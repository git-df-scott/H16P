# Strike 5: exact first-maximum reduction on a two-root fibre

2026-09-04. This note owns the independent Green analysis for the
two-root primitive case. It uses the inherited theorem N and the
two-anchor geometry established in this strike. No numerical search,
old shooting replay, tangency search, or endpoint asymptotic is used.

## 1. The dangerous fibre and its two comparison endpoints

Fix two primitive roots `0<r<s<1`. Let `B=H_{r,s,1}` be the
three-anchor boundary primitive, with `B(1)=0`, and let
`V=alpha*K0+beta*K1+K2` be the two-root variation with the same two
interior roots. Both have signs `+,-,+` away from their roots. Write
`ell=Y0(V)>0` and `Y0(B)<0`, and define

\[
 \lambda_c=-Y_0(B)/\ell>0,\qquad C=B+\lambda_c V.
\]

The two-anchor moment-curve proof gives `eta_C>0` and `Y0(C)=0`.
Every remaining dangerous fibre point is a convex combination

\[
 H_\theta=(1-\theta)B+\theta C,\qquad0<\theta<1.
\]

It has negative center value and positive initial primitive lobe.
All original and Green coordinates depend linearly on this combination.
On `0<t<=r`,

\[
 P_C(0)=-\eta_C/192<0,\qquad
 P_C'=-\Omega C<0,\qquad
 Z_C'=P_C/(py^2)<0,\qquad Z_C(0)=0.
\]

Hence `P_C(t)<0` and `Z_C(t)<0` throughout this interval.
The boundary version of theorem N gives `Phi_B(r)<0`. Since
`Phi_B'=Rcal*Omega*B>0` before `r`, also `Phi_B(t)<0` for `0<=t<=r`.

## 2. A determinant has only upward crossings

Define

\[
 \mathcal K(t)=P_B(t)Z_C(t)-Z_B(t)P_C(t).
\]

Its derivative is exactly

\[
 \mathcal K'(t)=\Omega(t)
       [C(t)Z_B(t)-B(t)Z_C(t)]. \tag{G5-1}
\]

The terms involving `Z'=P/(py^2)` cancel. At the center,
`K(0)=-Y0(B)*P_C(0)<0`.
If `K(t)=0` at a point where `P_B(t)>0`, then

\[
 Z_B(t)=\frac{P_B(t)}{P_C(t)}Z_C(t)>0.
\]

Since `B,C>0` and `Z_C<0` on the first primitive lobe, (G5-1) gives

\[
 \boxed{\mathcal K'(t)>0
   \quad\hbox{whenever }\mathcal K(t)=0,
          \ P_B(t)>0,\ 0<t<r.} \tag{G5-2}
\]

Thus every such zero is simple and crosses upwards. In particular there
can be at most one zero on any interval where `P_B>0`, because two
upward crossings would require a downward crossing between them.

At a zero `p` of `P_theta` in this first interval, one has `P_B(p)>0`
and

\[
 Z_\theta(p)
   =-\frac{1-\theta}{P_C(p)}\mathcal K(p).
\]

The prefactor is positive. The height of the first maximum of `Z_theta`
therefore has exactly the sign of `K` at that point.

## 3. PROVED exclusion when the boundary derivative reaches zero

Suppose `P_B(r)<=0`. If `P_B(0)<=0`, then the derivative of every
dangerous convex combination is negative on the whole first lobe, so
there is no required first maximum there.

Otherwise `P_B` has a unique first zero `p_B` in `(0,r]`. At that zero,

\[
 Z_B(p_B)=\Phi_B(p_B)<0,
\]

and consequently

\[
 \mathcal K(p_B)=-Z_B(p_B)P_C(p_B)<0.
\]

Together with `K(0)<0` and the one-way crossing property (G5-2), this
proves `K(t)<0` on `(0,p_B]`. Any zero of `P_theta` on the first
primitive lobe lies strictly before `p_B`, because `P_C<0`. Mixtures
with `P_theta(0)<=0` already fail the required first gate; those with
positive initial momentum have their first zero before `p_B`.
Its corresponding first
maximum therefore has strictly negative height.

\[
 \boxed{P_B(r)\le0\ \Longrightarrow\
 \text{the entire dangerous fibre fails its required first positive
 maximum.}} \tag{G5-3}
\]

This proof covers the equality `P_B(r)=0`; no transversal boundary
crossing is needed. It excludes the first-maximum condition needed for
four original zeros in this part of the two-root geometry.

## 4. Exact reduction of the residual case to one endpoint determinant

It remains to consider `P_B(r)>0`. Then `P_B>0` on `[0,r]`, so the
one-way crossing property holds throughout `(0,r)`.
The determinant can equally be written

\[
 \boxed{\mathcal K(r)=P_B(r)\Phi_C(r)
                         -\Phi_B(r)P_C(r),} \tag{G5-4}
\]

because the `P_B*P_C*Rcal` terms cancel.

The unique mixture whose first derivative minimum is zero at `r` is

\[
 \theta_e=\frac{P_B(r)}{P_B(r)-P_C(r)}\in(0,1).
\]

Its height there has the sign of `K(r)`. If `K(r)<=0`, property
(G5-2) implies `K(t)<0` for every `0<t<r`. To include equality at
the endpoint, observe that `K(r)=0` and `P_B(r)>0` force `Z_B(r)>0`.
Then (G5-1) is strictly positive immediately to the left of `r`, so
`K` is negative there; the one-way crossing argument applies as before.
Thus every interior first
maximum of every dangerous mixture is negative. Thus a nonpositive
endpoint determinant closes the entire fibre.

Conversely, if `K(r)>0`, there is a unique zero `t_0` of `K` in
`(0,r)`. Define

\[
 \theta(t)=\frac{P_B(t)}{P_B(t)-P_C(t)}.
\]

Its derivative is

\[
 \theta'(t)=\frac{\Omega(t)
       [B(t)P_C(t)-P_B(t)C(t)]}
                {[P_B(t)-P_C(t)]^2}<0.
\]

The interval `theta_e<theta<theta(t_0)` then gives a positive first
maximum before `r`. This is only the first Green-height condition;
the remaining derivative extrema, later heights, and original-integral
primitive signs remain independent requirements for four original zeros.

Consequently, in the residual case, the necessary and sufficient condition
for existence of a positive first maximum somewhere on the dangerous
fibre is exactly `K(r)>0`. Proving `K(r)<=0` for all such baselines would
complete this first-maximum route to exclusion. No sign for (G5-4) is
asserted without an additional proof.

## 5. A useful source-ratio order

The quotient `B/V` extends analytically through their two common simple
roots. Its derivative cannot vanish at an ordinary interior point:
`B-cV` would then have a double zero there in addition to both anchors,
contradicting the inherited primitive multiplicity bound three. At an
anchor, a stationary extended quotient would instead give a triple zero
there and the other anchor, again exceeding three. Thus the derivative
never vanishes on `(0,1)`. The quotient is positive away from its
removable zeros and tends to zero at one, so

\[
 \frac{d}{dt}(B/V)<0.
\]

The same holds for `B/C=(B/V)/(B/V+lambda_c)`. This source ordering may
help bound (G5-4), but alone it is not a proved sign for that determinant.

## 6. Where a residual determinant failure would have to occur

Since `C=B+lambda_c*V`,

\[
 \mathcal K=\lambda_c(P_BZ_V-Z_BP_V).
\]

With `Y0(V)=ell`, `eta_V=-1`, and
`P0=-C_a*Y0-eta/192`, cancellation of the lift-dependent terms gives

\[
 (P_BZ_V-Z_BP_V)(0)
 =-\frac{\eta_B\ell+Y_0(B)}{192}
 =-\frac{\ell\eta_C}{192}<0.
\]

This initial determinant is independent of the lift. Its derivative is
`Omega*(V*Z_B-B*Z_V)`. As long as `Z_V>=0` on the first lobe, the
strict inequality `Z_C=Z_B+lambda_c*Z_V<0` forces `Z_B<0`, so this
derivative is negative. A positive endpoint determinant is therefore
impossible unless `Z_V` first crosses from positive to negative before
`r`. Such a crossing, if present, is unique and simple: `P_V` strictly
decreases on that lobe and `Z_V(0)=ell>0`.

In the residual case `P_B(r)>0`, `Z_B` is strictly increasing throughout
the first lobe. If `K(r)>0`, its unique interior zero has `Z_B>0`, so
the baseline must also have crossed from negative to positive before
that determinant zero. In particular the residual first-maximum lead
requires both

\[
 Z_B(r)>0,\qquad Z_V(r)<0.
\]

These are additional exact necessary conditions, not a proof that they
can occur. A universal nonpositive bound for `Z_B(r)`, or a nonnegative
bound for `Z_V(r)`, would close the remaining first-maximum case.

## 7. Verification and precise stopping point

The reciprocal mathematical review found the determinant derivative,
mixture-height sign, one-way crossing argument, both equality cases,
and source-ratio argument correct. The tiny exact replay
`q4_fifth_green_checks.py` independently verifies the determinant
derivative, the critical height `K/(P_B-P_C)`, and the derivative of
`theta(t)`. It uses lowered priority, one numerical thread, and a
ten-second CPU fuse; it evaluates no parameter points.

**PROVED:** the entire dangerous fibre is excluded when `P_B(r)<=0`;
when `P_B(r)>0`, a nonpositive endpoint determinant excludes it as well.
A positive determinant is exactly equivalent to a positive first maximum
somewhere on the fibre, and does not imply the remaining four-zero gates.

**OPEN:** the universal sign of `K(r)` in the residual domain
`P_B(r)>0`. No parameter with positive determinant is certified here,
and no global sign is inferred from the algebra or the source ordering.
Accordingly this note does not settle the outside-lobe four-zero case or
improve the established global original distinct-zero bound from four
to three.
