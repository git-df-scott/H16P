# Strike 5: the complete two-primitive-root sign chain

2026-09-04. This note inherits the audited reconstruction, auxiliary
multiplicity bound, and Theorem N. It does not repeat their proofs. Its
new work is the exact two-root shooting classification, including endpoint
and multiple-root cases. No numerical search or sampled sign is used.

## 1. Homogeneous coefficient notation

To include the direction `beta1=0`, write the full universal primitive as

\[
H=\alpha K_0+\beta K_1+\gamma K_2+\delta K_3.
\]

Here `delta` is the coefficient called `beta1` before normalization. When
`delta=1`, the usual coordinates are `alpha=A-1`, `beta=B`, and
`gamma=-eta`. Multiplication of every coefficient by a nonzero scalar
preserves the original zero set. The exact linear center data are

\[
Y_0=\frac9{3080}\alpha+\frac{162}{85085}\beta
      +\frac3{560}\gamma+\frac{459}{170170}\delta,
\]
\[
Y_1=-\frac32(1+a)Y_0+\frac\gamma{192},\qquad
P_0=Y_1-r(a)Y_0.
\]

These formulas follow from the inherited normalized data by linearity;
in particular no division by `delta` is needed when `delta=0`.

Keep the inherited positive functions `y,p,Omega,Rcal` and write

\[
Z=Y/y,\qquad P=py^2Z',\qquad P'=-\Omega H,
\]
\[
\Phi=Z-P\mathcal R=Y_0+\int_0^t\mathcal R\Omega H,
\quad\mathcal R(t)=\int_0^t\frac{du}{p(u)y(u)^2}>0,
\]
\[
X(0)=0,\qquad X'=\frac{Y}{(1-at)^{3/2}},\qquad
I=-\frac{aC}{2}\sqrt{1-at}\,X.
\]

The domain here is a fixed finite original parameter `0<a<1`; all counted
zeros lie in the open interval `0<t<1`.

## 2. Distinct original zeros force derivative sign changes

**Anchored sign-change lemma.** If the nonzero analytic function `X` has
`n` distinct interior zeros `0<x1<...<xn<1`, then `Y` has at least `n`
sign changes in `(0,1)`.

Indeed put `x0=0`. On each disjoint open interval `(x[j-1],x[j])`,

\[
0=X(x_j)-X(x_{j-1})
 =\int_{x_{j-1}}^{x_j}\frac{Y(t)}{(1-at)^{3/2}}\,dt.
\]

The positive weight and analyticity imply that `Y` takes both signs on
that interval; it cannot vanish identically on an interval for a nonzero
original coefficient direction. Each interval therefore contains a sign
change. This argument counts distinct original zeros regardless of their
multiplicity. The nonvanishing multiplier between `I` and `X` preserves
their zero sets and multiplicities.

If `H` has `k` sign changes, then `P'` has exactly those sign changes,
so `P` has at most `k+1` sign changes. Since `Z'` has the sign of `P`,
`Z`, hence `Y`, has at most `k+2` sign changes. Thus if `H` has at most
one sign change, `I` has at most three distinct interior zeros.

This includes an `H` with two distinct zeros of which one is even-order.
The inherited bound `Z(H)<=3`, counting multiplicity, implies that two
sign-changing zeros of `H` must both be simple. A triple zero plus another
sign-changing zero would already exceed that bound. Consequently only
exactly two simple `H` zeros need the classification below. The case of
three simple `H` zeros is the already closed strict lobe case, up to
nonzero scaling. A nonzero original direction cannot have `H` identically
zero, by independence and the invertible inherited coefficient transport.

## 3. Endpoint behavior when H has two simple roots

Let `0<tau1<tau2<1` be the two simple roots, and define
`sigma=sign H(t)` before the first root. Then the successive signs of `H`
are `sigma,-sigma,sigma` and those of `P'` are
`-sigma,sigma,-sigma`.

The finite-parameter homogeneous asymptotics give positive constants
`cy,cp,cR,cOmega`, depending only on `a`, with

\[
y(t)\sim c_y\sqrt{1-t},\quad p(t)\sim c_p\sqrt{1-t},
\quad\mathcal R(t)\sim c_R(1-t)^{-1/2},
\quad\Omega(t)\sim c_\Omega(1-t)^{-1}.
\]

The universal `q` is bounded at one because `M(1-)=1`; `F=O(log(1/(1-t)))`.
Hence `H(1)` exists and

\[
H(t)-H(1)=O((1-t)\log(1/(1-t))).
\]

### Nonzero endpoint primitive

If `H(1)!=0`, its sign is `sigma`. Integration of `P'=-Omega H` gives

\[
P(t)\longrightarrow-\sigma\,\infty,
\qquad Z(t)\longrightarrow-\sigma\,\infty.
\]

More precisely `P` grows logarithmically and `Z` has the same sign as
`P`, with an additional factor of order `(1-t)^(-1/2)`. The resulting
`Y` grows at most logarithmically, so `X(1)` is finite.

### Zero endpoint primitive

If `H(1)=0`, then `Omega H=O(log(1/(1-t)))` is integrable, so

\[
P_e:=P(1)=P_0-\int_0^1\Omega H
\]

is finite. Also `Rcal*Omega*H=O((1-t)^(-1/2)log(1/(1-t)))` is integrable,
so `Phi(1)` exists.

If `Pe!=0`, then `Z=P Rcal+Phi` tends to infinity with the sign of `Pe`,
and `Y` has a finite nonzero endpoint limit of that same sign. If `Pe=0`,
then

\[
P(t)=\int_t^1\Omega(u)H(u)\,du
=O((1-t)\log(1/(1-t))),
\]

so `P Rcal->0` and `Z(1)=Phi(1)` is finite. In all of these cases `X(1)`
is finite. When `Pe=0`, the sign of `P` throughout the last `H` lobe is
`sigma`, since the displayed tail integral has that sign. Such an endpoint
zero of `P` is not an interior crossing.

## 4. Necessary and sufficient two-root shooting conditions

Use `Pe=-sigma*infinity` in the case `H(1)!=0`, and the finite value from
section 3 otherwise.

### First gate: three crossings of P

The exact conditions for three distinct simple interior zeros of `P` are

\[
\boxed{\sigma P_0>0,\qquad
\sigma P(\tau_1)<0,\qquad
\sigma P(\tau_2)>0,\qquad
\sigma P_e<0.} \tag{T1}
\]

The last condition is automatic if `H(1)!=0`, leaving precisely three
inequalities in that generic case. This follows from strict monotonicity
on the three consecutive `H` lobes. The roots have the order

\[
p_1<\tau_1<p_2<\tau_2<p_3.
\]

If `H(1)=0`, the finite endpoint inequality in (T1) is essential. Its
failure, including `Pe=0`, leaves at most two interior sign changes of `P`
and therefore excludes four distinct original zeros.

### Second gate: four crossings of Z and Y

Assume (T1). Then four sign changes of `Z`, equivalently four sign changes
of `Y`, occur precisely under the four strict conditions

\[
\boxed{\sigma Y_0<0,\qquad
\sigma Z(p_1)>0,\qquad
\sigma Z(p_2)<0,\qquad
\sigma Z(p_3)>0.} \tag{T2}
\]

The endpoint sign is `-sigma` by (T1) and section 3. Thus each of the four
strictly monotone pieces contributes exactly one crossing. Every such
zero is simple, since it lies away from a zero of `P`. Denote these four
zeros by `v1<v2<v3<v4`; the signs of `Y` are
`-sigma,sigma,-sigma,sigma,-sigma`.

This determines the needed sign of the previously unknown center datum:
**four original zeros require Y0 to have sign opposite to the first H
lobe.** The conditions `Y0=0` and `sigma*Y0>0` both exclude four, regardless
of the later primitive geometry.

### Third gate: four original zeros

Assume (T1) and (T2). The exact remaining conditions are

\[
\boxed{\sigma X(v_2)>0,\qquad
\sigma X(v_3)<0,\qquad
\sigma X(v_4)>0,\qquad
\sigma X(1)<0.} \tag{T3}
\]

The first anchored extremum already has `sigma X(v1)<0`. Successive strict
monotonicity therefore proves that (T3) is sufficient for four distinct
simple original zeros. It is also necessary for four distinct original
zeros: the four intervals between consecutive extrema, beginning after
`v1`, are the only possible crossing intervals. Equality at an extremum
can merge two neighboring crossings into one zero, but cannot add a
fourth distinct zero. Equality at the endpoint leaves the final strictly
monotone lobe without an interior crossing.

Combining the anchored sign-change lemma with the three-piece structure
of `P` proves that (T1)--(T3) are necessary for **any four distinct original
zeros**, not only for a generic four-simple-zero proposal. In that case
all inequalities must be strict, and the four original zeros are in fact
simple.

For `sigma=+1`, (T1) reads `P0>0,P(tau1)<0,P(tau2)>0,Pe<0`; (T2) reads
`Y0<0,Z(p1)>0,Z(p2)<0,Z(p3)>0`; and (T3) reads
`X(v2)>0,X(v3)<0,X(v4)>0,X(1)<0`. For `sigma=-1`, reverse every displayed
sign. This covers both `+,-,+` and `-,+,-` primitive orientations.

## 5. Non-generic cases do not evade the gates

If `sigma*P0<=0`, the first monotone piece cannot contain a crossing of
`P`; the two remaining pieces allow at most two. A zero of `P` at either
primitive root is a tangency, not a sign change, because that root of `H`
is simple. Thus equality in either interior extremum condition in (T1)
cannot replace the missing crossing.

If `Y0=0` and (T1) holds, the first piece of `Z` starts at zero and moves
with sign `sigma`, so it has no interior crossing; the remaining pieces
allow at most three. If (T1) fails, the same upper bound follows from the
number of sign changes of `P`. This also handles vanishing first or higher
center derivatives: they do not create a new interior monotone piece.

A zero of `Z` at one of the simple `P` roots is an even-order tangency and
blocks the two adjacent strict sign changes needed by (T2). Once four
`Y` sign changes exist, all `Y` roots are the four simple roots described
above. Consequently an additional multiple original zero cannot hide at
an uncounted `Y` zero or evade the four inequalities in (T3).

If an `H` root is multiple, the sign-change count in section 2 applies.
For instance a double and a simple primitive root give only one `H` sign
change, hence at most three original distinct zeros. No perturbation to a
generic coefficient tuple is required for any of these conclusions.

## 6. The complete beta1=0 sector is excluded

When `delta=beta1=0`, the primitive lies in `span{K0,K1,K2}`. If it has at
most one root/sign change, section 2 already excludes four original zeros.
If it has two roots, the inherited convex moment curve shows that their
initial sign is `sign gamma` and the center functional has the same sign:

\[
\operatorname{sign}Y_0=
\operatorname{sign}H(0+)=\operatorname{sign}\gamma.
\]

For `gamma>0`, this is exactly the audited inequality

\[
Y_0=\frac9{3080}
\left[\alpha+\frac{144}{221}\beta+\frac{11}{6}\gamma\right]
>\frac9{3080}\gamma\frac{25}{231}>0;
\]

for negative `gamma`, multiply by minus one. If `gamma=0`, the normalized
primitive is affine in the increasing coordinate `K1/K0` and has at most
one root. Thus every nonzero `beta1=0` direction has at most three distinct
original zeros. In the two-root case the center-sign condition in (T2)
is violated strictly. This is a complete lower-dimensional sector, not a
continuity assertion about the normalized `beta1=1` chart.

## 7. A useful first-extremum condition and the endpoint baseline

If four original zeros existed, (T1)--(T2) would give a first critical
point `p1<tau1` with `sigma Z(p1)>0`. At that point `P(p1)=0`, hence
`Z(p1)=Phi(p1)`. Since `sigma Phi'=W_a sigma H>0` before `tau1`,

\[
\boxed{\sigma\Phi(\tau_1)>0} \tag{T4}
\]

is necessary, together with the strict endpoint and center gates above.
This is a useful obstruction, not a sufficient four-zero condition.

There is also one complete two-root boundary sector inherited from N.
For fixed `0<r<s<1`, let `B` be the `K3`-normalized primitive with roots
`r,s` and `H(1)=0`, obtained as the limit of the strict three-anchor
primitive `T(r,s,z)` as `z->1`. Its signs are `+,-,+`, and the closed
moments make the coefficient and first-lobe integrals continuous in this
limit. The quantitative Theorem N comparison, independent of the third
anchor, gives

\[
\Phi_B(r)\le-\int_r^1 W_1H_*<0.
\]

Therefore this normalized endpoint baseline cannot have four distinct
original zeros, by (T4), for any finite lift parameter. This does not use
an incorrect assumption that only the last root moving to one forces the
all-roots corner.

One must distinguish this statement from `Z_B<0` everywhere before `r`.
N controls `Phi_B=Z_B-P_B Rcal`; if `P_B` stays positive on that interval,
`Z_B` might cross without a first critical point there. Thus a proposed
convex-combination argument involving endpoint-baseline values of `Z`
needs an additional joint `(P,Z)` comparison. The negative `Phi_B` bound
alone does not provide that comparison.

## 8. Scope of this note

The reduction proves that the only remaining four-zero possibility has
exactly two simple primitive roots, the opposite center sign, and the
strict conditions (T1)--(T3). It excludes every `beta1=0` direction,
`Y0=0`, insufficient derivative crossings, all multiple-primitive-root
cases, and the normalized two-root baseline with `H(1)=0`.

The interior two-root coefficient branches and any further Green comparison
are handled by the other Strike-5 notes. Nothing here asserts a numerical
candidate or a global three-zero theorem without that final step.

## 9. Independent review

The Green-analysis and geometry agents independently checked the complete
sign-chain classification in this note, including the endpoint zero,
zero-center, projective `beta1=0`, and multiple-root cases. This agent
reciprocally checked the determinant/equality cases in
`notes_fifth_green.md` and the convex-chord, center-order, ratio, and
normalized-baseline arguments in `notes_fifth_two_anchor.md`. No
unsupported global three-zero conclusion was inferred from those
reductions.

The exact original endpoint functional used by the last condition in (T3)
is recorded separately in `notes_fifth_loop_gate.md`; no sign for that
functional on the remaining fibre is asserted.
