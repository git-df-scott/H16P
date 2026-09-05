# Q4 Strike 5: two-root geometry and the exact remaining obstruction

> Follow-on, 2026-09-05: [Strike 6](Q4_SIXTH_BOUNDARY_REDUCTION.md) proves
> new necessary bounds on the first anchor and lift, and reduces a
> sufficient determinant sign proof to two second-anchor boundary
> functions. Their universal signs remain open; the four-zero question
> is not settled. The Strike-5 derivation below is preserved.

2026-09-04. Canonical inputs: main `7def4d5`, the accepted Theorem N,
and the audited reconstruction. **The four-zero question remains open.
No four-zero counterexample or global three-zero theorem is asserted.**

The new result is an exact reduction theorem. It determines the previously
unknown center sign, excludes complete coefficient sectors, and removes
the free coefficient from the question of whether a first positive Green
maximum is possible on a two-anchor fibre. The remaining test is one
explicit determinant depending on the lift and the two anchors.

## 1. The reduction theorem

Fix a finite lift `0<a<1` and anchors `0<r<s<1`. Let

\[
B=H_{r,s,1},\qquad V=K_2-bK_0-SK_1,
\]

where `b+Sx` is the chord through the two anchor values of the inherited
strictly convex curve `x=K1/K0`, `m=K2/K0`. Thus `B` has `K3` coefficient
one and roots `r,s,1`; `V` has `K2` coefficient one and roots `r,s`.
Both have signs `+,-,+` on the three open interior lobes, and `V(1)>0`.
For the transported center functional write

\[
Y_B=\mathscr Y(B)<0,\quad v=\mathscr Y(V)>0,\quad
\lambda_c=-Y_B/v,\quad C=B+\lambda_cV.
\]

The exact center-zero endpoint has `eta_C>0`, and its Green quantities
satisfy `P_C<0,Z_C<0` throughout `0<t<=r`.

**Theorem T.** If a nonzero original Q4 integral has four distinct
interior zeros, its primitive can be scaled to `K3` coefficient one and
has the form

\[
\boxed{H=(1-\theta)B+\theta C,\quad 0<\theta<1,}
\]

for two simple anchors `r,s`. In particular `H(1)>0`, its signs are
`+,-,+`, and `Y0<0`. Necessarily

\[
\boxed{P_B(r)>0,\qquad
\mathcal K(r):=P_B(r)Z_C(r)-Z_B(r)P_C(r)>0.} \tag{R5-1}
\]

Conversely, these two strict inequalities are sufficient for some
mixtures on that fibre to have a positive first local maximum of `Z`.
They are **not** sufficient for four original zeros: the later Green
heights and all four original primitive conditions in section 4 remain
necessary.

Every fibre with `P_B(r)<=0`, or with `P_B(r)>0` and `K(r)<=0`, has at
most three distinct original interior zeros. Equality cases are included.
No sign for `K(r)` on all remaining baselines is proved in this strike.

## 2. Complete coefficient geometry and the center sign

Every normalized primitive with these anchors is uniquely

\[
H_\lambda=B+\lambda V,\qquad
\lambda=H_\lambda(1)/V(1),\qquad
Y_0=Y_B+\lambda v,\quad\eta=\eta_B-\lambda.
\]

The quotient `R=B/V`, with removable singularities at the two anchors,
is strictly decreasing from `R(0)>0` to `R(1)=0`. A critical point of
the quotient would make `B-cV` have four zeros counted with multiplicity,
contradicting the inherited primitive bound. Its sign and endpoint limit
fix the direction of monotonicity.

Therefore the two-simple-root branches are `lambda>=0`, with signs
`+,-,+`, and `lambda<=-R(0)`, with signs `-,+,-`. The intervening
coefficient interval belongs to the strict lobe region except at the two
ordinary double-contact levels. Those contacts have only one H sign
change and do not evade the zero-count reduction.

The full four-term center functional, including the term absent from the
lower-dimensional (N3) lemma, is

\[
\boxed{\mathscr Y(\alpha K_0+\beta K_1+\gamma K_2+\delta K_3)
=\frac9{3080}\left[\alpha+\frac{144}{221}\beta
 +\frac{11}{6}\gamma+\frac{204}{221}\delta\right].} \tag{R5-2}
\]

The convex-chord argument at `X=144/221` gives the exact strict bounds

\[
\boxed{\frac{231}{50312}<v<\frac3{616},\qquad
\mathscr Y(B+\eta_BV)>\frac{27}{12578}>0.} \tag{R5-3}
\]

For the second inequality use the strictly convex moment graph
`n=K3/K0`, the weighted average of `tM(t)`, whose endpoint values are
`n(0)=0,n(1)=1819/9061`. The primitive `B+eta_B V` is `K3` minus the
chord of this graph times `K0`; its `K2` coefficient is zero. The chord
at X lies strictly below the chord joining the full endpoints. The same
argument for `m`, supplemented by its endpoint tangents, proves both
bounds on v. All constants are rational.

Since `Y_B+eta_B v>0`, one obtains

\[
0<\lambda_c<\eta_B,\qquad \eta_C=\eta_B-\lambda_c>0.
\]

Consequently the center sign is **MIXED**, with exactly one transition
on each positive branch. The negative branch has `Y0<0`, matching its
initial H sign. The positive branch has `Y0<0` for `lambda<lambda_c`,
zero at `lambda_c`, and positive center value above it.

Matching center and initial H signs cannot produce the four crossings
required below. At the zero-center endpoint, `P_C(0)=-eta_C/192<0`;
the positive first lobe forces `P_C` and `Z_C` to stay negative before r.
Thus the only possible interval is `0<lambda<lambda_c`, equivalent to
the mixtures in Theorem T.

Full proofs of the chart, double and center degeneracies, strict chord
bounds, and normalized anchor monotonicity are in
[q4/notes_fifth_two_anchor.md](q4/notes_fifth_two_anchor.md).

## 3. Proof of the determinant criterion

The fixed-anchor endpoint limit of Theorem N gives the strict estimate

\[
\Phi_B(r)\le-\int_r^1W_1H_*<0.
\]

Since `Phi_B'=W_a B>0` before r, `Phi_B(t)<0` on `[0,r]`. This controls
Phi, not Z: `Z_B=Phi_B+P_B Rcal` can only be bounded after the momentum
term has been accounted for.

Define `K=P_B Z_C-Z_B P_C`. The inherited equations give

\[
\mathcal K'=\Omega(CZ_B-BZ_C),\qquad
\mathcal K(0)=-Y_BP_C(0)<0.
\]

At a zero of K with `P_B>0`, one has
`Z_B=(P_B/P_C)Z_C>0`. Because `B,C>0` and `Z_C<0` on `(0,r)`,
`K'>0` at every such zero. Hence K has at most one zero while `P_B>0`,
and every such crossing is upwards.

At a first zero `p<r` of the mixture momentum `P_theta`,

\[
\boxed{Z_\theta(p)=-\frac{1-\theta}{P_C(p)}\mathcal K(p).} \tag{R5-4}
\]

Its prefactor is positive. If `P_B(r)<=0`, either B has no positive
initial momentum, or its first momentum zero `p_B<=r` satisfies
`Z_B(p_B)=Phi_B(p_B)<0`. Then `K(p_B)<0`. The one-way crossing property
forces K to stay negative before `p_B`. Every mixture with the required
`P_theta(0)>0` has its first momentum zero before `p_B`; mixtures without
positive initial momentum already fail the first gate. This proves the
entire-fibre exclusion.

If `P_B(r)>0`, the one-way crossing property applies on the full interval
`(0,r)`. When `K(r)<=0`, K is strictly negative in its interior. To handle
`K(r)=0`, note that `Z_B(r)>0`, so `K'` is positive immediately to the
left of r; K is negative there. Thus equality also excludes a positive
first maximum.

If instead `K(r)>0`, there is a unique `t0 in (0,r)` with `K(t0)=0`.
Set

\[
\theta(t)=\frac{P_B(t)}{P_B(t)-P_C(t)},\qquad
\theta_e=\theta(r).
\]

Direct differentiation gives

\[
\theta'(t)=\frac{\Omega(BP_C-P_BC)}{(P_B-P_C)^2}<0.
\]

The explicit interval `theta_e<theta<theta(t0)` therefore gives positive
first maxima before r. This proves precisely the converse asserted in
Theorem T. It provides no automatic sign for the later extrema.

For evaluation without cancellation of the large Rcal terms,

\[
\boxed{\mathcal K(r)=P_B(r)\Phi_C(r)-\Phi_B(r)P_C(r).} \tag{R5-5}
\]

The complete argument and equality cases are in
[q4/notes_fifth_green.md](q4/notes_fifth_green.md).

## 4. Exact counterparts of (S1)–(S3)

Let the two simple H roots be `r<s` and its initial sign be sigma.
If `H(1)!=0`, then `P(1),Z(1)` both have divergent sign `-sigma`.
The first gate is

\[
\boxed{\sigma P_0>0,\quad\sigma P(r)<0,\quad\sigma P(s)>0.}
\]

It produces precisely three simple momentum roots
`p1<r<p2<s<p3`. Four Y crossings require

\[
\boxed{\sigma Y_0<0,\quad\sigma Z(p_1)>0,\quad
\sigma Z(p_2)<0,\quad\sigma Z(p_3)>0.}
\]

If their four simple roots are `v1<v2<v3<v4`, the exact original gate is

\[
\boxed{\sigma X(v_2)>0,\quad\sigma X(v_3)<0,\quad
\sigma X(v_4)>0,\quad\sigma X(1)<0.}
\]

In the surviving normalized branch sigma is positive. All the displayed
inequalities are necessary and sufficient sequentially. A first positive
Green maximum proves only one of them.

The conclusions include **distinct** original zeros of arbitrary initial
multiplicity. If `X(0)=0` and X has n distinct later zeros, each of the n
disjoint zero-to-zero intervals has zero integral of the positive weight
times Y. Nonzero analyticity forces a Y sign change in each interval.
Thus n distinct original zeros require at least n Y sign changes. The
three-piece monotonicity of P then forces every gate to be strict; a
tangency does not replace a missing crossing.

The same argument excludes H with at most one sign change, including a
double and a simple H root. It also excludes `Y0=0`. The entire projective
sector `beta1=0` is excluded: (N3) makes the center and initial H signs
agree whenever a primitive in `span{K0,K1,K2}` has two roots. If `H(1)=0`,
the endpoint momentum is finite and a fourth first-gate condition
`sigma P(1)<0` is necessary. The endpoint baseline B fails the first
height condition by the strict inherited Phi bound. Thus that complete
boundary sector is also excluded.

Details, including all endpoint cases, are in
[q4/notes_fifth_sign_chain.md](q4/notes_fifth_sign_chain.md).
An exact rational integral for the original endpoint value is derived in
[q4/notes_fifth_loop_gate.md](q4/notes_fifth_loop_gate.md).

For the first required positive original extremum there is an additional
exact functional. Define

\[
u(t)=\int_0^t\frac{y(s)}{(1-as)^{3/2}}\,ds,\quad
j(t)=\int_0^t\frac{y(s)\mathcal R(s)}{(1-as)^{3/2}}\,ds,
\quad U=u-j/\mathcal R.
\]

Then `U'=j Rcal'/Rcal^2>0` and

\[
\mathscr F=X-\frac j{\mathcal R}Z,\qquad
\mathscr F'=U'\Phi,\qquad
\mathscr F(t)=U(t)Y_0+
\int_0^t\mathcal R(s)[U(t)-U(s)]\Omega(s)H(s)\,ds.
\]

The interior kernel is positive. At a Y zero, `F=X`; hence the first
original height condition is exactly `F(v2)>0`. This functional alone
also has mixed signs on the dangerous fibre, and no uniform negative
bound for it is assumed. The derivation and additional necessary
interlacing conditions are in
[q4/notes_fifth_original_height.md](q4/notes_fifth_original_height.md).

## 5. What remains and what has not been assumed

The global distinct interior bound remains **four**. A four-zero example
must belong to a fibre satisfying (R5-1), choose its coefficient within
the first-maximum interval, and satisfy all remaining gates in section 4.
No parameter point has been certified to satisfy them here.

A proof that `K(r)<=0` whenever `P_B(r)>0` would complete the global
three-zero theorem. Conversely, a rigorously positive determinant would
only justify investigating the later gates at that fibre. It would not
by itself refute the conjecture.

The direct analogue `Phi_lambda(r)<0` on the entire dangerous interval
is false: at `lambda_c`, `Y0=0` and
`Phi_C(r)=integral_0^r W_a C>0`, so continuity gives positive Phi already
at some coefficients with negative Y0. The momentum gate prevents this
observation from being promoted to an original-zero construction.

The new normalized baseline anchor order gives an additional exact
reduction: at fixed first anchor and lift, `Z_B(r)/|Y_B|` decreases with
the second anchor. Its upper endpoint is the finite confluent
interpolation `H_{r,r,1}`. No sign of that boundary value, and no
monotonicity at a moving first anchor, has been assumed. No confluent
shots or corner expansions were performed.

All new numerical execution is limited to small exact rational or
symbolic identity checks with one thread, reduced priority, and ten-second
CPU fuses. Analytic sign arguments supply the theorems; sampled signs
play no role. Claude's files and research lanes are unchanged.
