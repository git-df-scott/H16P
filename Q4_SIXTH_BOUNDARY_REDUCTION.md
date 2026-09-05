# Q4 route 4: exact exclusions and two remaining boundary signs

2026-09-05. Continued from main `3b94f34`, without redoing the earlier
campaign. **Route status: NARROWED, STILL OPEN.** No four-zero Q4 example,
global Q4 three-zero theorem, or solution of `H(2)=4` is claimed.

The task was to pursue the outside-lobe four-interior-zero route and reduce
the remaining work. This strike proves two further reductions. The new
proofs build on the repository's audited Theorem N and Strike-5 Theorem T;
they have exact algebraic/arithmetic replays but have not yet received an
external mathematical review.

## What is proved

1. **Whole parameter regions are excluded.** Any surviving four-zero
   candidate must have

   \[
   r>1-(7/22)^{3/2}=0.8205212489\ldots,
   \qquad\kappa>\kappa_*=2.89924108097\ldots.
   \]

   Here `r` is the first zero of the auxiliary primitive `H`, **not** the
   first original Abelian-integral zero or a cycle coordinate. The cutoff
   `kappa_*` is defined exactly by the polynomial in the proof. Equality
   in either displayed restriction is excluded too. Every survivor also
   has `eta/(-192Y0)>19/10`.

2. **The second anchor can be removed from a sufficient sign-exclusion
   proof.** At each fixed lift `a=1-1/kappa` and first anchor `r`, the
   determinant cannot have a positive interval hidden between two
   nonpositive boundary values. It is enough to prove nonpositivity on
   the two boundary choices `s -> r+` and `s -> 1-`.

   The exact boundary determinants are explicit in the original
   reconstructed derivative basis. Apart from an identically zero case,
   the determinant has at most one zero as the second anchor varies,
   and any interior zero is simple. The proof uses a cofactor primitive:
   two extra zeros would put it in the lobe closure, where its forced
   `Phi(r)=0` contradicts Theorem N's strict negative tail bound.

The first result uses a positive Green-function comparison and two exact
rational interval bounds at finite confluent anchors. The second result
is analytic; it uses no parameter sampling or assumed monotonicity of
the determinant.

## The precise remaining target

Let `k(t)=(K0,K1,K2,K3)(t)` and
`j_a(t)=(Y_K0,Y_K1,Y_K2,Y_K3)(t)`, with the audited original-family center
data. The following two signs, for every remaining `(a,r)`, would close
route 4:

\[
\det\begin{pmatrix}k(r)\\k'(r)\\j_a(r)\\j_a'(r)\end{pmatrix}\le0,
\qquad
\det\begin{pmatrix}k(r)\\k(1)\\j_a(r)\\j_a'(r)\end{pmatrix}\le0,
\]

\[
1-(7/22)^{3/2}<r<1,
\qquad 1-1/\kappa_*<a<1.
\]

All omitted normalization factors are proved positive. **Neither
two-variable sign inequality is proved here.** A positive boundary sign
would justify testing nearby strict fibres, with the separate
`P_B(r)>0` requirement and every later height gate still mandatory.

Thus the next research task is two two-variable sign problems on a
restricted domain. A new sweep of the free mixture coefficient or the
interior second-anchor coordinate is unnecessary for an exclusion proof.
This is a reduction of the remaining problem, not removal of the entire
Q4 route from the campaign list.

## Verification and numerical evidence

| Work | Domain / result | Evidence class |
|---|---|---|
| Symbolic replay | Center limit, gauge conversion, supersolution residual, determinant orientation, lift polynomial | Exact identities |
| Confluent certificates | `r=7/10`: slope ratio exceeds `167/90`; `r=4/5`: exceeds `19/10` | Exact rational interval arithmetic with bounded infinite-series tails |
| Initial diagnostic | 60 fixed strict `(r,s,a)` points; none passes `P_B(r)>0` and `K(r)>0` | Numerical only |
| New boundary diagnostic | 40 fixed values across both boundary choices; all `D<0` | Numerical only; includes eight `a=1` comparison limits |
| Independent control | `r=.95`, `s=1`, `a=1`; quadrature and determinant ODE differ by about `7.2e-18` | Numerical cross-check, not interval certification |

The floating initial diagnostic also found positive baseline heights at
11 points. Therefore `Z_B(r)<=0` everywhere is not a viable blanket
assumption; a proof would have to contradict those computations through
an independently found error. No such error was found. These heights
do not pass the actual determinant test.

All numerical scripts use one numerical thread, lowered process priority,
and explicit CPU ceilings. The retained runs used about ten process CPU
seconds in total; there was no large search. An initial boundary run
completed its calculations but failed JSON serialization of a NumPy integer;
that output-only error was corrected and the run replayed. Successful raw
records, scripts, and dependency versions are retained.

## Files and replay

- [Exact exclusion proof](q4/sixth/notes_exclusion_wedge.md)
- [Two-boundary reduction proof](q4/sixth/notes_boundary_reduction.md)
- [Replay instructions and audit boundaries](q4/sixth/README.md)
- [Exact check script](q4/sixth/check_exact.py) and [certificate output](q4/sixth/exact_checks.json)
- [Strict-fibre diagnostic](q4/sixth/explore_determinant.py) and [raw output](q4/sixth/determinant_exploration.json)
- [Boundary diagnostic](q4/sixth/boundary_diagnostic.py) and [raw output](q4/sixth/boundary_diagnostic.json)
- [Environment](q4/sixth/environment.json)

Foundational literature remains
[Gavrilov–Iliev](https://arxiv.org/abs/0811.4602) and
[Zhao](https://arxiv.org/abs/1011.2253). Zhao's published result is an upper
bound of five annulus cycles and an existence construction of three; the
new reductions above are repository research, not claims of an already
published stronger theorem. Interior-zero bounds alone do not resolve
additional cycles at the boundary graphic or general quadratic fields.
