# Audit of the Q4 pruning and zero certificate

Date: 2026-09-04. This audit concerns the reduced scalar functions and the
logic of certification. It does not construct five zeros or five cycles.

Primary reference: [Zhao, arXiv:1011.2253](https://arxiv.org/html/1011.2253),
especially Lemma 11, Proposition 13, Theorem 14, Proposition 17, and Section 6.
The displayed Theorem 14 and Proposition 17 disagree about a sign. The
calculations below resolve it directly. They also correct the differentiation
used in the paper's second closing comment.

## 1. Correct normalization and strip

The coefficient map in `q4_integrals.py` correctly performs the relabeling
after equation (20) before evaluating the reduced coefficients. If
`beta1 != 0`, divide **every** reduced coefficient by `beta1`, even when
`beta1 < 0`. A nonzero constant changes no zeros. The projective convention
on the original coefficients does not replace this normalization.

Write `k = kappa`, `L = k-1`, and `b = beta0/beta1` after normalization.
For `g(s)=P2(s)+(s-b)w(s)`, direct differentiation gives

\[
g'''(k)=3\frac{-25}{216L^2}
 +(k-b)\frac{775}{3888L^3}
 =-\frac{25(23k+31b-54)}{3888L^3}.
\]

Consequently the lower threshold is `(54-23*k)/31`. Proposition 17 uses
this sign; the opposite sign in Theorem 14 and Corollary 18 cannot be used
to prune parameters. The surviving strip is

\[
\frac{54-23k}{31}<b<1.
\]

It is nonempty for **every** `k>1`. In particular, it does not imply
`k<85/23`, a finite upper bound on `rho`, or an upper boundary chart at
`85/23`. Old searches made with that restriction sampled an arbitrarily
truncated parameter domain.

There is a useful covariance check. Set `r=(s-1)/L` and `B=(b-1)/L`.
The normalized Picard--Fuchs equations contain no `k`, and the strip becomes
`-23/31 < B < 0`. Equivalently, `eta=(k-b)/L=1-B` satisfies
`1 < eta < 54/31`. A supposed restriction on `k` derived solely from this
universal `g` problem would contradict that rescaling.

The code no longer rejects a small nonzero `beta1` merely because it is
less than a fixed floating tolerance. Its production coefficient map and
integral evaluator remain floating computations, so a search rejection is
still not a rigorous parameter-box exclusion.

## 2. A safe corrected polynomial filter

Assume `b<1`, and put `P=P2(b)` and `f(s)=g(s)/(s-b)`. Polynomial division
and two derivatives give the exact identity

\[
f''(s)=\frac{2P}{(s-b)^3}+w''(s).
\]

Both the factor `2` and the cubic denominator matter. Since `w'''>0`,
`w''(s)<w''(k)=-25/(216L^2)` in the open interval. If

\[
P\le\frac{25(1-b)^3}{432L^2},
\]

then `f''(s)<0` throughout that interval. For positive `P`, use
`s-b>1-b`; for nonpositive `P`, strict negativity is immediate.
Because `f(k)=0`, strict concavity bounds the number of interior zeros
by one, counting multiplicity. The zero-count chain then bounds `I` by
three. Thus a necessary survivor condition for five zeros is the strict
inequality

\[
P2(b)>\frac{25(1-b)^3}{432(k-1)^2}.
\]

The previous linear numerator `1-b` is unsupported. When `0<b<1`, it can
exclude parameters that the valid bound leaves undecided. The cubic bound
also has the correct scaling: both `P2(b)` and its threshold scale by `L`
when the normalized function is divided by `L`.

`zhao_reduced_filter` implements only this conservative bound and the
corrected strip. Fraction inputs permit exact arithmetic in the regression
tests. A passing return value means only that these particular necessary
tests have not rejected the input.

## 3. Stronger reductions for independent review

These observations are deliberately not additional operational filters.
For `b` in the corrected strip define

\[
K(s)=-\tfrac12(s-b)^3w''(s).
\]

Then

\[
f''(s)=\frac{2(P-K(s))}{(s-b)^3},\qquad
K'(s)=-\tfrac12(s-b)^2g'''(s).
\]

The positive-to-negative sign change of `g'''` makes `K` decrease to one
minimum and then increase. Its limit at `1` is positive infinity. If `g`
has three distinct interior roots, then `f` has those roots and the root
at `k`. Rolle's theorem forces two distinct interior zeros of `f''`.
Therefore the horizontal level `P` must cross both branches of `K`:

\[
\min_{1<s<k}K(s)<P<K(k)
 =\frac{25(k-b)^3}{432(k-1)^2}.
\]

The minimum occurs at the unique zero of `g'''`; its location depends only
on `k,b`, not on either coefficient of the polynomial part left free.
This is a potentially cheap one-dimensional refinement of the coarse
lower bound. Its covariance again becomes transparent after dividing by `L`.

The same shape argument gives `P2(1)<0` and `g'(k)<0` for three distinct
interior roots. Indeed, `f''` is negative, then positive, then negative;
`f'(1+)=+infinity`; and four roots of `f`, including `k`, force the last
decreasing branch to reach `f(k)=0` with negative derivative. Reading the
successive extrema backwards forces `f(1+)<0`. Equality `P2(1)=0` would
instead make `f` positive near `1` because the positive logarithmic `w`
term dominates the regular vanishing term.

Before making any of these stronger tests mandatory for **all** five-zero
integrals, explicitly carry the distinctness/multiplicity information
through the differential-operator part of the zero-count chain, or prove
the corresponding multiplicity version of the shape argument.

## 4. Six rigorous signs suffice for the scalar zero certificate

For a fixed exact nonzero coefficient vector and exact `k>1`, choose six
ordered interior points

\[
1<s_0<s_1<\cdots<s_5<k
\]

and rigorously enclose the six values of `I`. If all enclosures exclude
zero and alternate in sign, the intermediate value theorem supplies five
distinct roots in the five disjoint gaps. The multiplicity upper bound
`Z(I)<=5` then proves that these are the only interior roots and each has
multiplicity one. No derivative enclosure, interval Newton step, or
complement subdivision is logically needed for this conclusion.

The assumptions are essential: the function must be the same analytic
Abelian integral in all six evaluations, the selected component must be
the Q4 period annulus, and the multiplicity upper bound must apply to that
exact basis. Parameter boxes may be used if every sign enclosure is uniform
over the same box; six evaluations at six unrelated parameter choices do
not suffice. Floating sign changes establish none of these enclosures.

Derivative bounds and interval Newton steps remain useful for precise root
locations and quantitative perturbation control. They are optional extra
evidence at the Abelian-zero stage. The original statement that an
exactly-five claim needs the complement searched is unnecessarily strong
once the global multiplicity bound has been invoked.

This scalar certificate still does not realize a quadratic perturbation or
validate the five return-map fixed points for a specified nonzero epsilon.

## 5. Domain and historical-result cautions

A bounded domain with the strict strip, `beta1 != 0`, and other strict
conditions is not thereby compact: boundedness alone is insufficient, and
the canonical first-nonzero-positive representative is not a closed chart
of projective space. One needs explicit closed subcharts and margins for
any proof that depends on compactness.

Old stored searches retain their original numeric meaning, but their
pruning is not sound: some valid survivor regions were excluded and some
regions outside the corrected strip were admitted. They cannot support
coverage or nonexistence claims under the corrected parameterization.

Validation on 2026-09-04: the seven tests in `q4/test_q4.py` pass with
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1` and the workspace `q4-venv`
Python. The rational tests cover the corrected strip edges, absence of
the supposed `k=4` cutoff, the cubic threshold including equality, sign
and tiny nonzero projective rescaling, and exact `beta1=0` handling.
