# Stage 1: bounded independent K1 theory attempt

2026-09-05. Independent reviewer role; this is not a claim that the external
Fable model ran. One bounded session was performed, with zero orbit evaluations.
The requested two-working-session allowance has not elapsed.

**Result: K1 was not disproved.** A general order-one bound has not been proved
here. The full proof of the correctly identified order-two source was not
retrieved, so a complete line-by-line extension of that proof is still open.
The directly accessible Zhang–Cai comparison argument was tested, and an exact
failure condition plus a new amplitude restriction were established below.

## 1. Correct the source before trying to extend it

The council attribution conflates two results. [Zhang–Cai 1991,
*Quadratic systems with a weak focus*](https://doi.org/10.1017/S0004972700030008)
addresses weak/strong-focus distributions, principally the strong focus's nest.
It is not the general theorem bounding the weak focus's own nest at order two.
Its full primary PDF was available. Lemma 4 on page 515 requires a damping zero
strictly to the left of the restoring-force zero; the trace-zero origin does
not meet that hypothesis. Lemma 5 instead allows their common zero, but requires
the derivative of the damping/restoring ratio to have one strict sign on both
sides. Its proof uses the equal-potential involution and Rolle's theorem; the
contradiction uses precisely that derivative-sign hypothesis.

The general own-nest result is [Pingguang Zhang 1999,
*On the Uniqueness of the Limit Cycle of the Quadratic System with a 2nd-Order
Weak Focus*](https://doi.org/10.12386/A1999sxxb0031), Acta Mathematica Sinica 42(2),
289–304. The publisher's indexed abstract says at most one surrounding cycle,
with an additional no-cycle condition involving subsequent focus quantities.
The landing page's PDF control resolved to an icon, not the article, and direct
retrieval returned HTTP 403. I did not reconstruct an unseen proof from its
abstract. [Zhang 2002](https://doi.org/10.1007/BF02969414) independently lists this
1999 source and states the two-focus distribution result.

This distinction matters: a theorem about the other focus cannot be extended
to bound the weak focus's own nest simply by relabelling its hypotheses.

## 2. Exact test of the accessible comparison step

Use the existing KKL origin-side Liénard chart from `kkl/notes_lienard.md`:

\[
 \ddot z+f(z)\dot z+g(z)=0,\qquad
 z=(1-(1+x)^{-c})/c,\quad x>-1.
\]

Put \(u=1+x\), \(m=5(K+42)/(11c-5)\), \(d=16-10c\),
\(h=61/5-c\), and

\[
 W=m+(2m+10)x+(m+111/5)x^2+hx^3,
\]

\[
 N=\{du+(c+1)(21+dx)\}W-u(21+dx)W'.
\]

Direct differentiation gives

\[
 (f/g)_z=-\frac{u^{2c+1}N}{5W^2},\quad
 N(0)=5K,\quad [x^4]N=d(c-1)h.
\]

On the inherited box, \(W>0\) on \(x>-1\). If \(K>0\) and
\(1/2\le c<1\), the displayed leading coefficient is negative. Thus
\(N\) is positive at zero and negative sufficiently far right. The strict
one-sign derivative assumption needed for the accessible comparison proof
is false. This is a precise obstruction to **that proof**, not a construction
of three cycles and not a proof that no other comparison method works.

Changing from order two to order one restores the first nonzero focal term;
it cannot be treated as a harmless deletion of an equality. Nor can local
Bautin cyclicity bound finite cycles away from the focus.

## 3. New exact restriction near the terminal branch

For the whole rational rectangle

\[
 9/10\le c\le1,\qquad 0<K\le6/5,
\]

the following stronger statements hold:

1. \(N(x)\ge0\) for \(-1\le x\le0\).
2. The coefficients of \(x,x^2,x^3\) in \(N\) are strictly negative.
3. Its quartic coefficient is nonpositive. Hence \(N\) is strictly decreasing
   on \(x>0\), with exactly one positive root \(a(c,K)\).

The proof uses exact rational Bernstein coefficients, not a sample grid.
Multiply by the positive denominator \(5(11c-5)\), substitute
\(c=9/10+C/10\), \(K=6T/5\), and, for the left interval, \(x=-X\).
The tensor Bernstein basis is nonnegative and sums to one on its unit cube.
The script `staged_2026_09_05/theory_exact_checks.py` gives:

| Polynomial | Bernstein coefficient enclosure |
|---|---:|
| \(5(11c-5)[x]N\) | \([-2772,-8253/10]\) |
| \(5(11c-5)[x^2]N\) | \([-2664,-11627/50]\) |
| \(5(11c-5)[x^3]N\) | \([-792,-8113/200]\) |
| \(5(11c-5)N(-X)\) | \([0,945]\) |

The last enclosure uses 50 exact coefficients. Since the cubic coefficient
is strictly negative even at \(c=1\), \(N\to-\infty\) on the right throughout
the rectangle, proving the unique positive root assertion.

For a closed origin orbit its multiplier satisfies the inherited exact identity

\[
 \log M=\oint \frac{N(x)\dot x^2}{5(1+x)W(x)^2}\,dt.
\]

Consequently **every stable origin cycle and every multiplier-one origin
cycle in this rectangle must reach \(x>a(c,K)\)**. An origin-surrounding orbit
crosses \(x=0\) with nonzero velocity, giving a strictly positive contribution;
if its entire range had \(N\ge0\), its multiplier would exceed one.
The left side cannot provide the needed negative contribution here. This is
a concrete amplitude target for Stage 2 and an exact rejection test for an
alleged fold whose verified range stays at \(x\le a\).

This does not prove uniqueness: one sign change of this integrand does not
compare the orbit-dependent positive weights of distinct cycles. A new
comparison theorem or a separately verified suitable Dulac function would
be required to make that last step.

## 4. Corrections required in the stage logic

**Three cycles at trace zero need not contain a double cycle.** Three simple
positive return roots are compatible with \(D(r)=l_1r^3+O(r^4)\). A fold may
be on the boundary of the parameter region containing them; it need not be in
the field being sought. As a logic control only, the degree-nine planar field
with \(\dot\theta=1\), \(\dot r=r^3(1-r^2)(4-r^2)(9-r^2)\) has a first-order
weak focus and three hyperbolic cycles at radii 1, 2, 3. Their radial derivatives
are exactly \(-48,480,-6480\). It is explicitly **not quadratic** and supplies
no H16 counterexample. It shows why trace-zero local reasoning alone cannot
force a double cycle.

**Fold dimension gives no quick-search guarantee.** With four parameters
\(p\), the equations \(D=D_r=0\) impose two constraints in the five-dimensional
\((r,p)\)-space. At an ordinary fold their rank is two, giving dimension three.
Its parameter projection is generically a hypersurface with zero four-dimensional
volume. Direct random sampling almost surely misses the surface exactly;
Newton/continuation can still miss disconnected components or very small
basins. Ordinary-fold checks include \(D_{rr}\ne0\) and an effective parameter
derivative. Failure of one proof does not specify all possible components.

**Return-map monotonicity is not displacement monotonicity.** A planar first
return map is orientation preserving on a regular section, so \(P'>0\).
Folds require \(D'=P'-1=0\), which is compatible with \(P'>0\).
To exclude folds by monotonicity one needs a sign bound on \(D'\), or an
equivalent stronger property, over a covered return domain.

**Trace continuation is a probe, not a universal kill test.** If an innermost
cycle collapses at a nondegenerate Hopf point, the other hyperbolic cycles
persist locally. Observing that event in Shi and Chen–Wang seeds constrains
those continuations. It cannot exclude other families or disconnected fold
sheets. Degenerate equilibria, separatrix events, return-domain boundaries,
and escape to infinity are other ways a continuation can end.

**An order-one bound would have a specific scope.** It would kill the proposed
completion from three finite origin cycles at an order-one weak focus, and
the corresponding nondegenerate Hopf step. It would not, without a further
deformation or termination argument, bound every strong-focus nest or rule
out simultaneous global mechanisms or five cycles in a single nest.

**Zhang is a conditional remote-nest guard.** A two-focus candidate with at
least two cycles already in the origin nest and two in the remote nest
contradicts the distribution theorem. Two remote cycles in isolation are
not automatically a contradiction before the origin-nest premise is checked.

## 5. Stage decision and reproducibility

The bounded attempt did not establish the theory kill condition. It establishes
an exact failure of the available ratio-sign approach and a proven amplitude
restriction near the terminal branch. That permits the bounded numerical
diagnostic under an explicit **K1 unresolved** status. It does not satisfy the
requested two-session general-theorem investigation in full.

Reproduce with `python staged_2026_09_05/theory_exact_checks.py` (SymPy 1.14.0).
The machine record is `staged_2026_09_05/theory_exact_checks.json`. All assertions
pass using exact rationals; no shooting evaluations are consumed. The prior
Liénard identities were checked algebraically again; the whole-rectangle
Bernstein restriction in section 3 is new to this review.
