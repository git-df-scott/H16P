# Independent Astra review: KKL first-order Hopf completion

Council review only, 2026-09-04. No return-map construction, coefficient scan,
continuation, or numerical replay was run in this review. This is Astra's
independent assessment, not an attribution to Fable or another model.

## 1. Evidence and exact target

The recovered audit packet is at
`/Users/scottg/Documents/Codex/2026-09-04/your-auditing-now-for-h16-x20/work/audit-only-publication`,
commit `de39ea78d56208a2a3267b594ce5c117b6b14c1e`. This review inspected
`ASTRA_FIFTH_STRIKE_HANDOFF.md`, `FOUR_CYCLE_SEED_LEDGER.md`,
`H16_CERTIFICATION_PLAN.md`, and the scripts, JSON controls and red-team report
in `frontier_2026_09_04`. The recovered packet supersedes the earlier statement
that the post-Q4 handoff was unavailable. Its bounded box is an existing
experiment specification, not a box invented by this council.

The proposed field is exactly

    x' = y + x^2 + xy,
    y' = -10x^2 + (11/5)xy + c y^2 + alpha x + beta y.

The precursor has `beta=0`, `1/2 <= c <= 3/2`,
`-200 <= alpha <= -10`, and

    K = -alpha((11/5)c-1)-42 >= 1/64.

Require three **nonzero**, mutually separated hyperbolic origin cycles in
S/U/S order, plus a remote hyperbolic U cycle, at this same parameter pair.
These are four preexisting finite cycles. Merely finding three origin cycles
at one shape and a remote cycle at another supplies no precursor.

At the origin `J=[[0,1],[alpha,beta]]`. On `beta=0`, put
`omega=sqrt(-alpha)`, `u=x`, `v=-y/omega`, `tau=omega*t`. The cubic radial
coefficient in that normalization is

    l1 = K/(8 omega^3) > 0.

The real part of the physical eigenvalues has beta derivative `1/2`.
Consequently an ordinary Hopf bifurcation at a sufficiently small negative
beta produces an unstable inner cycle. Hyperbolicity and finite separation
make the four old cycles persist for sufficiently small beta. The resulting
origin order is U/S/U/S, with the remote U cycle unchanged. This conditional
completion is valid. Its unproved premise is the four-cycle precursor.

## 2. Published seed versus the proposed precursor

KKL's reduced equations are the family above. Its displayed numerical point
is `(-10,2.2,0.7,-72.7778,0.0015)`; the published construction uses finite
cycles plus local weak-focus bifurcations, and its first-order Hopf route
starts with three finite cycles **total**, then adds one. The decimal example
is visualized with three origin cycles and one remote cycle. This does not
supply four finite cycles on the Hopf plane. [KKL, 2013, equations (2.2),
Theorems 1–2 and Fig. 4](https://doi.org/10.1007/s12591-012-0118-6);
[author-uploaded manuscript actually inspected](https://www.researchgate.net/publication/257789723_Visualization_of_Four_Normal_Size_Limit_Cycles_in_Two-Dimensional_Polynomial_Quadratic_System).

The original author PDF endpoint returned an HTML access error in this
review; it was not a retrieved PDF. The accessible manuscript was read.
Yu–Zeng also numerically visualize this decimal field using an ordinary ODE
solver; this is independent published numerical support, not an interval
certificate. [Yu–Zeng, equations (8), (17) and Fig. 5](https://arxiv.org/html/2002.09987).

The exact incumbent rational tuple is
`(-10,11/5,7/10,-363889/5000,3/2000)`. Its beta-zero shape has
`K=-674997/250000<0`; the degeneracy at this c is `alpha=-700/9`.
The intended `K>0` precursor is on the other side of that degeneracy.

The recovered numerical record is now specific:

| Field | Saved numerical evidence |
|---|---|
| Incumbent, positive beta | Origin returns near 0.683210217, 2.183699825, 15.962783982; remote return near -3711.560806; S/U/S plus U |
| Same incumbent shape, beta=0 | Two detected nonzero origin roots near 3.068845425 and 15.064071451; U/S |
| New starting shape `(c,alpha,beta)=(7/10,-80,0)` | One detected origin S root near 64.555434341, multiplier 0.809691136; separate remote U control near -5391.141160, multiplier 12.168019 |

Here the new starting shape has exact `K=6/5>0`. The root controls are
ordinary floating-point DOP853 calculations, explicitly labelled
nonrigorous. The relevant raw files are `data/kkl.json`,
`data/hopf_probe.json`, and `data/hopf_start_remote.json`. They are not the
high-precision GT controls elsewhere in the packet.

The handoff describes loss of the inner incumbent cycle as the negative
control. Its saved `hopf_probe.py` computes twelve static beta-zero shapes;
it is not a recorded beta-continuation trajectory. The two detected roots
support the stated loss and agree with the local KKL mechanism, but do not
certify that exactly two finite origin cycles exist, or locate every event
along the path from beta=0.0015. Do not turn that evidence into a certified
branch-loss theorem or invent a measured fold location.

## 3. What the fold must accomplish

The coherent count starts from the new **one-origin-cycle** control, not
from pretending that the incumbent has four finite cycles at beta=0.
An ordinary fold can add a stable/unstable pair, giving three origin cycles
while the old stable orbit and remote unstable orbit persist. A fold creates
or removes two cycles. Reaching a fold that only destroys the tracked stable
cycle does not establish this coexistence event.

For a valid first-return map write `D(r;c,alpha)=P(r;c,alpha)-r`. A generic
fold requires

    D=0, D_r=0, D_rr != 0,
    a nonzero derivative of D in a transverse shape direction.

In `(r,c,alpha)` these two equations normally define a curve. On the
two-dimensional shape plane its image is a codimension-one fold curve.
Viewed in the restricted three-parameter family `(c,alpha,beta)`, requiring
both the Hopf plane and the fold is generically codimension two. The final
precursor lies off the fold, on its two-cycle side. Neither this dimension
count nor availability of two shape parameters proves that such a fold
exists in the box or that its new pair coexists with the other two cycles.

A continuation of the seeded root sheet is a defensible bounded experiment.
It cannot exclude disconnected root sheets or isolas that it never reaches.
Continuation must keep a cycle ledger at common parameters and orient each
actual fold. An S/U/S triple on one valid return itinerary, with all three
derivatives separated from one and a simultaneous remote U root, is the
first credible precursor signal.

## 4. Geometry and certification cannot be deferred

Use the handoff's exact equilibrium polynomial and Jacobian gates before
spending returns: unique real remote root left of x=-1, stable remote focus,
no multiple finite equilibrium, and strict K margin. The line x=-1 has
`x'=1` and cannot be crossed by a periodic orbit. This separates the nests
but does not certify a selected numerical path is a full return.

The c boundaries `241/250` and `1` change the infinity portrait. A branch
reaching either boundary needs a compactified itinerary account; no-return
at a finite integration horizon is not a cycle-disappearance proof. The
existing sparse remote scans contain changed return itineraries and cannot
be summed as additional cycles.

For the downward y=0 section the actual return derivative is

    P'(r) = Q(r,0)/Q(P(r),0) * exp(integral div(F) dt).

The divergence exponential alone is the multiplier at a periodic fixed
point. Using it as P' away from a root corrupts fold location and derivative
certification. The existing probe does not implement a fold solver and
must not be described as having already followed fold curves.

The handoff's section ranges are positive `r in [2^-12,2^10]` and remote
`r in [-2^20,-1]`, with remote `r<alpha/10`. They bound the experiment,
not all possible cycles. Prove full-return existence on each entire section
interval, fixed itinerary and orientation, strict endpoint displacement
signs, and a derivative enclosure excluding one or interval-Newton
inclusion. Order the origin sections and prove their orbits have the same
nesting assignment; prove separation from the remote tube and the future
Hopf neighborhood.

Only after the four finite gates coexist should one choose rational c and
alpha in their common open neighborhood and a sufficiently small explicit
negative rational beta. Revalidate all four old cycles at that same beta,
and certify the fifth by a Hopf remainder bound or validated small-section
shooting. The leading radius formula `r^2 ~ 4 alpha beta/K` is a locator,
not a fifth-cycle certificate. No global exact-count proof is needed to
establish at least five.

## 5. Fast-kill budget and honest cost basis

The inherited bounds are 256 continuation steps per seeded branch,
4096 total return/derivative evaluations, and 64 adaptive parameter cells
around actual fold or root-count events. Keep these as hard experimental
ceilings; count failed returns and derivative calls, save unresolved events,
and do not expand b or the box after exhaustion. Exhaustion without a
precursor means no precursor found on the explored sheets, not a theorem
about this whole box or H(2).

Before the full allowance, a small implementation gate should establish
that the corrected derivative and a bounded continuation preserve the
known K>0 origin control and its remote return on a common neighborhood.
A suggested pilot cap is 64 total evaluations or 16 accepted continuation
steps, whichever is reached first. This is a proposed engineering stop,
not a mathematical exclusion test. Stop early on an invalid return map,
wrong multiplier convention, untracked remote loss, or a shape leaving the
exact gates. Lack of a fold in this pilot is not an exclusion.

No KKL or Hopf-probe elapsed runtime is saved in the inspected JSON files.
The saved GT MPFR replay times, approximately 0.24–0.26 seconds per
order-112 endpoint and 0.34–0.35 seconds at order 128, are different fields
and algorithms and cannot price KKL continuation. The packet also reports
historical Galias–Tucker validated existence times of 4–14 seconds and
derivative/uniqueness times of 45 seconds–40 minutes per cycle on their
hardware. Those are certification benchmarks, not KKL estimates.

Therefore record actual wall and CPU time, successful/failed evaluations,
and solver work in the pilot, then estimate the remaining capped work from
that measured distribution. Price long nonreturns separately with explicit
per-call limits. A cost promise obtained by multiplying 4096 by a GT tiny-
cycle timing would be unsupported. The repository has numerical scripts,
not a completed interval verifier; its implementation and final five-cycle
certificate are separate conditional costs.

## 6. Council verdict

The first-order Hopf completion is a correct conditional construction and
the recovered K>0 one-cycle control makes its missing-pair fold experiment
well specified. The global coexistence premise remains open. The available
record contains no three-nonzero-origin-cycle K>0 precursor, no continuation
of its fold curves, and no new five-cycle certificate.

The least wasteful next commitment is the bounded common-parameter
return-sheet/fold pilot with an explicit remote-persistence ledger. Neither
generic codimension language, the original KKL four-cycle visualization,
nor a lost numerical return substitutes for that first new piece of
evidence. This council review ran no attack and makes no new construction
or exclusion claim.

## 7. Corrections to the subsequently received Fable council

The original `FASTRA_COUNCIL_2026_09_04.md` at the received main commit
`2136896` is historical council input. Its KKL statements need these
corrections in a separate addendum:

1. Section 6(a) says to seek the missing fold where there are already two
   origin cycles, not one. Restore the actual handoff: start from the K>0
   **one**-origin-cycle sheet; a new fold pair can give **three**. Keep the
   existing stable origin orbit and remote orbit in the coexistence ledger.
2. Section 1 calls the KKL replay high-precision numerical. It used ordinary
   double-precision DOP853. The separate MPFR GT replay does not change the
   KKL evidence level. Also replace the stated Shi/GT cycle scales: numbers
   such as 10^-13 and 10^-200 are coefficient scales. The exact GT section
   radii are approximately 7.07e-75, 2.25e-21 and 6.67e-8.
3. Sections 4 and 8 propose an every-sheet box conclusion from fold
   continuation. Bounded seeded continuation gives no complete cover of
   disconnected sheets, omitted initial radii or unresolved itineraries.
   Its negative output must describe only the explored branches. A global
   or box-wide bound needs a separate coverage proof.
4. Section 3.1(1) says a putative order-one weak-focus bound would settle
   H(2)=4 for two-focus systems. The stated bound alone concerns a Hopf
   stratum. Extending it to all strong-focus fields needs an additional
   deformation or termination theorem preserving the required finite
   cycles. That global implication is not established in this council.
   The `3-k` pattern is motivation for a conjecture, not a kill theorem.
5. “Finite size” does not itself prove good conditioning. The incumbent's
   first two multipliers are close to one and the remote orbit reaches
   large coordinates. “4096 evaluations, days” is a planning heuristic;
   no measured KKL runtime in the packet supports that wall-time estimate.
   The next-hour objective can attempt bounded continuation and record its
   first event; it cannot promise to find a fold within that hour.

These corrections leave the conditional Hopf completion intact and do not
assert that its precursor is impossible.
