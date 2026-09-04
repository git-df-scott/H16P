# Thought session: how would Euler or Newton hunt a fifth quadratic limit cycle?

2026-09-04. Written after Attack 1 was closed by theorem and Attack 2 was
found to be topologically defective. This is reasoning, not a result. Where
a statement is a computation done in this repository it says so; where it is
a literature fact it names the source; where it is a conjecture it says so.

## 1. What Newton would do first: find the one quantity that decides

Newton did not search; he found the single function whose zeros are the
phenomenon and then studied that function's shape. Here the function is
already known. Fix a quadratic field with a focus `F_1` and a transversal
ray from it. The displacement of the first return map,

\[
 d(r)=\Pi(r)-r,
\]

has one zero per limit cycle around `F_1`. Everything else is bookkeeping:

- A cycle around `F_1` cannot enclose any other singular point (Coppel), so
  the whole nest of `F_1` is one interval `(0,r_{\max})` of the ray, and
  `r_{\max}` is where the ray meets the boundary graphic `\Gamma` of the nest.
- Near `r=0`, `d(r)=\lambda r+\eta_1r^3+\eta_2r^5+\eta_3r^7+\dots`, and a
  quadratic focus has order at most three (Bautin). So at most three zeros can
  be pushed out of the focus, and only with the alternating hierarchy
  `\lambda,\eta_1,\eta_2,\eta_3`.
- Near `r=r_{\max}`, `d` has a Dulac-type expansion in the coordinate
  `w=r_{\max}-r` dictated by the type of `\Gamma` (finite saddle loop,
  graphic through infinite singularities, or a periodic orbit).
- Two foci at most; with two nests one of them holds at most one cycle
  (Zhang Pingguang, QTDS 2002). So **a fifth cycle means four zeros of `d`
  in one nest**, or five in a single nest.

That is the whole problem: make `d` change sign four times on `(0,r_{\max})`.
Newton would now ask what controls the sign of `d` at the two ends and how
many times a function with the known constraints can cross zero in between.

## 2. What Euler would do first: count parameters, then compute

Euler trusted counting and computed relentlessly. After affine changes and
time scaling a quadratic field with a focus at the origin has **five**
essential parameters (the Shi chart `\lambda,l,m,a,b` is one such slice).
A limit periodic set of cyclicity `k` needs, generically, a `k`-parameter
unfolding to realise `k` cycles. So a configuration producing five cycles
from one degenerate object must sit at a **codimension-five point**, an
isolated point of the five-dimensional parameter space. Every mechanism can
now be priced:

| Degenerate object | Codimension | Max cycles it can emit | Source |
|---|---:|---:|---|
| Weak focus of order `k` | `k` | `k` (`k\le3`) | Bautin |
| Hyperbolic saddle loop, saddle quantity `\ne0` | 1 | 1 | Leontovich |
| Saddle loop with zero saddle quantity | 2 | 2 (3 with further degeneracy) | Roussarie 1986; Han 1997 for integrable QS: 2 except one case |
| Semistable (double) limit cycle | 1 | 2 | classical |
| Center, codim-3 components `Q_3^H, Q_3^R, Q_3^{LV}` | 3 | open annulus: 2 for `Q_3^H` (Horozov–Iliev, Gavrilov), `Q_3^{LV}` generic 2 (Żołądek); **`Q_3^R` generic: almost nothing known** (Gavrilov, arXiv:1610.07582) | literature |
| Center `Q_4`, codim 4 | 4 | open annulus: **`\le4` distinct, `\le3` on the lobe region** (this repository, Theorem N) | Q4_THEOREM_N.md |
| Graphics through infinity, 2-saddle cycles | 1–3 | mostly `\le2` where proved (DRR program, 121 graphics) | literature |

Addition rule: two degenerate objects in the same nest can add their counts
only if they are *simultaneously* degenerate at one parameter point and
their unfolding directions are independent. Codimensions then add. The
only ways to reach five with five parameters:

1. order-3 focus (3) + boundary graphic of cyclicity 2 (2), same nest;
2. order-2 focus (2) + boundary object of cyclicity 3, or + two independent
   cyclicity-1 objects (loop and double cycle) in the same nest;
3. a codim-3 center whose *closed* annulus (loop included) has cyclicity 4,
   plus one cycle in the other nest (`4+1`);
4. a codim-4 center (`Q_4`) whose closed annulus has cyclicity 4, plus one in
   the other nest;
5. a codim-5 graphic of cyclicity 5 (none is known).

## 3. What the theorems and this repository's computations have already killed

- **Route 4 is dead.** `Q_4`'s open annulus carries at most four distinct
  zeros (Theorem N with (N1)), the loop adds nothing (Lane B), and the second
  finite singularity of every `Q_4` center is a strong **node**, not a focus
  (`audit/claude_center_identify.py`: trace `5, 8.4, 15.5` at
  `kappa=1.25, 2, 5`), so there is no second nest to put a fifth cycle in.
- **Route 1 as "focus + finite saddle loop" is dead on the order-3 stratum.**
  In the Shi chart the finite saddles exist only for `3a^2>l^2+2l`, and in
  that region the saddle loop around the origin occurs exactly on the curve
  `\eta_3=0`, to `10^{-14}` at four values of `a`
  (`audit/claude_laneC_splitting4.py`). There the origin is a **reversible
  center**, not a weak focus. Li Chengzhi (1986) and Cherkas proved no limit
  cycle surrounds an order-3 weak focus in any quadratic system; the
  computation says the loop, as a limit of cycles, is excluded as well.
- **Route 1 as "focus + graphic through infinity" is constrained by the same
  theorem.** Li–Cherkas forces the boundary graphic of the origin nest to
  have, at every point of the order-3 stratum, the same stability as the
  focus (otherwise Poincaré–Bendixson would give a cycle). After unfolding
  the focus, `d` has three small zeros and then the sign of `\eta_3`, which
  equals the sign at the boundary: no fourth crossing unless the boundary
  emits a **pair**. A pair requires the boundary graphic to have cyclicity
  `\ge2` at the stratum point. That is route 1 in its only surviving form.
- **Route 3 with `Q_3^H` or `Q_3^{LV}`** is dead by the known bounds (2). With
  `Q_3^R` it is **open**: the literature explicitly says almost nothing is
  known about the generic reversible case, and the conjecture (cyclicity 2,
  or 3 in special cases) is a conjecture.

## 4. Where a counterexample can still live

Two concrete places remain, and both are computable in the Euler sense.

### 4a. Reversible centers `Q_3^R` with a saddle loop and a second focus

The loop point found on the stratum is exactly such a center: origin a
reversible center, a zero-trace saddle on the symmetry axis with a
homoclinic loop bounding the annulus, a second antisaddle at `(0,1)`, and a
fourth finite equilibrium. A quadratic perturbation of this center produces
cycles in the annulus (first non-vanishing Melnikov function; first order
may vanish identically on reversible perturbations, which is why higher
Melnikov functions and the Bautin ideal enter) and possibly one cycle around
`(0,1)`. A fifth cycle needs the **closed** annulus of some reversible
center to have cyclicity four. Nobody has excluded this. Nobody has
computed it for the whole two-parameter family either. This is the direct
analogue of what Astra did for `Q_4`, and the tools (Picard–Fuchs, Stieltjes,
Chebyshev, the `\Phi`-functional trick) transfer, with the extra work of
higher-order Melnikov functions. Prediction from the way `Q_4` went: the
answer is probably three, but *probably* is not a theorem, and the
reversible family is the largest unexplored integrable component.

### 4b. The order-3 stratum's boundary graphics of cyclicity two

Llibre–Schlomiuk (Canad. J. Math. 2004) classify quadratic systems with a
third-order weak focus into 18 phase portraits, some with graphics. For each
stratum point the origin-nest boundary `\Gamma(l,a)` is a graphic through
infinite singular points (there is no finite saddle in the Shi box; in the
finite-saddle region the loop only appears at the center curve). The
questions are, in order:

1. What is `\Gamma(l,a)` for each of the 18 portraits, and which of its
   singular points are hyperbolic or semi-hyperbolic?
2. What is its first stability coefficient (the product of hyperbolicity
   ratios of the infinite saddles, adjusted by the regular transition)? Its
   sign must agree with `\eta_3` everywhere by Li–Cherkas; where it is
   *neutral* (equal to one) the graphic is a candidate for cyclicity two.
3. Is that neutrality curve nonempty on the stratum, and is the second
   coefficient there of the sign that lets a pair of cycles be born inward?
4. Are the two unfolding directions (graphic breaking and graphic stability)
   independent of the three focus unfolding directions `\lambda,\delta,\epsilon`?

If 3 and 4 are yes, there is a codimension-five point whose unfolding gives
`3+2=5` cycles around one focus. If the neutrality curve is empty, route 1
is dead and the order-3 stratum can only ever give three plus one.

## 5. What Euler and Newton would *not* do

- They would not run a random coefficient sweep. The count in section 2
  says five cycles are a codimension-five event; a sweep in five dimensions
  finds codimension-zero phenomena.
- They would not trust a plot of five long transients; Newton's whole
  method was to replace a picture with a function and a series.
- They would not stop at "the conjecture is probably true". Euler computed
  the sums nobody believed could be summed. The honest state is that
  `H(2)=4` rests on: Bautin, Li–Cherkas, Zhang, the `Q_3^H`/`Q_3^{LV}`
  bounds, Theorem N, and a large but incomplete DRR program. The gaps are
  exactly 4a and 4b.

## 6. Concrete next computations (Claude, own lanes)

1. **4b, step 1–2 on the Shi stratum.** For `(l,a)` on `\eta_1=\eta_2=0`,
   compactify, find the infinite singular points, identify `\Gamma(l,a)` and
   compute its stability coefficient as a function on the stratum; locate its
   neutrality set. Cost: elementary numerics plus the DRR formulas for
   hyperbolicity ratios at infinity.
2. **4a, step 1.** Write the reversible center at the loop point in a
   standard `Q_3^R` normal form, compute its period annulus, loop, and the
   first two Melnikov functions of a general quadratic perturbation; check
   whether the first is identically zero on the reversible subfamily and
   what the first non-vanishing one's dimension is. Then run the `Q_4`
   playbook: Chebyshev bound, `\Phi`-type necessary conditions, loop
   coefficients.

Astra's parallel assignment stays as written in FASTRA_H16_HANDOFF_5.md
(the outside-lobe `Q_4` four-zero question, which is now a pure-mathematics
completion, not a counterexample route).
