# Next strike: preserve four cycles and seek an additional outer pair

**Status: proposed experiment, not executed. No counterexample found.**

The question is whether ONE real planar polynomial field of total degree at most two has at least five distinct isolated periodic orbits. Five or six rigorously verified cycles would meet that target. Sampled return brackets, center ovals, auxiliary integral zeros, and cycles from different parameter values do not.

## Why this experiment

The overnight Q4 continuation and extreme-scale reversible searches repeatedly produced three numerical crossings. Increasing radius alone did not add cycles. The next experiment should start with the saved four-cycle numerical positive control and explicitly preserve its three upper and one lower cycles while seeking a new outer fold in the upper annulus. A nondegenerate fold could create two additional upper cycles, giving six in one field. This is a testable construction hypothesis, not evidence that such a fold exists.

Source of the seed: [original-field control](reversible_reseed/verify_control.py) and its [data directory](reversible_reseed/data/). Its four-cycle evidence is numerical, not interval certified. Reproduce the finite-field returns before treating this as a working seed.

Use the family

```
P = (b-2)/4 + (1-b)y + a*x^2 + b*y^2 + e1*x + e2*x*y
Q = e0 - 2*x*y
```

Start at a=-7/4, b=1/3, e0=-4*tau/11,
e1=-31379*tau/25000, e2=-7517*tau/5000, with tau=1e-4.
The independent control also used tau=5e-5. Each is a separate field; their cycles must never be added together. Allow all five coefficients (a,b,e0,e1,e2) to vary during continuation. This family is a restricted subset of quadratic systems, so failure would not settle H16P.

## Bounded execution sequence

1. **Re-establish the four-cycle seed.** Recompute the six saved section witnesses using the Cartesian solver and the logarithmic solver. Refine four disjoint brackets, identify the complete return itinerary, and estimate the derivative at each root. Stop if either solver fails to reproduce the four brackets or the complete returns. Do not start a fold search on a failed control.

2. **Track the four existing cycles while moving the coefficients.** Use predictor-corrector continuation with event-corrected parameter and section derivatives. Maintain disjoint isolating section intervals for the three upper roots and the lower root. A branch is rejected if a root is lost, intervals merge, an equilibrium changes the itinerary, or a return gate fails. Keep e0 nonzero; the invariant-line center limit is not an acceptable candidate.

3. **Solve for an additional upper outer double zero.** Beyond the outermost tracked upper cycle and within the same valid return domain, define D(s;theta) as the forward/backward log-height mismatch at section coordinate s=log(T), theta=(a,b,e0,e1,e2). Seek D=0 and partial_s D=0. An initial formulation fixes the four tracked root positions and solves their four mismatch equations together with these two fold equations for the five coefficients and s. This six-equation formulation is a local search device: fixed positions can miss solutions. If its Jacobian is ill-conditioned, release one anchor position and continue the resulting curve while retaining all four brackets.

4. **Test the putative unfolding in actual fields.** Require a nonzero second section derivative and a nonzero transverse parameter derivative, checked independently. Move to both sides of the candidate fold while tracking the four old roots. The productive side must show two additional disjoint upper brackets in the SAME coefficient vector. A small optimizer residual alone is not a result.

5. **Freeze and independently verify any candidate with at least five brackets.** Save exact rational coefficients or exact decimal coefficients interpreted as rationals. Replay complete returns with an independent arbitrary-precision implementation. Then use validated integration and interval return-map/derivative bounds to establish existence and isolation on disjoint sections. Check that brackets correspond to different geometric cycles, with no equilibria or itinerary changes in the certified domains. This certification machinery still needs implementation; tighter floating-point tolerances do not substitute for it.

## Cheap falsification before expensive runs

- Use displacement normalization fixed from the seed and a documented local scale. Do not let a denominator grow to manufacture a small objective. Report raw displacement, scaled displacement, derivatives, and precision sensitivity together.
- Use variational derivatives checked against independent differences; event timing must be included.
- Reject an almost-center field whose entire displacement vanishes numerically. Require distinguishable, isolated tracked roots.
- Record every failed return. An unresolved end of the section is a missing domain, not a zero or evidence of absence.
- Do not continue the exact third-order weak-focus-plus-surrounding-cycle route. The [primary-source correction](research_2026_09_08/outputs/H16P-focus-route-correction.md) explains why that precursor was abandoned.
- The prior [compact-orbit compatibility argument](research_2026_09_08/outputs/H16P-center-locus-compatibility.md) applies to a particular boundary-canceling direction. Do not claim it excludes this entire finite-parameter experiment, or ignore it when a branch approaches that direction.

## Stop rules and useful negative outcome

First pass: one independently reproduced seed and at most two locally distinct continuation branches, with a predeclared cap of 50 accepted continuation steps per branch. Stop early on failed controls, loss of an existing cycle, persistent domain failure, or residuals that do not survive increased precision. Such a cap bounds the experiment; it is not an exhaustive parameter search.

If neither branch produces a credible fold, publish the explored coefficient ranges, surviving root intervals, closest fold residuals, and failure reasons. Do not spend the remaining budget merely extending radii. The useful negative result would be a documented obstruction along these branches, not a global nonexistence theorem.

## Evidence ledger

| Claim | Evidence | Status / next falsification |
|---|---|---|
| The seed has a 3+1 cycle pattern | Saved finite-field control and moment witnesses | Numerical; reproduce with independent full itinerary checks |
| Current extreme-radius searches do not establish five cycles | [Extreme boundary report](research_2026_09_08/outputs/H16P-extreme-boundary-search.md) | Saved results; incomplete return regions remain |
| An outer fold could add two cycles while four persist | Local fold mechanism under nondegeneracy and persistence assumptions | Construction hypothesis; solve and test the unfolding |
| This proposed strike will succeed tonight | No supporting evidence | No such prediction is justified |

## Budget and execution status

The user assigned a 30% usage budget. The shared main meter was 6% when it was recorded and is 33% when this plan is written: a 27-point account-wide increase, not measured task-only consumption. Earlier work stopped at 29% total under a more conservative interpretation. The remaining allowance cannot be measured exactly from that shared meter. This document and its push do not start the proposed numerical search or resume the paused automation. Any execution must remain bounded by the user's budget and retain capacity for verification and reporting.

The most valuable next computation is Step 1, followed by the constrained outer-fold test. There is no honest guarantee of a counterexample tonight.
