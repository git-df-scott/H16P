# H16P counterexample hunt — executable plan

Date: 9 September 2026. Planning only; no new ODE search or certification was run.

**Recommendation:** reproduce the saved four-cycle field, then continue its four cycles while seeking a separate outer fold in the upper nest. Use moving cycle coordinates and explicit return-domain checks. Reserve the endpoint-asymptotic route for a specific obstruction identified by this experiment.

The target remains one real planar polynomial vector field of degree at most two with at least five distinct isolated periodic orbits. There is no counterexample in the records reviewed, and this plan offers no guarantee of one.

## Repository basis and precedence

Pinned main: [`498b58fc1378486601f348f67e0ab3a5362e7933`](https://github.com/git-df-scott/H16P/commit/498b58fc1378486601f348f67e0ab3a5362e7933), published 9 September. I inspected the branch list, main tree, recent commit history, current next-strike document, September 8 status/corrections and selected supporting reports, the seed implementation and its complete saved JSON, the logarithmic return implementations, and the certification protocol. This was a targeted planning review, not a full repository or theorem audit. Separate branch histories were not comprehensively reread.

Source precedence matters: the root README and STATUS still contain September 5 checkpoints and an older blanket Q4 closure. The newer archive explicitly rejects extending interior results to closure of the original boundary problem. Its corrections take precedence. The current next-strike proposal is explicitly **unexecuted**.

| Evidence from the repository | Consequence for the hunt |
|---|---|
| Saved reversible control has a numerical 3+1 sign pattern at each of two separate parameter points | Start from one of these fields and reproduce it; never combine their cycles |
| Extreme-scale searches reached a lower bracket near 10^249, with three brackets per field | Larger radius alone is not a productive objective |
| Six finite-Q4 outer-fold fits retained nonzero normalized residuals | Do not repeat those same anchor/modulus searches and call convergence a fold |
| Exact higher-order-focus precursor routes were retracted using published restrictions | Do not restart those prohibited exact-focus constructions |
| Center-divided compatibility argument applies to a particular direction on compact regular ovals | Treat that direction as constrained; do not infer a uniform endpoint or whole-family exclusion |
| The new four-cycle-preserving outer-fold search has not been executed | This is the immediate experiment worth resolving |

## 1. Recover one trustworthy starting field

Use the five-parameter family from the current handoff:

\[
\dot x=(b-2)/4+(1-b)y+ax^2+by^2+e_1x+e_2xy,
\qquad \dot y=e_0-2xy.
\]

At the starting point set

\[
a=-7/4,\quad b=1/3,\quad e_0=-1/27500,\quad
e_1=-31379/250000000,\quad e_2=-7517/50000000.
\]

These are the exact coefficients of the recorded arc at tau=10^-4. The tau=5×10^-5 point is a separate secondary control, not an additional source of cycles.

The following are saved Cartesian half-return mismatches, not freshly reproduced values and not interval bounds:

| Nest | Starting section heights y at x=0 | Saved signs of forward_y − backward_y |
|---|---|---|
| Upper | 1.65072993, 3.48477147, 7.31241834, 28.17882337 | +, −, +, − |
| Lower | −3161.39476, −20733.33317 | +, − |

Two upper witnesses have mismatch magnitudes approximately 2.06×10^-9 and 2.31×10^-9. Parse coefficients and full witness strings accurately; use the abbreviated table only for orientation. The JSON's 65-digit precision describes its moment quadrature, while its finite-field shooting uses floating-point SciPy. Do not present it as 65-digit ODE evidence.

Reproduce the six witnesses with Cartesian integration and the latest logarithmic implementation, then refine four disjoint root brackets and compute complete first-return itineraries and multipliers. Cross-check displacements after accounting for the coordinate transformation: log absolute height preserves the upper mismatch sign but reverses the Cartesian mismatch sign on the lower section.

**Pass gate:** both coordinate implementations agree on four brackets and valid itineraries; higher precision supports the smallest signs and simple-root derivatives. Failure stops continuation until explained. Passing remains numerical evidence unless validated integration is used.

## 2. Build the derivatives needed for a meaningful search

Use a log section coordinate s=log|y| and mismatch D between forward and backward half passages. Establish the local correspondence between D=0 and a periodic orbit; retain a full first-return map for isolation and Floquet checks.

Integrate state, section, and parameter variations. For an endpoint on x=0, include the event correction

\[
\delta y_{\mathrm{hit}}=\delta y-(Q/P)\delta x,
\]

then apply the log-height chain rule. Check transversality before dividing by P. Initialization kicks and any parameter-dependent section positions also require consistent derivatives. Check selected derivatives against independent high-precision differences at decreasing step sizes.

Keep the logarithmic solver's actual gates in view: its height thresholds 1/2 and (2-b)/(2b) are convenient section partitions, not the perturbed equilibrium locations when e0 is nonzero. Track actual equilibria, section tangencies and crossing order. If the chosen section ceases to work, establish a valid replacement; a failed gate is unresolved geometry, not cycle disappearance.

Normalize residuals using documented seed-based scales. Always retain raw D, scaled D, derivatives, integration tolerances and precision sensitivity. Never obtain a small objective by driving the perturbation to a center or enlarging a denominator.

## 3. Continue moving cycles and locate a genuinely new outer extremum

Let theta=(a,b,e0,e1,e2), with tracked root coordinates s1<s2<s3 in the upper nest and s4 in the lower nest. Continue the four equations D_i(si;theta)=0 while allowing the roots to move. These equations define a locally five-dimensional set when all four roots are simple; they do not by themselves choose a path. Choose a direction in scaled coefficient coordinates, solve for the root motion, and use trust-region prediction/correction with the explicit gates below.

The inherited proposal initially fixes all four root positions. That is a legitimate narrow diagnostic, but it can miss folds whose approach requires root motion. Use moving roots from the start; fixed anchors may supply a comparison, not a coverage claim.

Search outside s3, within a verified numerical return domain, for a separate stationary point s* satisfying D_s(s*;theta)=0. Then continue this stationary point while driving D(s*;theta) toward zero. If no outer stationary point is initially present, use a bounded adaptive derivative profile on at most two selected coefficient paths. Do not launch an unconstrained double-zero solve from an arbitrary large radius.

Three existing roots already require stationary points between them by Rolle's theorem. Their presence is baseline geometry. The useful signal is a separate outer extremum approaching zero while the existing four roots persist.

An augmented continuation formulation is

\[
D_i(s_i;\theta)=0\ (i=1,2,3,4),\qquad D_s(s_*;\theta)=0.
\]

There are ten unknowns (five coefficients, four roots, one stationary point) and five equations. A one-dimensional predictor-corrector path therefore needs four additional explicitly recorded slice constraints plus a pseudo-arclength condition. This dimension accounting prevents treating an underdetermined system as a uniquely specified search. Choose the slices using scaled sensitivities and compare at most two locally distinct choices.

Accept a step only while all four roots remain simple and separated, their full return itineraries remain valid, and the outer stationary point stays distinct from the old roots. Keep e0 nonzero and record approach to center or invariant-line limits. Reduce the step on a threatened gate. If it cannot be resolved within the local budget, retain the last valid state and classify the obstruction.

## 4. Unfold the fold in one actual field

A candidate must satisfy

\[
D(s_*;\theta_*)=D_s(s_*;\theta_*)=0,
\quad D_{ss}(s_*;\theta_*)\ne0,
\quad D_\lambda(s_*;\theta_*)\ne0
\]

for a transverse parameter direction lambda permitted by the family. Verify the last condition on the admissible continuation set, not merely for a parameter direction that violates imposed constraints.

Test the two sides of the fold using complete return maps. Four persistent old cycles plus a new nondegenerate pair would give six; six also refutes H(2)=4. We are not assuming such a fold exists. At the fold itself an isolated double cycle together with the old four could already give five, but certifying a nearby hyperbolic pair is the preferred numerical route.

Stop discovery as soon as one fixed coefficient vector supports at least five distinct credible brackets. Record the entire vector for every count; never aggregate across continuation steps or use a remote cycle from a different field. If an old cycle is lost before the fold, the preserve-four branch has failed its stated gate. Save the event; do not quietly reinterpret it as a successful continuation.

## 5. Certify the candidate rather than extend the search

Freeze exact rational coefficients and independently replay the candidate with arbitrary precision. Rationalization must be followed by a fresh replay of that exact field. Seek five disjoint validated isolating regions and interval Poincare maps proving existence and isolation. For hyperbolic roots, require strict endpoint signs plus a derivative interval excluding zero for P−id, or a validated interval-Newton inclusion.

Prove complete oriented first returns, absence of equilibria in the relevant flow boxes, and geometric distinctness; several section intersections may belong to one orbit. Archive coefficient encodings, return-time and derivative enclosures, dependencies, logs, and replay instructions. Exactly counting all cycles is unnecessary.

The main handoff says this candidate-certification machinery still needs implementation. Any separate branch that claims a ready verifier must be inspected and reproduced before reuse. Tighter floating-point tolerances alone are not a certificate. Begin with a seed-cycle validation smoke test before trusting the pipeline with a new candidate.

## 6. Bounds, outcomes and fallback

First pass: one primary seed, at most two locally distinct continuation branches, and at most 50 accepted steps per branch, as in the current handoff. Also impose an attempted-evaluation cap, because accepted-step limits do not bound repeated Newton failures. A proposed initial cap is 500 paired mismatch evaluations per branch, including rejected trials and root refinements; each pair costs two half passages. Derivative integrations, solver function evaluations, arbitrary-precision replays and wall time must also be logged. This is a planning cap, not an assertion about current available account usage. Historical account percentages and paused automations are not a new allowance or execution instruction.

Reserve capacity for replay, certification and reporting. Stop a branch early when signs fail precision checks, a center limit destroys isolation, an old cycle cannot be retained, or no reliable return domain can be established. Any negative outcome applies only to the explored paths and bounds.

Deliverables for execution: a seed replay ledger; moving-root and event histories; raw/scaled outer-extremum residuals; exact coefficient snapshots; explicit unexplored gaps; and either a candidate certification package or a scoped negative report.

If both branches end without a credible fold, the next task is analytic: derive the parameter-dependent endpoint return map at the observed limiting organizer, including the compact-orbit displacement, the boundary terms and a uniform remainder in one stated joint scaling. The key question is whether the same coefficients can preserve the known roots and create the missing boundary roots. The compact compatibility argument does not answer that endpoint question. Do not restart extreme-radius sweeps without such a new balance and a testable sign prediction.

Marín–Villadelprat explicitly distinguish individual and simultaneous cyclicity; their hemicycle analysis is relevant background for that fallback. It does not itself supply five cycles here. [Primary paper](https://arxiv.org/html/2501.16924v1).

## Source links

- [Current next-strike proposal](https://github.com/git-df-scott/H16P/blob/498b58fc1378486601f348f67e0ab3a5362e7933/H16P_NEXT_STRIKE.md)
- [Latest archived status](https://github.com/git-df-scott/H16P/blob/498b58fc1378486601f348f67e0ab3a5362e7933/research_2026_09_08/outputs/H16P-current-status.md)
- [Control implementation](https://github.com/git-df-scott/H16P/blob/498b58fc1378486601f348f67e0ab3a5362e7933/reversible_reseed/verify_control.py) and [saved data](https://github.com/git-df-scott/H16P/blob/498b58fc1378486601f348f67e0ab3a5362e7933/reversible_reseed/data/verified_control.json)
- [Latest logarithmic return solver](https://github.com/git-df-scott/H16P/blob/498b58fc1378486601f348f67e0ab3a5362e7933/research_2026_09_08/outputs/reversible_log_return.py)
- [Finite-Q4 failed fold searches](https://github.com/git-df-scott/H16P/blob/498b58fc1378486601f348f67e0ab3a5362e7933/research_2026_09_08/outputs/H16P-Q4-outer-fold-search.md)
- [Focus-route correction](https://github.com/git-df-scott/H16P/blob/498b58fc1378486601f348f67e0ab3a5362e7933/research_2026_09_08/outputs/H16P-focus-route-correction.md)
- [Compact center-locus compatibility](https://github.com/git-df-scott/H16P/blob/498b58fc1378486601f348f67e0ab3a5362e7933/research_2026_09_08/outputs/H16P-center-locus-compatibility.md)
- [Extreme-scale report](https://github.com/git-df-scott/H16P/blob/498b58fc1378486601f348f67e0ab3a5362e7933/research_2026_09_08/outputs/H16P-extreme-boundary-search.md)
- [Implementation audit](https://github.com/git-df-scott/H16P/blob/498b58fc1378486601f348f67e0ab3a5362e7933/research_2026_09_08/outputs/H16P-implementation-audit.md)
- [Certification protocol](https://github.com/git-df-scott/H16P/blob/498b58fc1378486601f348f67e0ab3a5362e7933/RIGOROUS_CERTIFICATION.md)
