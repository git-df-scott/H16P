# Five-cycle mechanisms after Q4

Audit cutoff: September 4, 2026. This file is the attack map; detailed theorem hypotheses and independent algebra are in the linked appendices.

## Configuration first

Zegeling 2024, [Theorem 1.2](https://doi.org/10.1515/anona-2024-0012), permits only (n,0) or (n,1). Therefore the five-cycle target is **(5,0) or (4,1)**. Configurations (3,2), (3,1,1), (2,2,1) are excluded. Four distinct cycles in (2,2) are also excluded. No accepted (4,0) seed was located. At four distinct real finite equilibria, Theorem 5.4 further permits only (n,0) or (1,1).

A cycle surrounds exactly one focus; a cycle surrounding an entire collection of foci and saddles is not an admissible quadratic architecture. A graphics cyclicity theorem and a focus cyclicity theorem cannot be added unless both occur in the same quadratic unfolding and all counted cycles persist simultaneously.

## Scores, 0–5

Higher is favorable in every column. Open-status/novelty measures an identifiable unresolved problem, not evidence of existence. The scores are research judgments; a 3/5 is not a probability.

| Route | Plausibility | Parameter freedom | Open status | Computation | Certificate | Near a known four | Verdict |
|---|---:|---:|---:|---:|---:|---:|---|
| A. Shi maximal focus + finite outer loop | 0 | 1 | 0 | 3 | 3 | 5 | CLOSED as stated; broader Shi4+1 needs reformulation |
| B. Original Q4 endpoint + interior cycles | 3 | 3 | 4 | 2 | 2 | 1 | UNKNOWN; real two-saddle infinity transport missing |
| C. Generic Q4: M1=0 enlarges basis | 0 | 0 | 0 | 3 | 3 | 0 | CLOSED as a distinct escape |
| C2. Exceptional higher orders / multiplicity splitting | 2 | 3 | 4 | 2 | 2 | 0 | UNKNOWN; identify actual center stratum and common arc |
| D. Raw finite-distance Galias–Tucker continuation | 3 | 4 | 3 | 1 | 3 | 5 | LIVE; terrible original scale separation |
| E. Generic infinity/polycycle exploration | 2 | 3 | 3 | 2 | 2 | 1 | UNKNOWN until a precise graphic is specified |
| E2. 2025 resonant hemicycle a=-1 | 2 | 3 | 4 | 3 | 2 | 0 | UNKNOWN; individual bounds are not additive |
| F. KKL first-order Hopf completion | 3 | 3 | 4 | 4 | 4 | 4 | **Best next bounded strike**, LIVE but no precursor |
| Hamiltonian two-saddle / alien control | 0 | 2 | 0 | 4 | 4 | 0 | Five locally excluded under GI hypotheses; excellent control |

The recommendation is F because the candidate implication can be checked with four persistent cycles and one local Hopf coefficient, and the incumbent numerical seed has much stronger return-map signals than original Songling. It does not win because this repository already contains code. The competing Q4 endpoint problem has substantive new geometric justification and is the best analytic fallback.

## What survives and what does not

**Shi/Bautin.** An elementary quadratic weak focus has local cyclicity three. The exact QW3 stratum is m=5a, b=3l+5. The second-focus condition forces the other finite pair to be complex. Li excludes a preexisting cycle around the order-three focus; QW3 infinity graphics lie around the other focus. Adding two cycles there while three small ones persist would be (3,2). Leaving this stratum is necessary for a new mechanism. [Exact focus, finite and infinity calculations](frontier_2026_09_04/SHI_TOPOLOGY_AUDIT.md).

**Q4 endpoint.** Theorem N is a distinct-interior-zero theorem. In original coordinates the relevant endpoint has two infinity saddle directions, although one transformed elliptic saddle represents them. Neither 3 interior+2 endpoint nor 4+1 is constructed or excluded by the present audit. The latter also needs a four-zero interior integral that has not been supplied. In the strict lobe the logarithmic endpoint coefficient is nonzero; this obstructs a particular two-shadowed-zero degeneration, not every alien mechanism. [Coordinate map and precise Lane B defects](frontier_2026_09_04/Q4_GRAPHICS_AUDIT.md).

**Higher Melnikov.** For fixed generic Q4, every first nonzero Melnikov function lies in the same four-dimensional generating space. Distinct-zero bounds do not automatically bound multiplicity splitting or nonuniform endpoint limits. Symmetric/center-intersection strata have different essential orders. No degree-three perturbation of a quadratic Hamiltonian is admitted as a quadratic counterexample. [Bautin ideals, formulas and essential orders](frontier_2026_09_04/Q4_GRAPHICS_AUDIT.md#3-higher-order-melnikov-exact-finite-dimensional-statement).

**Hamiltonian graphics.** The nondegenerate quadratic Hamiltonian two-saddle closed annulus has cyclicity <=3 under arbitrary quadratic deformation, and explicit alien-two controls exist. The associated theorem gives mutually exclusive endpoint/interior alternatives; adding their separate maxima is wrong. Finite homoclinic cyclicity two and non-Hamiltonian exceptions must retain their original hypotheses. [Graphic ledger, formulas and restrictions](frontier_2026_09_04/Q4_GRAPHICS_AUDIT.md#4-quadratic-hamiltonian-graphics-the-additive-loophole-is-explicitly-obstructed).

**Infinity.** For a≠−1, the 2025 two-hemicycle family has individual cyclicity two and simultaneous cyclicity only two or three. On a=−1 the individual lower bound is two; the sharp upper bound remains unresolved in that paper. The resonance coalesces exponents; it is not a five-cycle existence argument. Saddle-node, nilpotent and semihyperbolic graph results often give only existential finiteness. The H^3_14 2026 preprint is still a claim about one graphic, not H(2). Exact normal forms and retained uncertainties are in the [full appendix](frontier_2026_09_04/Q4_GRAPHICS_AUDIT.md#5-accepted-2025-infinity-result-and-its-exact-exceptional-lane).

**First-order Hopf completion.** The possible escape from Bautin is one local cycle plus three *preexisting*, nonzero cycles in the same nest. For the selected KKL Hopf slice, K>0 gives a repelling weak focus. If the required S/U/S nest and remote cycle coexist hyperbolically, a sufficiently small beta<0 adds an unstable inner cycle. The remote cycle must survive at the same coefficients. No audited theorem excludes this order-one precursor, but its existence is exactly the missing mathematics. [Complete bounded task](ASTRA_FIFTH_STRIKE_HANDOFF.md).

## Red-team rules for every route

- Check the full vector field remains degree at most two under every coordinate/perturbation transport.
- Count isolated orbits, not center ovals, repeated traversals or disconnected signs across a return-map pole.
- Require one coefficient vector and simultaneous survival; never add cycles at different parameter values.
- Use the repaired distribution theorem and the stronger four-real-equilibrium restriction before optimizing.
- Ordinary folds can create two cycles. Six hyperbolic cycles would exceed the lower bound too; do not discard them merely because the motivating question says fifth. The selected strike uses a single Hopf to avoid certifying a semistable fold orbit.
- Every bounded failure is a statement about the specified box/branches only. Missing samples, nonreturns and failed validation remain unresolved, not exclusions.
