# H16P: a rational control at the degenerate-infinity organizer

## Status and target

No counterexample. The target is one real planar polynomial vector field of degree at most two with at least five distinct isolated periodic orbits. The advance in this pass is a simpler, more robust numerical four-bracket control in the repaired a=-2 chart, plus finite-field preflights spanning a change in the equilibria at infinity. Neither fifth-cycle existence nor an interval certificate has been obtained.

Repository inspected: main 498b58fc1378486601f348f67e0ab3a5362e7933; PR #6 head 31d8a41a5c7f93ad9fb2d5955f269557abf2d419. The recent report and certification code are on the PR. Earlier repo narratives contain disputed closure statements; this plan does not rely on those statements.

## What is new in this pass

The inherited repaired-chart moment search at a=-2,b=1/3 saved the first preferred 3+1 pattern rather than maximizing its sign margin. A bounded search over the sampled line arrangements gave m approximately 0.88903654 and c approximately 1.98749729. Nearby exact rational values m=8/9, c=159/80 retain the pattern.

Independent mpmath quadrature at 50 and 80 decimal digits, using direct reciprocal-height integration rather than the inherited logarithmic Gauss quadrature, gives these normalized first-order endpoint values (rounded):

| Side | Energy bracket | Endpoint values |
|---|---|---|
| Upper | (-1/6,-1/12) | +0.000509730, -0.000500080 |
| Upper | (3/5,4/5) | -0.000641347, +0.000700659 |
| Upper | (1,4/3) | +0.001157808, -0.002893247 |
| Lower | (5,11/2) | -0.116184456, +0.273987461 |

These are numerical values, not interval bounds. The minimum displayed absolute endpoint value exceeds 5e-4, versus about 3.8e-9 in the previously replayed inherited control. Different witness intervals are used, so this is a comparison of the available controls, not a condition-number theorem.

## Exact proposed family

For positive rational delta and rational A use

    x' = -5/12 + (2/3)y + (-2+A delta^2)x^2 + y^2/3 - (159/80)delta^3 x
    y' = -delta^3/3 - 2xy + (8/27)delta^3 x^2.

At delta=1/50 this is

    x' = -5/12 + (2/3)y + (-2+A/2500)x^2 + y^2/3 - (159/10000000)x
    y' = -1/375000 - 2xy + x^2/421875.

The preflight used delta=1/50 and 1/100, and A=-1,-1/2,0. All six fields retained the four numerical endpoint sign brackets at relative tolerances 3e-12 and 3e-13. This is not proof of four isolated cycles or persistence throughout an interval of A.

## Why the A direction is specific

For the repaired family, set u=1/x,v=y/x and use the signed desingularization d/ds=u d/dt. On u=0 the exact angular equation is

    v' = gamma - (a+2)v - b v^3.

With v=delta z and the displayed coefficients it is delta^3 times

    8/27 - A z - z^3/3.

The discriminant vanishes at 4A^3+64/81=0, hence A*=-(16/81)^(1/3), between -1 and -1/2. The compact numerical brackets survive at tested values on both sides. The cubic counts equilibria at infinity, not periodic orbits. It supplies an organizer for the experiment, not a cycle-birth theorem.

## Attack and decisive gates

1. Establish a controlled four-cycle region near delta=1/50 and A in [-1,-1/2]. Continue the four roots, compare section partners, and check nondegeneracy. Replay with independent high precision. Validate the whole return domains and endpoint signs before claiming a certified count. The leading moment signs alone are insufficient at finite delta.

2. Determine the actual outer itineraries. Track the saddle branches and regular transitions on both sides of A*. Construct the full return from these passages with consistent time orientation. The old fixed-sign crossing assumption on y=0 cannot be imported: here y'=(delta^3/27)(8x^2-9) on that line. A return-domain failure is not an extra zero. A changed infinity phase portrait is not yet a cycle birth.

3. Test the outer full-return displacement for an additional simple zero or a distinct nondegenerate fold while the four old cycles survive. Search across the organizer and on both sides, with log-height or compactified coordinates. Use the endpoint expansion to select scales; do not merely increase a Cartesian radius ceiling. A simple fifth cycle suffices; a fold creating a pair also qualifies if all four survive.

   The minimal sign target is particularly concrete. Using the base-energy labels only to identify initial section points (they are not conserved energies of the perturbed field), the finite preflight has upper signs +,-,+,- at labels -1/6,3/5,1,4/3, and lower signs -,+ at labels 5,11/2. Seek one further upper section point, exterior to all those upper points, with positive displacement. A continuous regular return domain connecting the five upper witnesses would give four upper zeros by IVT; the lower pair would give a fifth. Isolation and distinctness remain separate obligations. Thus a new fold is not the only target: one additional outer sign witness can suffice.

4. If the proposed one-direction family cannot meet the matching condition, test a bounded neighborhood in the exact controls m,c,b while maintaining the compact sign inequalities. An explicit compatibility obstruction can terminate this particular route. Do not assume separate local cyclicities add or interpret an optimizer residual as a construction.

5. Freeze exact rational coefficients at any credible five-or-more numerical candidate. Verify existence, isolation, and geometric distinctness for those same coefficients by validated integration. An exact global cycle count is unnecessary.

## Computational formulation and its limits

The finite preflight uses r=1/abs(y), U=x/y, R=log(r), side sigma=sign(y), and the positive time change d/ds=r d/dt on each half-plane. Writing alpha=a+2 gives

    R' = sigma(2U-e0 r^2-gamma U^2)
    U' = sigma b+(1-b)r+sigma(b-2)r^2/4+sigma alpha U^2
         +e1 r U-sigma e0 U r^2-sigma gamma U^3.

The coordinate formulas were checked symbolically for both signs. Half-return events and launch/terminal transversality are numerical gates only. This chart cannot be used through y=0; a complete outer itinerary needs overlapping charts. The preflight has no validated error bounds or derivative certificates and uses the same DOP853 implementation at both tolerances.

## Reproduction

- h16p_rational_boundary_control.py and its JSON: independent first-order quadrature with eight rational energy witnesses at two precisions.
- h16p_finite_boundary_preflight.py and its JSON: six finite fields, eight witnesses per field, two tolerances, 192 half-return computations. Each half has a 200-unit chart-time cap.
- h16p_boundary_scaling_check.py and its JSON: exact infinity-chart identities and the earlier inherited-control replay. That file retains the old m,c values as historical evidence; the improved control above supersedes them for the proposed attack.

Most valuable next result: a full outer-return matching calculation showing either an extra zero compatible with the four compact brackets, or a specific obstruction. No such result has yet been obtained.
