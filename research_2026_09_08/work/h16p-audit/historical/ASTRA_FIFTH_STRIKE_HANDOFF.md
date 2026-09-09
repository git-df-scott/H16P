# Astra Strike #5 — first-order Hopf completion of a (3,1) precursor

## Decision and honest starting point

Find a smooth real quadratic **first-order** weak focus with three preexisting, nonzero hyperbolic surrounding cycles and a fourth hyperbolic cycle at the other focus. A sufficiently small Hopf perturbation then supplies cycle five. This precursor has not been found. The target lies outside the published KKL four-cycle wedge; survival of four cycles along a path is an obligation, not an assumption.

This strike was chosen after the Shi maximal-focus/finite-saddle construction failed exact topology checks. Q4/Theorem N remains a successful closed five-interior-zero mechanism. This is not another Q4 strike.

## Prompt-ready specification

Carry out the following bounded construction or exclusion task. Work in exact rational quadratic fields

    x' = y + x² + xy,
    y' = −10x² + (11/5)xy + c y² + alpha x + beta y.

The precursor slice has beta=0 and

    1/2 <= c <= 3/2,
    −200 <= alpha <= −10,
    K = −alpha(11c/5−1)−42 >= 1/64.

The particular lower margin 1/64 is an experimental conditioning choice. Failure with that margin does not exclude K closer to zero. Fixing b=11/5 likewise restricts the experiment, not all quadratic fields.

At the origin det(J)=−alpha>0, tr(J)=beta, and on beta=0 the normalized cubic radial coefficient is

    l1 = K / [8(−alpha)^(3/2)] > 0.

For the final unfolding use one exact rational beta in [−2^−10,0), with |beta| made smaller (beta moved closer to zero) if required by the computed persistence/Hopf bounds. In normalized time the radial linear coefficient is beta/[2 sqrt(−alpha)], with nonzero derivative at beta0. Thus a small negative beta creates an **unstable** inner cycle. The predicted normalized radius squared is 4 alpha beta/K to leading order; this is only an initial guess, never a certificate.

### Exact geometric gates

At beta=0 the nonzero equilibria solve

    T(x) = (c−61/5)x³ + (alpha−111/5)x²
           + (2alpha−10)x + alpha = 0,
    y = −x²/(1+x).

Require disc(T)<0, its unique real root x*<−1, and at E*=(x*,y*)

    trace(J*) < 0,
    det(J*) > 0,
    trace(J*)²−4det(J*) < 0.

Use exact polynomial/root isolation and interval signs. The additional exact pruning condition is c>11/20. In the one-real-remote-root regime, its negative trace requires m=−alpha greater than 21(1000c²+1021c+481)/[50(2c+1)²(8−5c)]; keep the Jacobian gate as the authoritative check. These conditions retain one repelling weak focus at origin, one attracting strong focus remotely, and a nonreal additional pair. At final beta repeat the equilibrium check using

    T_beta(x)=(c−61/5)x³+(alpha−111/5−beta)x²
              +(2alpha−10−beta)x+alpha.

The line x=−1 has x'=1 everywhere. A periodic orbit cannot cross it; the two nests stay in different half-planes. The direct finite connection between the antipodal vertical infinity saddles is not available. Do not resurrect a finite saddle or that vertical connection as the missing cycle.

In the vertical compactification x=u/z, y=1/z, the equations are

    u' = 10u³−(6/5)u²+(1−c)u + z(1−alpha u²−beta u),
    z' = z(10u²−(11/5)u−c) − z²(alpha u+beta).

The vertical eigenvalues are (1−c,−c). Extra real infinity directions appear at c=241/250; the vertical direction degenerates at c=1. Treat the three open c-regions and their two boundary strata separately. A numerical branch that reaches either boundary requires its compactified itinerary to be checked before continuation.

### Bounded numerical task

1. Reproduce the incumbent rational seed c=7/10, alpha=−363889/5000, beta=3/2000. The positive downward y=0 return coordinates are approximately 0.6832102174, 2.1836998249, 15.9627839816; the remote downward crossing is approximately −3711.5608064. Their multipliers are S/U/S and U. Reproduce actual returns and stability, not trajectory spirals.
2. Continue beta down to zero at the incumbent shape. Expect the inner stable cycle to disappear. This audit finds only two nonzero origin roots, approximately 3.0688454246 and 15.0640714513. Here K<0. This expected loss is a negative control: retaining a purported third nonzero root here is an error.
3. Start the K>0 component from c=7/10, alpha=−80, beta=0, where this audit locates an origin return root near 64.55543434 with multiplier 0.80969114. This is only a **one-cycle control**, not the desired precursor. The remote unstable return is numerically near −5391.14116 (multiplier12.1680); certify it separately. Sparse negative-section scans include return-itinerary changes and must not be counted as extra cycles. Continue this return-root sheet in the two shape variables within the fixed rectangle and continue its fold curves using D=0 and partial_r D=0; use second derivatives to orient folds. Track existing branches and compactified escape boundaries. The sought event is coexistence of **three nonzero** simple roots in this K>0 slice; a fold can add the missing pair. Do not count the origin root r=0.
4. Use the downward section y=0 with positive origin coordinates r in [2^−12,2^10], and remote r in [−2^20,−1] subject to r<alpha/10. Work in log(r) or bounded section coordinates to avoid a uniform raw grid. Each sign change must belong to a single continuous full-return itinerary. At most 256 continuation steps per seeded branch and 4096 total return/derivative evaluations in this strike; at most 64 adaptive parameter cells around actual folds or changes of return-root count. No unseeded global coefficient optimization. These bounds limit the experiment and do not produce an exhaustive parameter cover.
5. At any candidate precursor require three ordered, disjoint positive return intervals J1<J2<J3 with stability S/U/S, and a disjoint remote interval J4 with stability U. Require all four returns and their derivatives to persist on one common parameter neighborhood at beta=0. Check every cycle's enclosing focus and isolation. If no precursor survives, report exact explored branches, unresolved folds/nonreturns, and the box-specific negative result. Do not expand b or the box in this strike.
6. Only after a precursor is validated, choose rational c,alpha within its validated open neighborhood and a sufficiently small explicit negative rational beta. Use an analytic Hopf estimate with a rigorous remainder, or validated shooting on a rescaled small section, to propose an innermost interval J0 disjoint from J1. Preserve all four original interval gates at this SAME beta. The final field must have five independently verified limit cycles; an asymptotic “sufficiently small beta” statement alone is not the final artifact.

The derivative for fold continuation must be the actual return derivative. On y=0,

    P'(r) = Q(r,0)/Q(P(r),0) * exp(integral div(F) dt).

The bare divergence exponential equals the multiplier only at a periodic fixed point. Use a validated variational equation plus the section projection when proving the bound.

### Exact certification gate

For one rational tuple (c,alpha,beta), construct five disjoint compact transverse section intervals. For each Ji:

- Prove the complete first return exists for every initial point in Ji, has the selected crossing orientation, and stays on one itinerary with no singularity or escape.
- Enclose P at the endpoints and show strict opposite signs of D=P−identity.
- Enclose P'(Ji) and exclude 1, or verify interval Newton N(Ji) lies strictly inside Ji.
- Prove distinctness by common-section order and the two half-plane/nesting assignments. Repeated traversals and two sections of one orbit do not count twice.

Save exact coefficients as rational strings, outward-rounded interval endpoints, return-time enclosures, derivative bounds, precision/order, all accepted flow boxes, software/source hashes and a clean replay command. Obtain five successful gates using MPFR-backed interval integration. A failed enclosure is unresolved, not a rejected mathematical candidate. See [H16_CERTIFICATION_PLAN.md](H16_CERTIFICATION_PLAN.md).

### Hostile kill conditions

- K=0 or a multiple equilibrium: wrong generic Hopf target; do not call it first-order.
- Three roots supplied only by a local truncated polynomial: no precursor established.
- Three origin cycles at one parameter and the remote cycle at another: invalid sum.
- Extra finite saddle while (3,1) remains hyperbolic: conflicts with the four-real-equilibrium theorem.
- The remote cycle splits into two while the three survive: forbidden (3,2).
- A return sign jumps across an escape/nonreturn boundary: no intermediate-value argument.
- Loss of an existing cycle before beta<0 reaches the proposed fifth: stop or reduce beta and revalidate.
- Exhausted bounded continuation without a candidate: write a reproducible negative/unknown report and stop. No claim of H(2)<=4 follows.

### Required outputs

`STRIKE5_PRECURSOR.md`, exact coefficient JSON, a continuation-event ledger, four precursor and five final return certificates if successful, replay logs, and a clear mathematical or box-specific obstruction if not. Preserve Q4 and Theorem N unchanged.

## Why this survived the hostile audit

Bautin limits *local* cyclicity: this proposal creates one local cycle and requires three preexisting nonzero ones. Li's order-three exclusion and order-two uniqueness do not apply at K>0. The (4,1) configuration is permitted by the published distribution restrictions. Targeted checks of Zhang–Cai 1991, Zhang–Zhao 2001 and Zegeling 2024 found no theorem excluding the order-one precursor. This is absence of a located obstruction, not a proof it exists. The main practical risk is that all reachable three-cycle regions lie in K<0 and lose a cycle before K changes sign.

The audit's 12-point beta0 probe found no precursor. Its only positive-root detections were the two incumbent-shape roots above and the one K>0 root at alpha=−80. Unsampled intervals, large cycles and nonreturns remain unresolved. The experiment is a control and a starting sheet, not a parameter-space exclusion.

Independent hostile review: [STRIKE_REDTEAM.md](frontier_2026_09_04/STRIKE_REDTEAM.MD).
