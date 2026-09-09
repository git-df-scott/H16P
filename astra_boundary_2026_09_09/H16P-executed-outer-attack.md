# Executed outer attack: no fifth-cycle candidate

## Outcome

The proposed rational a=-2 boundary-organizer family did not produce the missing positive upper return witness in this bounded experiment. All resolved outer samples had negative displacement. This is numerical evidence against the tested construction, not a proof of absence between samples, at other coefficients, or outside the sampled domain.

The stronger rational first-order control (m=8/9,c=159/80) and its six finite-field four-bracket preflights remain numerical positive controls. They are not certified isolated cycles.

## Executed searches

- Initial upper outer probe: 11 A values (-2,-1,-0.65,-0.60,-0.585,-0.5824,-0.58,-0.55,-0.5,0,1), 13 log-heights from 4.04 through 5120; 143 attempted displacement evaluations. The original maximum step was 1. All resolved values were negative.
- The initial range guard could reject an internal Runge-Kutta stage. Such a rejection does NOT show the actual trajectory leaves the chart. The historical failures are retained with this correction.
- Reduced-step followup: A=-1,-0.5824,-0.5,0; log-heights 320,640,1280,2560,5120,10000; maximum step 0.05. All 24 displacement evaluations resolved and were negative.
- Different stepper: RK45 at the same four A values and log-heights 40,640,10000, relative tolerance 3e-12. All 12 values were negative. This changes the numerical method but shares the coordinate derivation and SciPy implementation; it is not a validated or fully independent proof.
- Compact-control feasibility: a linear program on six frozen quadrature witnesses required the 3+1 signs with margin 1e-5. The resulting extrema for m were approximately 0.88062217 and 0.90432578. This is a numerical fixed-witness constraint, not a theorem covering all four-cycle fields.
- Finite-control followup: controls halfway from the rational seed to four LP extrema were tested at A=-1,-0.5824,-0.5. The extrema include near-duplicates, so these are 12 records, not 12 substantially distinct directions. All retained the six sign witnesses for four brackets, and both additional outer evaluations at log-heights 40 and 640 remained negative.

At log-height 10000 the reduced-step upper displacements were approximately -0.15198025, -0.25873537, -0.28921732, and -0.57658675 in the A order above. None is a near-zero fifth-cycle candidate.

## What the boundary algebra actually says

On r=0 in the reciprocal chart, write U>0 for the forward branch and reflect the backward branch by U -> -U. Put B=b+alpha U^2, alpha=a+2. Wherever the relevant denominators are nonzero, their formal dR/dU difference is exactly

    (2U-gamma U^2)/(B-gamma U^3) - (2U+gamma U^2)/(B+gamma U^3)
      = 2 gamma U^2 [(2-alpha)U^2-b] / [B^2-gamma^2 U^6].

This identity was checked symbolically. It shows explicitly where the directional asymmetry enters: gamma, while alpha modifies the even part. The sign can change across zeros of the denominator, so the identity is NOT a global monotonicity proof and cannot be integrated across a saddle without a separate passage analysis. Also R itself is infinite at r=0; this is the limiting coefficient calculation for the reciprocal equations, not a complete physical return on r=0.

The compact-control LP retains positive m and hence positive gamma. Simply reversing that asymmetry is not available while retaining these particular witness inequalities. Larger control changes, moving witnesses, and different b values remain untested by that calculation.

## Limitations and next useful step

The logarithmic solver evaluates exp(R). Below approximately R=-745 it underflows in double precision. The omitted terms are extremely small locally, but their accumulated error is not rigorously bounded. Hence the largest-height calculations are diagnostics, not enormous-radius certificates. Time orientation, global return domains, isolation, distinctness, and interval rounding still require proof.

The tested A variation should not receive a longer blind radius scan. A further run should first establish a new boundary-matching mechanism or change controls within a verified four-cycle region. The present fixed-sign asymptotic identity alone does not exclude a fold pair.

## Concurrent repository update

During this run PR #6 advanced from 31d8a41 to f0e8420e668e73bf413f31bad98ec6084934a423. Its new saved step5 data reports 634 attempted fold solves, none converged. The backward branch completed 19 accepted steps and zero rejections before a wall-time stop. That branch is no longer unrun. I inspected the updated report and machine-read the summary and final backward state; I did not independently replay its 172,980 evaluations.

The report identifies attraction of fold iterations toward existing roots. A stationary-point-first or carefully deflated formulation is worth investigating, but deflation based on approximate moving roots needs its own error and degeneracy checks. The failure statistics do not prove the outer region empty.

## Budget and execution state

The shared main account meter was 58% when queried during this pass. The handoff described a historical 6% baseline and a 30-percentage-point allowance. These cannot establish task-specific consumption, but they do not justify assuming an unused overnight allowance. A fresh cap has been requested before extended execution. All bounded computations described above completed; no search process or monitoring automation was left running.

Raw results are preserved in h16p_outer_probe.json, h16p_outer_refined.json, h16p_outer_rk45.json, h16p_compact_control_constraints.json, and h16p_outer_control_variation.json. The current h16p_outer_return_probe.py uses the repaired maximum step and writes a separate replay file, preserving the historical initial failures.
