# Outer-fold strike, 2026-09-09

Execution of Steps 1 and a bounded Step 2/3 reconnaissance from
[../H16P_NEXT_STRIKE.md](../H16P_NEXT_STRIKE.md).

Run from the repository root:

```bash
python3 outer_fold_2026_09_09/step1_seed.py        # first pass, superseded
python3 outer_fold_2026_09_09/step2_outer_scan.py  # first pass, superseded
python3 outer_fold_2026_09_09/step3_compensated.py # repaired machinery
```

- `machinery.py`: the repaired evaluator — itinerary gates, event-corrected
  variational derivatives in the section coordinate and in the five
  coefficients, fixed residual scales, and a log-coordinate cross-check.
  The Cartesian and logarithmic formulations share DOP853 and the same event
  machinery, so their agreement is a coordinate cross-check, not independent
  verification.
- `step3_compensated.py`: repaired seed audit with derivative uncertainty,
  sensitivity-predicted root tracking inside identity-preserving guards, and
  the compensated continuation (`a` driven, `(b, e0, e1, e2)` compensating)
  searching for a double zero `D = 0`, `dD/ds = 0` beyond the outermost upper
  cycle. Writes `data/step3_compensated.json`.
- `step1_seed.py` (first pass, superseded by `machinery.py`): two coordinate
  formulations of the same field (Cartesian `(x, y)`, logarithmic
  `(x, log|y|)`), sharing the stepper and, in this script, a Cartesian launch;
  section
  `{x = 0}` with coordinate `s = log|y|`, displacement
  `D(s) = s_forward - s_backward`. Refines the four saved sign brackets to
  isolated roots, estimates the section derivative at each, and gates every
  half return on a complete itinerary: exactly one transverse crossing, no
  interior crossing, no approach within `1e-3` of an equilibrium.
  Writes `data/step1_seed.json`. Exits nonzero if the control fails.
- `step2_outer_scan.py`: Part A scans `D(s)` outward beyond the outermost
  tracked upper root at the seed; Part B runs two continuation branches
  (capped at 50 accepted steps) that re-correct all four roots at every
  step and re-scan the outer region for a fold signature.
  Writes `data/step2_outer_scan.json`.
- `data/outer_extension.json`: the outward scan extended to `s = 11` at the
  seed and at the branch-A endpoint.

Evidence class is **NUM** throughout: floating-point ODE integration, Newton
correction and bisection. No interval arithmetic, no validated integration, no
enclosure of any kind. The itinerary gates are numerical safeguards that can
reject a bad return; passing them establishes nothing about existence or
isolation of a periodic orbit. No fifth cycle is asserted, and sampled
profiles are never read as absence statements.

Results and their limits: [../OUTER_FOLD_RESULT_2026_09_09.md](../OUTER_FOLD_RESULT_2026_09_09.md).
