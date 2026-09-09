# Outer-fold strike, 2026-09-09

Execution of Steps 1 and a bounded Step 2/3 reconnaissance from
[../H16P_NEXT_STRIKE.md](../H16P_NEXT_STRIKE.md).

Run from the repository root:

```bash
python3 outer_fold_2026_09_09/step1_seed.py
python3 outer_fold_2026_09_09/step2_outer_scan.py
```

- `step1_seed.py`: re-establishes the four-cycle seed. Two independent
  solvers on the same field (Cartesian state `(x, y)`, logarithmic state
  `(x, log|y|)` — a different ODE, not a rescaling of the first), section
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

Evidence class is **NUM** throughout: floating-point ODE integration and
bisection, no interval arithmetic, no validated integration. Nothing here
establishes existence or isolation of any cycle in the rigorous sense, and
no fifth cycle is asserted.

Results and their limits: [../OUTER_FOLD_RESULT_2026_09_09.md](../OUTER_FOLD_RESULT_2026_09_09.md).
