# Binary128 precision handoff

One call, Python 3 plus GCC/libquadmath (no Python numerical package needed by the interface):

```python
from full_return128 import full_returns
result = full_returns(
    coefficients=['0','0','1','1','1','0',
                  '0','-37','0','-10','11/5','969/1000'],
    ray_angle='0',
    log_radius_grid=['-20','1.9','2.3','37','44'],
    tolerance='2e-27',
)
```

CLI: `python full_return128.py request.json > results.json`, with the same keyword
names as JSON keys. `example_request.json` is ready to call.

Coefficient order: P and then Q, each `[1,x,y,x^2,xy,y^2]`. Angles are radians,
counterclockwise from positive x. Log radii are natural logs of **physical**
Euclidean distance from `center`, default `(0,0)`. The center must be an exact
rational equilibrium; translation is done with Python `Fraction` before rounding
to binary128. Each trajectory is integrated forward in time through one full turn;
the clockwise/counterclockwise sense is selected from the initial vector field.

All precision-sensitive inputs must be integer, Fraction or string. Floats are
rejected. Arbitrary rational coefficients are rounded directly to binary128,
not via double. Output displacements are 36-digit decimal **strings**; converting
them to float is the caller's decision. They are numerical approximations, not
36 certified digits. Results preserve the input rational coefficient vector.

Returned quantities:

- `log_displacement = log(R_return/r_initial)`.
- `radial_displacement = r_initial * expm1(log_displacement)`.
- `status`, integration direction and evaluation/step/rejection counters.

Statuses: `OK_NUMERICAL`, `ANGULAR_CHART_UNRESOLVED`, `EVALUATION_LIMIT`,
`NONFINITE_OR_RANGE`. Every failure has a null displacement. A failure does NOT
prove absence of a return or a cycle. This is a full angular-return routine,
not two-sided half-map matching; it requires angular monotonicity around the
chosen equilibrium. It does not provide a Cartesian fallback for turning orbits.

`y_scale='auto'` conditions the linear off-diagonal terms using
`sqrt(-Q_x/P_y)` when positive, otherwise 1. You can override it with a positive
rational/decimal string. This only changes internal coordinates (`y=scale*Y`),
not the physical ray or radii. Default evaluation budget is 500,000 per radius.
The internal log-radius guard is 2000. Tolerances below 1e-32 are rejected.

Modified-midpoint extrapolation uses binary128 throughout and explicitly rejects
nonfinite intermediate values, including error estimates. Repeat relevant signs
at a tighter tolerance. Tolerance is a local error target, not a global enclosure.
No outward rounding or interval certification is performed.

Validation: `python validate_full_return128.py`. Analytic controls cover a
1e-17 linear-focus displacement, log radii above 36, coefficients of magnitude
1e14, an anisotropic center on a nonzero ray, the nonlinear rational reversible
center from Proposition A's correction, and a chart failure. These establish
instrument controls, not coverage of every nonlinear field in the three regimes.

## Reproduction and provenance

- `audit_proposition_a.py`: independently derived exact elimination and reflection identities.
- `run_center_sign_map.py`: 24 fixed K values, folded pair continuation and full returns.
- `extend_domain_edge.py`: extends the returning endpoint on exactly those same fields.
- `inherited_half_quad.cpp` and `inherited_center_events.jsonl`: copied without changing
  their mathematical content from commit `17e54e5d136ac5b044f4a53b76d91d4808cb04e2`,
  paths `fold_surface_2026_09_05/half_quad.cpp` and
  `fastra_d1_2026_09_05/events_positive_center.jsonl`.
  The inherited routine is used only for fold correction. All reported sign-map
  displacements come from the new general full-return interface.
- Source branch reviewed: `claude/conjecture-progress-report-ixsmgv`,
  commit `e3594aa` (2026-09-05).

The drivers append JSONL records so partially completed numerical work survives.
To repeat a complete run, move aside previous output files first.
