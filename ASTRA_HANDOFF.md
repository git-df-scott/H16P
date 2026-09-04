# Astra handoff

## Mission

Determine whether one of three mathematically reduced quadratic families
contains an explicit field with at least five distinct limit cycles. The
campaign is **YELLOW**: it is not authorized to search generic coefficients.

## Binding facts

- `H(2) >= 4`; no finite uniform upper bound is known.
- `H(2)=4` is open as of 2026-09-04.
- Bautin gives at most three small cycles from one quadratic focus/center.
- With two foci, one nest has at most one cycle; a five-cycle target is `4+1`,
  not `3+2`.
- Galias--Tucker prove one Songling field has exactly four by interval methods.
- Candidate verification is tractable; global exhaustive discovery is not.

## Execute in this order

1. `Q4` elliptic-integral zero hunt.
2. Shi third-order weak focus plus outer separatrix continuation.
3. Kuznetsov four-cycle boundary continuation.

Exact families, budgets, and stop rules are in
[`ATTACK_MATRIX.md`](ATTACK_MATRIX.md). Do not silently expand them.

## Required pipeline

```text
exact family
  -> symbolic bifurcation constraints
  -> multiprecision scalar candidate screen
  -> five persistent first-return roots
  -> rationalize coefficients
  -> interval Poincare maps
  -> five distinctness certificates
  -> independent replay
```

## Candidate record

For every promoted candidate, save:

- exact coefficient vector and normalization chart;
- parameter provenance (branch, box, continuation step);
- all equilibria and their interval classifications;
- Poincare sections and first-return conventions;
- root brackets at 128/256/512 bits;
- return times, displacement values, derivative/Floquet estimates;
- reason each pair of cycles is distinct;
- exact failure or proof status.

## Promotion gates

`NUM-CANDIDATE` requires five roots that persist under precision doubling,
integrator change, and section perturbation. `CAP-CANDIDATE` requires five
interval fixed-point proofs. `COUNTEREXAMPLE` requires a clean independent
replay of all five and a human-readable mathematical argument connecting the
machine inequalities to five isolated periodic orbits.

## Automatic rejection rules

Reject and label, do not rescue by longer integration, when:

- the system is piecewise, discontinuous, degree greater than two, or
  three-dimensional;
- a reported return is an iterate or starts at the event surface;
- two apparent cycles are the same orbit on different sections;
- a root disappears under precision doubling;
- a multiplier interval contains one and no multiple-root proof is supplied;
- a `3+2` two-focus picture is reported;
- coefficients are only rounded decimals without outward enclosures;
- evidence consists only of a plot or truncated series.

## Stop discipline

At each attack's budget:

1. freeze raw logs and exact boxes;
2. state `CE FOUND`, `FAMILY EXCLUDED`, or `UNRESOLVED`;
3. record the sharpest obstruction or unresolved boundary;
4. do not roll spare compute into a wider search.

Failure in a narrow family says nothing about global `H(2)`. The next campaign
is justified only by a new analytic reduction, a new bifurcation intersection,
or an independently reproduced five-root candidate.

## Audit checkpoint

This repository contains only the feasibility audit and a four-cycle
floating-point regression. **No five-cycle hunt has been run.**
