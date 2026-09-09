# General-quadratic cycle counting probe, 2026-09-09

A side probe, separate from the reversible-reseed continuation in
[../outer_fold_2026_09_09/](../outer_fold_2026_09_09/). Same evidence class:
floating point, no enclosures.

## Why this target

Two published results narrow what a five-cycle quadratic field can look like:

- **Zhang Pingguang**: the limit cycles of a quadratic system with two foci
  are distributed `(0,1)` or `(1,i)` — one nest holds *exactly one* cycle.
  Five cycles therefore require **four around a single focus**.
- **Li Chengzhi (1986)**: no limit cycle surrounds a third-order weak focus of
  a real quadratic system. Since three small-amplitude cycles need exactly
  that focus (Bautin caps small-amplitude cyclicity at three), a fourth cycle
  in the same nest must be a large one, appearing at finite distance from that
  stratum.

Both were taken from the repository's existing literature audit and the
[focus-route correction](../research_2026_09_08/outputs/H16P-focus-route-correction.md);
neither proof was independently audited here.

## Contents

```bash
python3 five_cycle_2026_09_09/validate.py   # counter vs. a published field
python3 five_cycle_2026_09_09/search.py     # bounded parameter probe
```

- `counter.py`: nest-aware cycle counting for
  `x' = lam x - y + l x^2 + m xy + n y^2`, `y' = x + lam y + a x^2 + b xy + c y^2`.
  Equilibria by resultant, foci by eigenvalues, and a ray section from each
  focus whose crossing direction is taken from the rotation sense of the
  linearisation rather than assumed. Counts sign changes of
  `D(r) = log(R(r)/r)`. An unresolved return ends the domain and is recorded
  as a missing domain, never as a zero.
- `validate.py`: run against the Chen-Wang visualisation parameters
  (`l,m,a,b = -3, 0.99, 2/9, -3`, `lam = -2e-5`). Finds **three** cycles about
  the origin, outermost at `r = 0.33380`, against the `0.332839` recorded in
  [../STAGED_SHI_2026_09_05.md](../STAGED_SHI_2026_09_05.md). The two inner
  radii differ from that document's and the difference is not explained here.
  The remote cycle about the second focus was not detected at these
  parameters.

## Search outcome, and what it is not

576 parameter samples, structured perturbations of the Chen-Wang and a
conditioned Shi point, 1658 s. **Zero** samples with four cycles in a nest —
and also **zero** samples with even three cycles in a nest.

That second number is the informative one, and it is a verdict on the
sampler, not on the family. The three-cycle configurations sit in a
fantastically thin region (Chen-Wang needs `0 < d2 << d1 << 1`), and
perturbations of 2-90% leave it immediately. **This probe therefore covers
essentially none of the region where a fourth cycle could be found, and
supports no conclusion about quadratic fields at all.** A useful version
would have to move along the weak-focus stratum rather than jump off it.
