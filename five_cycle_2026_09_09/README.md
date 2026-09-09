# General-quadratic cycle counting probe, 2026-09-09

A side probe, separate from the reversible-reseed continuation in
[../outer_fold_2026_09_09/](../outer_fold_2026_09_09/). Same evidence class:
floating point, no enclosures.

## Why this target

Two published results narrow what a five-cycle quadratic field can look like:

- **Huang and Reyn (1995)**, Bull. Austral. Math. Soc. 52, 461-474,
  [DOI 10.1017/S0004972700014945](https://doi.org/10.1017/S0004972700014945):
  in a quadratic system with two nests of limit cycles, one nest contains
  exactly one cycle — **but the abstract restricts to systems where the sum of
  the multiplicities of the finite critical points equals three.** An earlier
  version of this file attributed an unrestricted form of this to Zhang
  Pingguang; that was wrong and is withdrawn. The hypothesis must be checked
  before the result is applied to any particular family, and it is *not*
  verified for the reversible seed family used elsewhere in this repository.
  Only the abstract was read.
- **Li Chengzhi (1986)**: no limit cycle surrounds an **exact** third-order
  weak focus of a real quadratic system. This says nothing about a positive
  parameter-distance neighbourhood of that stratum, and must not be read as
  excluding one.

Both were taken from the repository's existing literature audit and the
[focus-route correction](../research_2026_09_08/outputs/H16P-focus-route-correction.md).
Neither proof was audited here, and in the Huang-Reyn case not even the full
paper — only the abstract, which is where the multiplicity hypothesis appears.

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
