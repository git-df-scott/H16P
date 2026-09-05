# Reversible re-seed evidence

Main result and scope: [../REVERSIBLE_RESEED_2026_09_05.md](../REVERSIBLE_RESEED_2026_09_05.md).

Run from the repository root:

```bash
python reversible_reseed/check_geometry.py
python reversible_reseed/check_arrangement.py
python reversible_reseed/moment_search.py
python reversible_reseed/boundary_search.py
python reversible_reseed/verify_control.py
```

- `check_geometry.py`: exact conjugacy, equilibrium determinants, first
  integrals, infinity eigenvalues, the full unfolding determinant and its
  repair at `a=-2`, and the flux identity.
- `check_arrangement.py`: synthetic five-root controls in both annuli and
  a parabola negative control. These are not quadratic vector fields.
- `moment_search.py`: 54 finite shape samples, each with two 41-point
  moment curves; same-parameter line-arrangement search.
- `boundary_search.py`: 10 finite compact-moment probes at `a=-2` and
  `a=0`, with the missing `gamma x^2` direction at `a=-2`.
- `verify_control.py`: independent 65-digit quadrature at six fixed
  energies, followed by 24 original-field return differences at two
  rational perturbation sizes and two tolerances. Asserts all four
  numerical sign brackets persist.

The geometry identities are exact symbolic calculations. All moment-search
absence statements and ODE results are **NUM**, not interval certificates
or uniform exclusions. No five-cycle field was found and the full
reversible route remains open.

Dependencies are listed in `requirements.txt`; actual versions are recorded
in the output. `data/MANIFEST.json` records hashes of this strike's files.
Timing fields will change on replay. The files require no network access,
credentials, or inherited KKL/Chen–Wang calculations.
