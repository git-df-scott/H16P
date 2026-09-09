# Validated integration layer

A validated (interval) Taylor integrator for planar quadratic fields, and a
validated return map for the reversible-reseed seed family. This is the
certification machinery whose absence blocked every rigorous claim elsewhere
in this repository.

```bash
python3 certify/certify_seed.py
```

- `taylor.py`: interval arithmetic on `mpmath.iv` (outward rounded, arbitrary
  precision), Taylor coefficient recurrences for quadratic fields, a
  Picard-Lindelof rough enclosure (`z0 + [0,h] f(B)` contained in `B`, which
  proves the solution exists on `[0,h]` and stays in `B`), and an order-K step
  whose remainder is bounded by re-running the recurrence over `B`.
- `poincare.py`: validated half returns for
  `P = (b-2)/4 + (1-b)y + a x^2 + b y^2 + e1 x + e2 xy`, `Q = e0 - 2xy`, on the
  section `{x = 0}`. The crossing is localised rigorously: the step's Taylor
  polynomial in `t` is used with its validated remainder, `x'` is checked to
  hold one sign across the whole step so the crossing is unique, and the
  crossing time is then bisected with interval evaluations.
- `certify_seed.py`: displacement enclosures at the endpoints of the seed's
  candidate brackets. Writes `data/certify_seed.json`.

## Measured behaviour

| Check | Result |
|---|---|
| Harmonic oscillator, one full period, thin initial data | enclosure width `4e-36` (315 steps, dps 40) |
| Chen-Wang field, one revolution, thin initial data | width `7.6e-18` |
| Interval initial data, one revolution | wrapping amplification about `267x` |
| Seed field, one half return | width `3.8e-21` in about 4.4 s (673 steps) |

## What it establishes, and what it does not

At `s0 +- 1e-3` around each of the three upper candidate roots, the enclosure
of `D` is strictly signed, with opposite signs at the two endpoints:

| bracket | `D(s0 - 1e-3)` | `D(s0 + 1e-3)` | enclosure width |
|---|---|---|---|
| upper inner | positive | negative | `4.9e-29` |
| upper middle | negative | positive | `1.7e-25` |
| upper outer | positive | negative | `7.8e-21` |

**Those six signs are rigorous.** Concluding that a root lies between each
pair additionally requires `D` to be defined and continuous across the whole
bracket, which needs the return map validated for *every* initial condition in
`[s0-1e-3, s0+1e-3]`, not just the two endpoints. **That step is not done**, so
this is not yet an existence proof, and none of the four candidate cycles is
certified. Completing it means propagating the segment as interval initial
data (subdivided, since the measured wrapping is about `267x` per revolution).

The lower bracket is not attempted at all: at `|y| ~ 9355` the field has
`x' ~ b y^2 ~ 3e7`, which drives the validated step below any usable floor.
That needs a validated logarithmic reformulation, which is not implemented.
Its omission is a limitation of this integrator, not a statement about the
orbit.
