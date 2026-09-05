# Resonant compact/endpoint compatibility

Main proof and scope: [../RESONANT_JOINT_2026_09_05.md](../RESONANT_JOINT_2026_09_05.md).

```bash
python resonant/check_resonant.py
python resonant/shoot_control.py
```

`check_resonant.py` checks the first integral, weighted divergence,
resonant invariant, first normal Dulac formula, and the mixed compact
generator. It cross-checks area moments and the mixed derivative by
independent quadrature.

`shoot_control.py` evaluates three fixed sign brackets on each of three
rational members of an explicit parameter arc, using the original
quadratic field at two tolerances. This is NUM evidence only. It does not
prove isolation, give interval enclosures, or count all cycles.

Dependencies used: Python 3.12, SymPy 1.14, mpmath 1.3, SciPy 1.17, NumPy.
Machine-readable logs in `data/` include actual versions.

The proof excludes the stated two-compact mechanism. It does not compute
the full broken-connection resonant second-order expansion or solve H(2).
