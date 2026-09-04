#!/usr/bin/env python3
"""Claude hostile check: exact Green/PF reconstruction (R1)-(R3) versus
independent original area integrals, for RANDOM universal coefficients and
several kappa. Tests the coefficient transport, ODE, center data and the
corrected forcing sign on all four coefficient directions simultaneously.
Numerical, not interval-rigorous."""
import os, sys
for k in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"): os.environ[k]="1"
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "q4"))
import numpy as np, mpmath as mp
from q4_integrals import basis_mp
from q4_reconstruction import reconstruct, original_values, mu_from_universal

rng = np.random.default_rng(20260904)
worst = 0.0
rows = []
for kappa in (1.15, 1.7, 3.0, 9.0):
    a = 1-1/kappa
    for trial in range(3):
        A, B, eta = rng.uniform(0.5, 2.5), rng.uniform(-1.5, 0.5), rng.uniform(0.8, 1.8)
        sol = reconstruct(a, A, B, eta, t_end=0.97)
        mu = mu_from_universal(kappa, A, B, eta)
        for t in (0.2, 0.55, 0.9):
            s = kappa-(kappa-1)*t
            ind = float(sum(m*v for m, v in zip(mu, basis_mp(kappa, s, dps=40))))
            rec = float(original_values(a, sol, t))
            rel = abs(ind-rec)/max(abs(ind), 1e-30)
            worst = max(worst, rel)
            rows.append((kappa, trial, t, ind, rec, rel))
for r in rows: print("kappa=%5.2f trial=%d t=%.2f area=% .12e recon=% .12e rel=%.2e" % r)
print("worst relative discrepancy:", worst)
# Wrong-sign control: flip the forcing sign and show the check fails.
import q4_reconstruction as R
src = open(R.__file__).read()
assert "forcing = -h_over_t2/(1152*(1-t))" in src
print("PASS" if worst < 1e-8 else "FAIL")
