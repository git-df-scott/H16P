#!/usr/bin/env python3
"""Step 5b: near-center hierarchy on the boundary face. Impose X(1)=0, take
Y0<0 tiny and eta<0 small with |Y0| << |eta|/384 so the two leading
center coefficients alternate against the O(1) quadratic coefficient
(Y2 ~ -q0/2304 < 0). Count interior zeros (dense near both endpoints)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "q4"))
import numpy as np
from q4_reconstruction import reconstruct, original_values
TE = 1-1e-7
ts = np.concatenate([np.logspace(-8, -1, 400), np.linspace(0.1, 0.9, 200), 1-np.logspace(-1, -6.9, 300)])
def X1(a, A, B, eta): return reconstruct(a, A, B, eta, t_end=TE).sol(TE)[3]
def count(a, co):
    sol = reconstruct(a, *co, t_end=TE); I = original_values(a, sol, ts); s = np.sign(I)
    idx = np.nonzero(s[:-1]*s[1:] < 0)[0]
    return len(idx), [float(ts[i]) for i in idx]
for a in (0.3, 0.6, 0.9):
    x0 = X1(a, 0, 0, 0); xA = X1(a, 1, 0, 0)-x0; xB = X1(a, 0, 1, 0)-x0; xE = X1(a, 0, 0, 1)-x0
    for e2 in (-1e-1, -1e-2, -1e-3):
        for e1 in (-1e-5, -1e-7, -1e-9):
            M = np.array([[1326.0, 864.0], [xA, xB]]); rhs = np.array([e1*1361360/3+2431*e2+102, -x0-xE*e2])
            A, B = np.linalg.solve(M, rhs); q0 = A-1-e2/6
            n, locs = count(a, (A, B, e2))
            print(f"a={a} eta={e2:+.0e} Y0={e1:+.0e}: q0={q0:+.3f} interior zeros={n} at {np.round(locs,6)}{'  <==' if n>=3 else ''}")
