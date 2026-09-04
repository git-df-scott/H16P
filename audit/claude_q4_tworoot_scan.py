#!/usr/bin/env python3
"""Fable lane, step 4: closed-interval first-order zero count in the two-root
region. Primitive H with K3-coefficient one and roots r<s only: a line in
(A,B,eta); parametrize by eta. For each (r,s,eta) with exactly two interior
primitive zeros, and each lift a, count interior zeros of I (log-dense near
the loop) and record X(1). Target: interior count 3 with X(1) near zero
(closed count four), or interior count 4 anywhere."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "q4"))
import numpy as np, mpmath as mp
from q4_threshold_path import primitive_basis_closed, primitive_value_closed
from q4_reconstruction import reconstruct, original_values
mp.mp.dps = 30
ts = np.concatenate([np.linspace(1e-3, 0.9, 300), 1-np.logspace(-1, -7, 300)])
def two_root_coeffs(r, s, eta):
    # (A-1)K0 + B K1 - eta K2 + K3 = 0 at r and s -> solve for A-1, B
    rows = [primitive_basis_closed(mp.mpf(t)) for t in (r, s)]
    M = mp.matrix([[rw[0], rw[1]] for rw in rows]); rhs = mp.matrix([eta*rw[2]-rw[3] for rw in rows])
    sol = mp.lu_solve(M, rhs)
    return float(sol[0]+1), float(sol[1]), float(eta)
def primitive_count(co):
    A, B, eta = co
    grid = np.concatenate([np.linspace(0.005, 0.9, 200), 1-np.logspace(-1, -5, 100)])
    H = np.array([float(primitive_value_closed(mp.mpf(t), (A, B, eta))) for t in grid])
    s = np.sign(H); return int(np.sum(s[:-1]*s[1:] < 0))
best = []
for r in (0.15, 0.3, 0.5, 0.7):
    for s in (r+0.15, r+0.3, 0.9, 0.98):
        if s >= 1: continue
        for eta in (0.3, 0.7, 1.0, 1.2, 1.5, 1.8, 2.5):
            co = two_root_coeffs(r, s, eta)
            if primitive_count(co) != 2: continue
            for a in (0.2, 0.5, 0.8, 0.9, 0.95, 0.99):
                sol = reconstruct(a, *co, t_end=1-1e-7)
                I = original_values(a, sol, ts); sg = np.sign(I)
                n = int(np.sum(sg[:-1]*sg[1:] < 0)); X1 = float(sol.sol(1-1e-7)[3])
                best.append((n, r, s, eta, a, X1))
                if n >= 3: print(f"  ** interior count {n}: r={r} s={s} eta={eta} a={a} X(1)={X1:+.3e} co={np.round(co,5)}")
best.sort(reverse=True)
print("max interior count found:", best[0][0], " sample:", best[0])
from collections import Counter
print("histogram of interior counts:", Counter(b[0] for b in best))
