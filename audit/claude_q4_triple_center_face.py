#!/usr/bin/env python3
"""Fable lane, step 5: combine Zhao's near-center triple zero with the boundary
face. At fixed a, X(1) is affine in (A,B,eta). Impose Y0 = e1 (small), eta = e2
(small), X(1) = 0; solve for (A,B). Then count interior zeros of I with dense
sampling near both endpoints. Closed count four (three interior + boundary
zero) or four interior would be new."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "q4"))
import numpy as np
from q4_reconstruction import reconstruct, original_values
TE = 1-1e-7
ts = np.concatenate([np.logspace(-6, -1, 300), np.linspace(0.1, 0.9, 200), 1-np.logspace(-1, -7, 300)])
def X1(a, A, B, eta): return reconstruct(a, A, B, eta, t_end=TE).sol(TE)[3]
def count(a, co):
    sol = reconstruct(a, *co, t_end=TE); I = original_values(a, sol, ts); s = np.sign(I)
    idx = np.nonzero(s[:-1]*s[1:] < 0)[0]
    return len(idx), [float(ts[i]) for i in idx]
for a in (0.3, 0.6, 0.9, 0.97):
    # affine model X1 = x0 + xA*A + xB*B + xE*eta
    x0 = X1(a, 0, 0, 0); xA = X1(a, 1, 0, 0)-x0; xB = X1(a, 0, 1, 0)-x0; xE = X1(a, 0, 0, 1)-x0
    for e1 in (0.0, 3e-5, -3e-5):
        for e2 in (0.0, 1e-2, -1e-2):
            # Y0 = 3(1326A+864B-2431 eta-102)/1361360 = e1 ; X1 = 0
            M = np.array([[1326.0, 864.0], [xA, xB]])
            rhs = np.array([e1*1361360/3+2431*e2+102, -x0-xE*e2])
            A, B = np.linalg.solve(M, rhs)
            n, locs = count(a, (A, B, e2))
            flag = " <==" if n >= 3 else ""
            print(f"a={a} Y0={e1:+.0e} eta={e2:+.0e}: A={A:.5f} B={B:.5f} X(1)={X1(a,A,B,e2):+.1e} interior zeros={n} at {np.round(locs,5)}{flag}")
