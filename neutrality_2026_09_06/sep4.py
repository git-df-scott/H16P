#!/usr/bin/env python3
"""Fate of the two transverse infinite separatrices, in the Poincare disc.
A = infinite saddle with equator-stable/transverse-unstable (orbit leaves into
the plane); B = its antipode (orbit enters the plane).  A heteroclinic
connection A->B closes a graphic together with an equator arc B->A."""
import numpy as np
from scipy.integrate import solve_ivp
from sphere3 import make, inf_points, tangent_eigs

def separatrix(l, a, which, T=800.0, eps=1e-9):
    pts = inf_points(l, a)
    A = [p for u, p in pts if p[0] > 0][0]
    B = [p for u, p in pts if p[0] < 0][0]
    s0, tgt, sgn = (A, B, +1) if which == "A_out" else (B, A, -1)
    w, v = tangent_eigs(l, a, s0)
    k = int(np.argmax(np.abs(v[2, :]))); d = v[:, k]/np.linalg.norm(v[:, k])
    if d[2] < 0: d = -d
    s = s0 + eps*d; s /= np.linalg.norm(s)
    F = make(l, a, sgn=sgn)
    def hit(t, y): return np.linalg.norm(y - tgt) - 1e-5
    hit.terminal = True; hit.direction = -1
    sol = solve_ivp(F, (0, T), s, rtol=1e-11, atol=1e-13, events=hit, dense_output=True)
    return sol, A, B, tgt

def summarise(l, a):
    row = []
    for which in ("A_out", "B_in"):
        sol, A, B, tgt = separatrix(l, a, which)
        d = np.linalg.norm(sol.y - tgt[:, None], axis=0)
        e = sol.y[:, -1]; nd = abs(np.linalg.norm(sol.y, axis=0) - 1).max()
        xy = e[:2]/e[2] if e[2] > 1e-9 else None
        row.append(f"{which}: T={sol.t[-1]:7.2f} drift={nd:.1e} s3end={e[2]:+.6f} "
                   f"dmin_to_antipode={d.min():.3e} "
                   f"xy={'inf' if xy is None else np.round(xy,4)}")
    return row

for (l, a) in [(-10.0, 1.0), (-12.0, 0.8), (-8.0, 1.2), (-8.0, 0.8), (-12.0, 1.2), (-9.0, 1.0)]:
    print(f"=== l={l} a={a} ===")
    for r in summarise(l, a): print("   ", r)
