#!/usr/bin/env python3
import numpy as np
from scipy.integrate import solve_ivp
from sphere2 import make, inf_points, tangent_eigs

def run(l, a, T=2000.0):
    F = make(l, a); pts = inf_points(l, a)
    A = [p for u, p in pts if p[0] > 0][0]      # equator-stable, transverse-unstable
    B = [p for u, p in pts if p[0] < 0][0]
    out = {}
    for name, s0, sgn in (("A_unstable_fwd", A, +1), ("B_stable_bwd", B, -1)):
        w, v = tangent_eigs(l, a, s0)
        k = int(np.argmax(np.abs(v[2, :])))     # the branch leaving the equator
        d = v[:, k]/np.linalg.norm(v[:, k])
        if d[2] < 0: d = -d
        s = s0 + 1e-8*d; s = s/np.linalg.norm(s)
        sol = solve_ivp(F, (0, sgn*T), s, rtol=1e-11, atol=1e-13, max_step=0.02)
        y = sol.y; nrm = np.linalg.norm(y, axis=0)
        out[name] = (sol, nrm)
    return A, B, out

for (l, a) in [(-10.0, 1.0), (-12.0, 0.8), (-8.0, 1.2), (-8.0, 0.8), (-12.0, 1.2)]:
    A, B, o = run(l, a)
    print(f"=== l={l} a={a} ===")
    for name, (sol, nrm) in o.items():
        e = sol.y[:, -1]
        tgt = B if name.startswith("A") else A
        d = np.linalg.norm(sol.y - tgt[:, None], axis=0)
        r = np.hypot(*(e[:2]/e[2])) if abs(e[2]) > 1e-12 else np.inf
        print(f" {name}: |s|dev={abs(nrm-1).max():.2e}  end s3={e[2]:+.6f} finite-r={r:.4g}"
              f"  min-dist-to-other-saddle={d.min():.3e}")
