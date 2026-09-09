#!/usr/bin/env python3
"""Where do the separatrices of the two infinite saddles go?
A = (+dir) has equator-stable / transverse-unstable: one orbit leaves into s3>0.
B = (-dir) has equator-unstable / transverse-stable: one orbit enters from s3>0.
If A's unstable branch has omega-limit B, the boundary graphic of the origin
nest is [that connection] + [an equator arc B->A]."""
import numpy as np
from scipy.integrate import solve_ivp
from sphere import make, jac, inf_points

def branch(l, a, s0, forward, T=400.0):
    F = make(l, a)
    J = jac(l, a, s0); w, v = np.linalg.eig(J); w = w.real; v = v.real
    # tangent eigenvectors: drop the (radial) one most aligned with s0
    idx = [i for i in range(3)]
    rad = max(idx, key=lambda i: abs(np.dot(v[:, i]/np.linalg.norm(v[:, i]), s0)))
    tang = [i for i in idx if i != rad]
    sgn = +1 if forward else -1
    # pick the eigendirection that leaves the equator (s3 component nonzero)
    k = max(tang, key=lambda i: abs(v[2, i]))
    d = v[:, k]/np.linalg.norm(v[:, k])
    if d[2] < 0: d = -d                     # into the upper hemisphere s3>0
    out = []
    for eps in (1e-7,):
        s = s0 + eps*d; s = s/np.linalg.norm(s)
        sol = solve_ivp(F, (0, sgn*T), s, rtol=1e-12, atol=1e-14, dense_output=False,
                        max_step=0.02)
        out.append(sol)
    return out[0], w[tang], v[:, tang]

def report(l, a):
    print(f"=== l={l} a={a} ===")
    pts = inf_points(l, a)
    for uu, s0 in pts:
        J = jac(l, a, s0); w = np.linalg.eigvals(J).real
        sol, wt, vt = branch(l, a, s0, forward=(s0[0] > 0))
        end = sol.y[:, -1]
        xy = end[:2]/end[2] if abs(end[2]) > 1e-12 else None
        print(f" saddle at s={np.round(s0,4)} tangent eigs {np.round(wt,4)}")
        print(f"   integrated {'fwd' if s0[0]>0 else 'bwd'} to t={sol.t[-1]:.2f}: s={np.round(end,6)}"
              f"  (x,y)={None if xy is None else np.round(xy,6)}  |s3|={abs(end[2]):.3e}")
        # closest approach to the other infinite point
        other = [p for _, p in pts if np.dot(p, s0) < 0][0]
        d = np.linalg.norm(sol.y - other[:, None], axis=0)
        print(f"   min distance to antipodal saddle along the branch: {d.min():.3e} at t={sol.t[d.argmin()]:.2f}")
        # closest approach to (0,1)
        n = np.sqrt(2.0); p01 = np.array([0.0, 1.0, 1.0])/n
        d2 = np.linalg.norm(sol.y - p01[:, None], axis=0)
        print(f"   min distance to the focus (0,1):                   {d2.min():.3e}")

for (l, a) in [(-10.0, 1.0), (-12.0, 0.8), (-8.0, 1.2)]:
    report(l, a)
