#!/usr/bin/env python3
"""Sequence of crossings of the section {y=0, x>0} on the Poincare sphere,
for (i) the separatrix leaving infinite saddle A, (ii) the separatrix entering
its antipode B (backward), (iii) a small circle around the weak focus."""
import numpy as np
from scipy.integrate import solve_ivp
from sphere3 import make, inf_points, tangent_eigs

def crossings(l, a, s_start, sgn, T=4000.0, nmax=400):
    F = make(l, a, sgn=sgn)
    def ev(t, y): return y[1]
    ev.direction = +1 if sgn > 0 else -1
    sol = solve_ivp(F, (0, T), s_start, rtol=1e-12, atol=1e-14, events=ev,
                    dense_output=False)
    xs = []
    for p in sol.y_events[0]:
        if p[0] > 0 and p[2] > 1e-12:
            xs.append(p[0]/p[2])            # finite x at the crossing
    e = sol.y[:, -1]
    return np.array(xs[:nmax]), e, sol.t[-1]

def sep_start(l, a, which, eps=1e-9):
    pts = inf_points(l, a)
    A = [p for u, p in pts if p[0] > 0][0]; B = [p for u, p in pts if p[0] < 0][0]
    s0, sgn = (A, +1) if which == "A_out" else (B, -1)
    w, v = tangent_eigs(l, a, s0)
    k = int(np.argmax(np.abs(v[2, :]))); d = v[:, k]/np.linalg.norm(v[:, k])
    if d[2] < 0: d = -d
    s = s0 + eps*d
    return s/np.linalg.norm(s), sgn

def eta3(l, a):
    return -25*a*(2*a**2+l+2)*(5*a**2*l+6*a**2-3*l**3-12*l**2-15*l-6)/64

for (l, a) in [(-10.0, 1.0), (-12.0, 0.8), (-8.0, 1.2)]:
    print(f"=== l={l} a={a}   eta3={eta3(l,a):+.4g} ===")
    for which in ("A_out", "B_in"):
        s, sgn = sep_start(l, a, which)
        xs, e, T = crossings(l, a, s, sgn)
        head = np.round(xs[:6], 6); tail = np.round(xs[-4:], 8)
        print(f"  {which}: {len(xs)} crossings of y=0,x>0; first {head} ... last {tail}")
    # small orbit around the origin, forward
    x0 = 0.02; s = np.array([x0, 0.0, 1.0]); s /= np.linalg.norm(s)
    xs, e, T = crossings(l, a, s, +1, T=20000.0)
    print(f"  origin r=0.02 fwd: x-crossings {np.round(xs[:3],8)} ... {np.round(xs[-3:],8)}"
          f"  ({'outward=unstable' if len(xs)>2 and xs[-1]>xs[0] else 'inward=stable'})")
