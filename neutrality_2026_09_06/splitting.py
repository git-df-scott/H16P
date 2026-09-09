#!/usr/bin/env python3
"""Splitting function for the candidate boundary graphic.
A = infinite saddle whose transverse separatrix leaves into the plane;
B = its antipode, whose transverse separatrix enters from the plane.
Both are integrated to their FIRST crossing of the section {s1=0, s2>0}
(the positive y-axis, lifted to the sphere).  The connection A->B exists
exactly where the two crossing points coincide."""
import numpy as np
from scipy.integrate import solve_ivp
from sphere3 import make, inf_points, tangent_eigs

def start(l, a, which, eps=1e-9):
    pts = inf_points(l, a)
    A = [p for u, p in pts if p[0] > 0][0]; B = [p for u, p in pts if p[0] < 0][0]
    s0, sgn = (A, +1) if which == "A" else (B, -1)
    w, v = tangent_eigs(l, a, s0)
    k = int(np.argmax(np.abs(v[2, :]))); d = v[:, k]/np.linalg.norm(v[:, k])
    if d[2] < 0: d = -d
    s = s0 + eps*d
    return s/np.linalg.norm(s), sgn

def cross(l, a, which, T=60.0):
    s, sgn = start(l, a, which)
    F = make(l, a, sgn=sgn)
    def ev(t, y): return y[0]
    ev.terminal = True
    sol = solve_ivp(F, (0, T), s, rtol=1e-12, atol=1e-14, events=ev)
    if len(sol.t_events[0]) == 0: return None, sol
    p = sol.y_events[0][0]
    return (p[1], p[2]), sol           # (s2,s3) on the s1=0 section

def split(l, a):
    ca, _ = cross(l, a, "A"); cb, _ = cross(l, a, "B")
    if ca is None or cb is None: return None
    # compare on the section by the finite ordinate y = s2/s3 (or 1/y near infinity)
    ya = ca[1]/ca[0] if abs(ca[0]) > 1e-14 else 0.0     # 1/y, regular at infinity
    yb = cb[1]/cb[0] if abs(cb[0]) > 1e-14 else 0.0
    return ya - yb, ya, yb

print(" l      a      1/y(A_out)   1/y(B_in)    splitting")
for l in (-12.0, -11.0, -10.0, -9.0, -8.0):
    for a in (0.8, 1.0, 1.2):
        r = split(l, a)
        print("%6.2f %5.2f   " % (l, a) +
              ("no section crossing" if r is None
               else "%+.9f  %+.9f  %+.3e" % (r[1], r[2], r[0])))
