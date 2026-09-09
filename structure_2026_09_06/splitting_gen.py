#!/usr/bin/env python3
"""Splitting of the two-saddle graphic on the order-2 + neutral family.
A = infinite saddle whose transverse separatrix LEAVES into s3>0.
B = infinite saddle whose transverse separatrix ENTERS from s3>0.
Both are run to the section {s1=0, s2<0 or >0}; the graphic exists where the
two crossing points coincide."""
import numpy as np
from scipy.integrate import solve_ivp
from sphere_gen import make, inf_points, tangent_eigs

def params(lv, av):
    mv = av*(av**2 + 3*lv + 1)/(av**2 + lv + 1)
    return mv, -2*lv + mv*(lv + 1)/av

def branch(l, m, a, b, s0, eps=1e-9, T=60.0):
    """Return (sgn, solution) for the transverse separatrix that lives in s3>0."""
    w, v = tangent_eigs(0.0, l, m, a, b, s0)
    k = int(np.argmax(np.abs(v[2, :]))); d = v[:, k]/np.linalg.norm(v[:, k])
    if d[2] < 0: d = -d
    out = []
    for sgn in (+1, -1):
        F = make(0.0, l, m, a, b, sgn=sgn)
        s = s0 + eps*d; s = s/np.linalg.norm(s)
        def ev(t, y): return y[0]
        ev.terminal = True
        sol = solve_ivp(F, (0, T), s, rtol=1e-12, atol=1e-14, events=ev)
        grew = sol.y[2].max()
        out.append((grew, sgn, sol))
    out.sort(reverse=True)
    return out[0][1], out[0][2]

def split(lv, av):
    mv, bv = params(lv, av)
    pts = inf_points(lv, mv, av, bv)
    sad = [(u, sv) for (u, sv, k, le, lt) in pts if k == "saddle"]
    us = sorted(set(round(u, 9) for u, _ in sad))
    if len(us) != 2: return None
    res = {}
    for u, sv in sad:
        sgn, sol = branch(lv, mv, av, bv, sv)
        if len(sol.t_events[0]) == 0: continue
        p = sol.y_events[0][0]
        key = (round(u, 9), sv[0] > 0)
        res[key] = (sgn, p[1], p[2])       # (direction, s2, s3) on {s1=0}
    return mv, bv, res

for (lv, av) in [(-3.0, 0.4), (-6.0, 0.4), (-12.0, 1.0), (-12.0, 0.4), (-3.0, 1.0)]:
    r = split(lv, av)
    if r is None: print("l=%g a=%g: not two saddle directions" % (lv, av)); continue
    mv, bv, res = r
    print("l=%-6g a=%-5g m=%-11.6g b=%-11.6g" % (lv, av, mv, bv))
    for key, (sgn, s2, s3) in sorted(res.items()):
        yv = s2/s3 if abs(s3) > 1e-13 else np.inf
        print("    saddle u=%+.6f s1%s  branch=%s : hits {x=0} at y=%+.8g   (s3=%.3e)"
              % (key[0], ">0" if key[1] else "<0", "fwd" if sgn > 0 else "bwd", yv, s3))
