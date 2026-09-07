#!/usr/bin/env python3
"""Identify the boundary of the origin nest and which infinite singularities it
limits on, at a point of the order-2 + neutral-graphic family."""
import sys, numpy as np
from scipy.integrate import solve_ivp
from sphere_gen import make, inf_points

def params(lv, av):
    mv = av*(av**2 + 3*lv + 1)/(av**2 + lv + 1)
    return mv, -2*lv + mv*(lv + 1)/av

def returns(l, m, a, b, x0, T=80.0):
    F = make(0.0, l, m, a, b, sgn=+1)
    s = np.array([x0, 0.0, 1.0]); s = s/np.linalg.norm(s)
    def ev(t, y): return y[1]
    ev.direction = +1
    sol = solve_ivp(F, (0, T), s, rtol=1e-10, atol=1e-12, events=ev)
    for t, p in zip(sol.t_events[0], sol.y_events[0]):
        if t > 1e-7 and p[0] > 0 and p[2] > 1e-12: return True
    return False

def follow(l, m, a, b, x0, sgn, pts, T=150.0):
    F = make(0.0, l, m, a, b, sgn=sgn)
    s = np.array([x0, 0.0, 1.0]); s = s/np.linalg.norm(s)
    sol = solve_ivp(F, (0, T), s, rtol=1e-11, atol=1e-13)
    res = sorted((np.linalg.norm(sol.y - sv[:, None], axis=0).min(), u, kind, sv[0] > 0)
                 for (u, sv, kind, le, lt) in pts)
    return res[:2], sol.y[:, -1]

lv, av = float(sys.argv[1]), float(sys.argv[2])
mv, bv = params(lv, av)
print("l=%g a=%g m=%.10g b=%.10g" % (lv, av, mv, bv), flush=True)
lo, hi = 0.02, 1.5
if not returns(lv, mv, av, bv, lo):
    print("no inner return"); sys.exit()
for i in range(28):
    mid = 0.5*(lo+hi)
    if returns(lv, mv, av, bv, mid): lo = mid
    else: hi = mid
print("nest boundary crosses y=0 at x* = %.10f" % lo, flush=True)
pts = inf_points(lv, mv, av, bv)
for xt, tag in ((lo*(1-1e-7), "inside "), (hi*(1+1e-7), "outside")):
    for sgn, nm in ((+1, "fwd"), (-1, "bwd")):
        top, e = follow(lv, mv, av, bv, xt, sgn, pts)
        s = "  ".join("u=%+.5f(%s,%s) d=%.2e" % (u, k, "+" if p else "-", d) for d, u, k, p in top)
        xy = e[:2]/e[2] if e[2] > 1e-9 else None
        print("  %s %s: %s | end s3=%.2e xy=%s"
              % (tag, nm, s, e[2], "inf" if xy is None else np.round(xy, 3)), flush=True)
