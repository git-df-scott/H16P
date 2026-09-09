#!/usr/bin/env python3
"""Identify the boundary of the origin's period annulus/nest: bisect the last
returning point on {y=0,x>0}, then follow the boundary orbit forward and
backward and record which infinite singularities it limits on."""
import numpy as np
from scipy.integrate import solve_ivp
from sphere_gen import make, inf_points

def params(lv, av):
    mv = av*(av**2 + 3*lv + 1)/(av**2 + lv + 1)
    return mv, -2*lv + mv*(lv + 1)/av

def returns(lam, l, m, a, b, x0, T=200.0):
    F = make(lam, l, m, a, b, sgn=+1)
    s = np.array([x0, 0.0, 1.0]); s = s/np.linalg.norm(s)
    def ev(t, y): return y[1]
    ev.direction = +1
    def esc(t, y): return y[2] - 1e-11
    esc.terminal = True; esc.direction = -1
    sol = solve_ivp(F, (0, T), s, rtol=1e-12, atol=1e-14, events=[ev, esc])
    for t, p in zip(sol.t_events[0], sol.y_events[0]):
        if t > 1e-7 and p[0] > 0 and p[2] > 1e-12:
            return True
    return False

def limits(lam, l, m, a, b, x0, pts, T=300.0):
    out = {}
    for sgn in (+1, -1):
        F = make(lam, l, m, a, b, sgn=sgn)
        s = np.array([x0, 0.0, 1.0]); s = s/np.linalg.norm(s)
        sol = solve_ivp(F, (0, T), s, rtol=1e-12, atol=1e-14)
        best = None
        for (u, sv, kind, le, lt) in pts:
            d = np.linalg.norm(sol.y - sv[:, None], axis=0).min()
            if best is None or d < best[0]: best = (d, u, kind, sv)
        e = sol.y[:, -1]
        out["fwd" if sgn > 0 else "bwd"] = (best, e)
    return out

for (lv, av) in [(-3.0, 0.4), (-6.0, 0.4), (-12.0, 1.0)]:
    mv, bv = params(lv, av)
    lo, hi = 0.05, 2.0
    if not returns(0, lv, mv, av, bv, lo): print("l=%g a=%g: no inner return" % (lv, av)); continue
    for _ in range(50):
        mid = 0.5*(lo + hi)
        if returns(0, lv, mv, av, bv, mid): lo = mid
        else: hi = mid
    print("=== l=%g a=%g m=%.8g b=%.8g ===" % (lv, av, mv, bv))
    print("  nest boundary crosses y=0 at x* = %.12f" % lo)
    pts = inf_points(lv, mv, av, bv)
    for xt, tag in ((lo*(1-1e-9), "just inside"), (hi*(1+1e-9), "just outside")):
        L = limits(0, lv, mv, av, bv, xt, pts)
        for k, (best, e) in L.items():
            d, u, kind, sv = best
            xy = e[:2]/e[2] if e[2] > 1e-9 else None
            print("   %-12s %s: nearest infinite sing u=%+.6f (%s) at distance %.3e ; end s3=%.2e xy=%s"
                  % (tag, k, u, kind, d, e[2], "inf" if xy is None else np.round(xy, 3)))
    print()
