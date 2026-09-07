#!/usr/bin/env python3
"""Positive control for the independent ray-shooting cycle detector.

FC.remote_cycles counts cycles in the SECOND nest and has returned 0 everywhere.
A zero from an unvalidated detector is worth nothing, so point the same
scipy-based detector at the FIRST nest, where the binary128 cusp engine has
independently certified the cycle count, and check it reproduces it.
"""
import json
import numpy as np
from scipy.integrate import solve_ivp

ROWS = json.load(open("ledger_opus/unfold.json"))

def cycles_about_A(mu, dmax=3.0, n=40, tend=80.0):
    """Same algorithm as FC.remote_cycles, but on the focus A=(1,-1),
    shooting along the ray y=-1, x>1 (the cusp engine's own section)."""
    a, a20 = float(mu["a"]), float(mu["a20"])
    a11, a01, a10 = float(mu["a11"]), float(mu["a01"]), float(mu["a10"])
    a00 = a01 + a11 - a10 - a20 - a
    def F(t, z):
        x, y = z
        return [1 + x*y, a00 + a10*x + a20*x*x + a01*y + a11*x*y + a*y*y]
    def disp(d):
        z0 = [1.0 + d, -1.0]
        def ev(t, z): return z[1] + 1.0
        ev.direction = -1
        sol = solve_ivp(F, (0, tend), z0, rtol=1e-9, atol=1e-12, events=ev, max_step=0.05)
        for tt, p in zip(sol.t_events[0], sol.y_events[0]):
            if tt > 1e-6 and p[0] > 1.0: return (p[0] - 1.0) - d
        return None
    ds = [dmax*k/n for k in range(1, n+1)]
    vals = [(d, disp(d)) for d in ds]
    vals = [(d, v) for d, v in vals if v is not None]
    roots = []
    for i in range(len(vals)-1):
        if vals[i][1]*vals[i+1][1] < 0:
            lo, hi = vals[i][0], vals[i+1][0]; f0 = vals[i][1]
            for _ in range(24):
                mid = (lo+hi)/2; fm = disp(mid)
                if fm is None: break
                if fm*f0 > 0: lo = mid
                else: hi = mid
            roots.append((lo+hi)/2)
    return roots, len(vals), (vals[-1][0] if vals else None)

print("POSITIVE CONTROL: same detector, first nest (engine-certified counts)\n")
for r in ROWS[:2]:
    mu = {k: v for k, v in r["mu"].items()}
    roots, nv, dmax = cycles_about_A(mu)
    eng = [float(x) - 1.0 for x in r["local"]]
    print("eta_T=%-6s engine cycles at s-1 = %s" % (r["etaT"], ["%.6f" % v for v in eng]))
    print("            detector  at s-1 = %s   (samples %d, reach %.3f)"
          % (["%.6f" % v for v in roots], nv, dmax if dmax else -1))
    ok = len(roots) == len(eng) and all(min(abs(v-e) for e in eng) < 2e-3 for v in roots)
    print("            -> %s\n" % ("AGREES" if ok else "DISAGREES"))
