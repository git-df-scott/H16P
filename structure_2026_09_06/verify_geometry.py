#!/usr/bin/env python3
"""(1) r=1 recomputed from the sphere Jacobian, independent of the chart formulas.
(2) At the collision point l=-10, a=2, m=10, b=-25 the origin is a centre, so it
    has a genuine period annulus. Find its outer boundary and check that the
    boundary orbit limits on BOTH infinite saddles (i.e. is the two-saddle graphic)."""
import numpy as np
from scipy.integrate import solve_ivp
from sphere_gen import make, inf_points, tangent_eigs

lv, av, mv, bv = -10.0, 2.0, 10.0, -25.0
print("=== (1) neutrality from the sphere Jacobian ===")
pts = inf_points(lv, mv, av, bv)
rows = []
for (u, sv, kind, le, lt) in pts:
    if kind != "saddle": continue
    w, _ = tangent_eigs(0.0, lv, mv, av, bv, sv)
    w = np.sort(w)
    rows.append((u, sv[0] > 0, w[0], w[1], abs(w[0]/w[1])))
    print("  u=%+.8f s1%s  sphere eigs = (%+.8f, %+.8f)  |ratio| = %.10f"
          % (u, ">0" if sv[0] > 0 else "<0", w[0], w[1], abs(w[0]/w[1])))
pos = [r for r in rows if r[1]]
print("  product over the two non-antipodal saddles (one from each pair):")
# The graphic runs equator -> S1 -> plane -> S2 -> equator, so it uses the two
# saddles on the SAME side, not a saddle and the antipode of the other:
# at S1 the ratio is |lam_eq/lam_tr|, at S2 it is the reciprocal.
r1 = [r for r in rows if r[0] < 0 and r[1]][0]
r2 = [r for r in rows if r[0] > 0 and r[1]][0]
print("    |ratio| at the two saddle directions: %.10f and %.10f" % (r1[4], r2[4]))
print("    r = %.12f  (chart formula predicts exactly 1)" % (r1[4]*(1/r2[4])))

print("\n=== (2) the period annulus of the centre and its boundary ===")
F = make(0.0, lv, mv, av, bv, sgn=+1)
def closes(x0, T=200.0):
    s = np.array([x0, 0.0, 1.0]); s = s/np.linalg.norm(s)
    def ev(t, y): return y[1]
    ev.direction = +1
    sol = solve_ivp(F, (0, T), s, rtol=1e-13, atol=1e-15, events=ev)
    for t, p in zip(sol.t_events[0], sol.y_events[0]):
        if t > 1e-7 and p[0] > 0 and p[2] > 1e-12:
            return abs(p[0]/p[2] - x0) < 1e-7, p[0]/p[2] - x0
    return False, None
for x0 in (0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4):
    ok, d = closes(x0)
    print("   x0=%-6g closed orbit? %-5s  displacement=%s" % (x0, ok, "n/a" if d is None else "%+.2e" % d))
lo, hi = 0.05, 0.6
for _ in range(40):
    mid = 0.5*(lo+hi)
    if closes(mid)[0]: lo = mid
    else: hi = mid
print("   annulus boundary crosses y=0 at x* = %.10f" % lo)
for sgn in (+1, -1):
    G = make(0.0, lv, mv, av, bv, sgn=sgn)
    s = np.array([lo*(1-1e-8), 0.0, 1.0]); s = s/np.linalg.norm(s)
    sol = solve_ivp(G, (0, 200.0), s, rtol=1e-13, atol=1e-15)
    res = sorted((np.linalg.norm(sol.y - sv[:, None], axis=0).min(), u, kind)
                 for (u, sv, kind, le, lt) in pts)
    print("   boundary orbit %s: closest approaches -> %s"
          % ("fwd" if sgn > 0 else "bwd",
             ", ".join("u=%+.4f(%s) d=%.2e" % (u, k, d) for d, u, k in res[:3])))
