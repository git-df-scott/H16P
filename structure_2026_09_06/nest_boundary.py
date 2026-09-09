#!/usr/bin/env python3
"""Where does the origin's period annulus end, and what bounds it?
Return map on the ray {y=0, x>0}, plus the fate of the transverse separatrix
of each infinite saddle."""
import numpy as np
from scipy.integrate import solve_ivp
from sphere_gen import make, inf_points, tangent_eigs

def params(lv, av):
    mv = av*(av**2 + 3*lv + 1)/(av**2 + lv + 1)
    return mv, -2*lv + mv*(lv + 1)/av

def first_return(lam, l, m, a, b, x0, T=120.0):
    F = make(lam, l, m, a, b, sgn=+1)
    s = np.array([x0, 0.0, 1.0]); s = s/np.linalg.norm(s)
    def ev(t, y): return y[1]
    ev.direction = +1
    def esc(t, y): return y[2] - 1e-11
    esc.terminal = True; esc.direction = -1
    sol = solve_ivp(F, (0, T), s, rtol=1e-12, atol=1e-14, events=[ev, esc])
    for t, p in zip(sol.t_events[0], sol.y_events[0]):
        if t > 1e-7 and p[0] > 0 and p[2] > 1e-12:
            return p[0]/p[2], t
    return None, None

def sep_fate(lam, l, m, a, b, s0, sgn, T=200.0):
    w, v = tangent_eigs(lam, l, m, a, b, s0)
    k = int(np.argmax(np.abs(v[2, :]))); d = v[:, k]/np.linalg.norm(v[:, k])
    if d[2] < 0: d = -d
    s = s0 + 1e-9*d; s /= np.linalg.norm(s)
    F = make(lam, l, m, a, b, sgn=sgn)
    sol = solve_ivp(F, (0, T), s, rtol=1e-12, atol=1e-14)
    return sol

for (lv, av) in [(-3.0, 0.4), (-6.0, 0.4), (-12.0, 1.0)]:
    mv, bv = params(lv, av)
    print("=== l=%g a=%g  m=%.8g b=%.8g ===" % (lv, av, mv, bv))
    for x0 in (1e-3, 3e-3, 0.01, 0.03, 0.06, 0.1, 0.2, 0.4, 0.8):
        r, t = first_return(0.0, lv, mv, av, bv, x0)
        if r is None: print("   x0=%-7g NO RETURN" % x0)
        else:         print("   x0=%-7g R=%.10g  d=%+.4e" % (x0, r, r-x0))
    pts = inf_points(lv, mv, av, bv)
    for (u, s0, kind, le, lt) in pts:
        if kind != "saddle": continue
        sgn = +1 if lt > 0 else -1      # follow the branch that lives in s3>0
        sol = sep_fate(0.0, lv, mv, av, bv, s0, sgn)
        e = sol.y[:, -1]
        xy = e[:2]/e[2] if e[2] > 1e-9 else None
        # nearest other infinite singularity
        dmin, uarg = 1e9, None
        for (u2, s2, k2, _, _) in pts:
            if np.allclose(s2, s0): continue
            d = np.linalg.norm(sol.y - s2[:, None], axis=0).min()
            if d < dmin: dmin, uarg = d, (u2, k2)
        print("   saddle u=%+.6f s=(%+.3f,%+.3f) dir=%s: end s3=%.2e xy=%s"
              % (u, s0[0], s0[1], "fwd" if sgn > 0 else "bwd", e[2],
                 "inf" if xy is None else np.round(xy, 4)))
        print("        closest approach %.3e to infinite singularity u=%+.6f (%s)"
              % (dmin, uarg[0], uarg[1]))
    print()
