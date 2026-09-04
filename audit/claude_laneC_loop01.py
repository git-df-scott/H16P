#!/usr/bin/env python3
"""Lane C: loop through the saddle (0,1) around the origin on the stratum
(exists as a saddle iff l>-2, trace = 5a != 0 so hyperbolic). Precise
splitting: crossing angles of the true stable branch (backward) and the
returning unstable branch on a circle of radius RHO around (0,1)."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
RHO = 0.05
def field(l, a):
    m, b = 5*a, 3*l+5
    return lambda t, u: [-u[1]+l*u[0]**2+m*u[0]*u[1]+u[1]**2, u[0]+a*u[0]**2+b*u[0]*u[1]]
def sigma(l, a):
    m, b = 5*a, 3*l+5
    pt = np.array([0.0, 1.0]); Jn = np.array([[m, 1.0], [1.0+b, 0.0]])
    if np.linalg.det(Jn) >= 0: return None
    f = field(l, a)
    w, v = np.linalg.eig(Jn); w = w.real; v = v.real
    vu = v[:, np.argmax(w)]/np.linalg.norm(v[:, np.argmax(w)]); vs = v[:, np.argmin(w)]/np.linalg.norm(v[:, np.argmin(w)])
    def ev_out(t, u): return np.hypot(*(u-pt))-RHO
    ev_out.direction = +1; ev_out.terminal = True
    def ev_in(t, u): return np.hypot(*(u-pt))-RHO
    ev_in.direction = -1; ev_in.terminal = True
    best = None
    for ss in (+1, -1):
        sb = solve_ivp(f, (0, -30), pt+1e-8*ss*vs, rtol=1e-12, atol=1e-15, max_step=0.005, events=ev_out)
        if len(sb.t_events[0]) == 0: continue
        es = (sb.y_events[0][0]-pt)/RHO
        for su in (+1, -1):
            ub = solve_ivp(f, (0, 300), pt+1e-8*su*vu, rtol=1e-12, atol=1e-15, max_step=0.005, events=ev_in)
            if len(ub.t_events[0]) == 0 or np.linalg.norm(ub.y[:, -1]) > 100: continue
            ang = np.unwrap(np.arctan2(ub.y[1], ub.y[0])); wind = (ang[-1]-ang[0])/(2*np.pi)
            if abs(wind) < 0.5: continue          # must wind around the origin
            eu = (ub.y_events[0][0]-pt)/RHO
            if np.dot(eu, es) < 0.3: continue      # return on the same side as the stable branch
            best = (es[0]*eu[1]-es[1]*eu[0], wind, ss, su)
    return best
for a in (0.5, 0.7, 1.0, 1.5, 2.0):
    hi = -1+np.sqrt(1+3*a*a)
    ls = np.linspace(-1.95, hi+1.0, 24)
    vals = [sigma(l, a) for l in ls]
    print(f"a={a} (finite saddles x<0 exist for l<{hi:.3f}): " + " ".join("  .  " if v is None else f"{v[0]:+.3f}" for v in vals))
    for i in range(len(ls)-1):
        if vals[i] is not None and vals[i+1] is not None and vals[i][0]*vals[i+1][0] < 0:
            lstar = brentq(lambda l: sigma(l, a)[0], ls[i], ls[i+1], xtol=1e-10)
            eta3 = -25*a*(2*a*a+lstar+2)*(5*a*a*lstar+6*a*a-3*lstar**3-12*lstar**2-15*lstar-6)/64
            print(f"   LOOP through (0,1) at l*={lstar:.10f}: eta3={eta3:+.5f}  saddle trace={5*a}  winding={sigma(lstar,a)[1]:+.2f}")
