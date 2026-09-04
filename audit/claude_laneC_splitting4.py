#!/usr/bin/env python3
"""Lane C: precise splitting. On the circle |u-S|=RHO compare the crossing
angle of the true stable branch (integrated backward from S+1e-8*vs) with the
crossing angle of the returning unstable branch. Locate the loop l*(a) and
compare with the center curve eta_3=0 on the stratum:
  5a^2 l + 6a^2 = 3l^3 + 12l^2 + 15l + 6."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from claude_laneC_splitting3 import saddle, field
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
RHO = 0.05
def crossings(l, a):
    sd = saddle(l, a)
    if sd is None: return None
    pt, Jn = sd; f = field(l, a)
    w, v = np.linalg.eig(Jn); w = w.real; v = v.real
    vu = v[:, np.argmax(w)]/np.linalg.norm(v[:, np.argmax(w)]); vs = v[:, np.argmin(w)]/np.linalg.norm(v[:, np.argmin(w)])
    if np.linalg.norm(pt+RHO*vs) > np.linalg.norm(pt-RHO*vs): vs = -vs
    def ev_out(t, u): return np.hypot(*(u-pt))-RHO
    ev_out.direction = +1; ev_out.terminal = True
    # stable branch backward: leaves the circle
    sb = solve_ivp(f, (0, -20), pt+1e-8*vs, rtol=1e-12, atol=1e-15, max_step=0.005, events=ev_out)
    es = (sb.y_events[0][0]-pt)/RHO
    def ev_in(t, u): return np.hypot(*(u-pt))-RHO
    ev_in.direction = -1; ev_in.terminal = True
    for su in (+1, -1):
        ub = solve_ivp(f, (0, 300), pt+1e-8*su*vu, rtol=1e-12, atol=1e-15, max_step=0.005, events=ev_in)
        if len(ub.t_events[0]) == 0 or np.linalg.norm(ub.y[:, -1]) > 100: continue
        ang = np.unwrap(np.arctan2(ub.y[1], ub.y[0])); wind = (ang[-1]-ang[0])/(2*np.pi)
        if abs(wind) < 0.5: continue
        eu = (ub.y_events[0][0]-pt)/RHO
        if np.dot(eu, es) < 0: continue
        return es[0]*eu[1]-es[1]*eu[0]   # sine of angle from stable crossing to unstable crossing
    return None
def center_l(a):
    r = np.roots([-3, -12, 5*a*a-15, 6*a*a-6]); r = r[np.isreal(r)].real
    return r
for a in (1.0, 1.5, 2.0, 3.0):
    lc = [v for v in center_l(a) if -1-np.sqrt(1+3*a*a) < v < -1+np.sqrt(1+3*a*a)]
    # bracket near the previous estimate
    l0 = -1.19
    ls = np.linspace(l0-0.15, l0+0.15, 13); vals = [crossings(l, a) for l in ls]
    br = [(ls[i], ls[i+1]) for i in range(12) if vals[i] is not None and vals[i+1] is not None and vals[i]*vals[i+1] < 0]
    if not br: print(f"a={a}: no bracket; vals={vals}"); continue
    lstar = brentq(lambda l: crossings(l, a), *br[0], xtol=1e-12)
    eta3 = -25*a*(2*a*a+lstar+2)*(5*a*a*lstar+6*a*a-3*lstar**3-12*lstar**2-15*lstar-6)/64
    sd = saddle(lstar, a)
    print(f"a={a}: loop l*={lstar:.12f}  center-curve l_c={[round(v,12) for v in lc]}  l*-l_c={lstar-min(lc, key=lambda v: abs(v-lstar)):+.3e}  eta3(l*)={eta3:+.3e}  trace(saddle)={np.trace(sd[1]):+.3e}")
