#!/usr/bin/env python3
"""Lane C: continuous separatrix-splitting on the stratum m=5a, b=3l+5.
Saddle S: the finite saddle with x<0 (continuous in (l,a) in the region
3a^2>l^2+2l). Orientation: vs points from S toward the origin side.
Splitting sigma(l,a): the unstable branch that winds around the origin is
integrated until it re-enters the circle |u-S|=rho from outside; sigma is
the sine of the angle between vs and (hit-S), i.e. the signed transverse
offset of the returning branch from the stable direction. A sign change in
l at fixed a locates a homoclinic loop around the third-order weak focus."""
import numpy as np, sympy as S
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
x, y = S.symbols('x y')
RHO = 0.08
def saddle(l, a):
    m, b = 5*a, 3*l+5
    P = -y+l*x**2+m*x*y+y**2; Q = x+a*x**2+b*x*y
    J = S.Matrix([[S.diff(P,x),S.diff(P,y)],[S.diff(Q,x),S.diff(Q,y)]])
    cands = []
    for s in S.solve([P,Q],[x,y],dict=True):
        px, py = complex(s[x]), complex(s[y])
        if abs(px.imag)>1e-9 or abs(py.imag)>1e-9: continue
        Jn = np.array(J.subs(s).evalf(), dtype=float)
        if np.linalg.det(Jn) < 0 and px.real < -1e-6: cands.append((np.array([px.real, py.real]), Jn))
    return cands[0] if cands else None
def field(l, a):
    m, b = 5*a, 3*l+5
    return lambda t, u: [-u[1]+l*u[0]**2+m*u[0]*u[1]+u[1]**2, u[0]+a*u[0]**2+b*u[0]*u[1]]
def sigma(l, a, detail=False):
    sd = saddle(l, a)
    if sd is None: return None
    pt, Jn = sd; f = field(l, a)
    w, v = np.linalg.eig(Jn); w = w.real; v = v.real
    vu = v[:, np.argmax(w)]/np.linalg.norm(v[:, np.argmax(w)]); vs = v[:, np.argmin(w)]/np.linalg.norm(v[:, np.argmin(w)])
    if np.linalg.norm(pt+RHO*vs) > np.linalg.norm(pt-RHO*vs): vs = -vs   # toward origin side
    def ev(t, u): return np.hypot(*(u-pt))-RHO
    ev.direction = -1; ev.terminal = True
    best = None
    for su in (+1, -1):
        u0 = pt+1e-7*su*vu
        # leave the circle first
        sol1 = solve_ivp(f, (0, 5), u0, rtol=1e-11, atol=1e-14, max_step=0.01, events=lambda t,u: np.hypot(*(u-pt))-RHO*1.01)
        sol = solve_ivp(f, (0, 200), u0, rtol=1e-11, atol=1e-14, max_step=0.01, events=ev)
        if len(sol.t_events[0]) == 0 or np.linalg.norm(sol.y[:, -1]) > 100: continue
        ang = np.unwrap(np.arctan2(sol.y[1], sol.y[0])); wind = (ang[-1]-ang[0])/(2*np.pi)
        if abs(wind) < 0.5: continue
        hit = sol.y_events[0][0]; e = (hit-pt)/RHO
        s = vs[0]*e[1]-vs[1]*e[0]   # sine of angle from vs to e
        if np.dot(e, vs) < 0: continue  # returned from the wrong side
        best = (s, su, wind, pt)
    if detail and best: print("      ", best)
    return None if best is None else best[0]
if __name__ == "__main__":
    for a in (1.0, 1.2, 1.5, 2.0, 2.5, 3.0):
        lo = -1-np.sqrt(1+3*a*a); ls = np.linspace(lo+0.02, -0.4, 28)
        vals = [sigma(l, a) for l in ls]
        print(f"a={a}: l in [{ls[0]:.2f},{ls[-1]:.2f}]: " + " ".join("  .  " if v is None else f"{v:+.3f}" for v in vals))
        for i in range(len(ls)-1):
            if vals[i] is not None and vals[i+1] is not None and vals[i]*vals[i+1] < 0:
                lstar = brentq(lambda l: sigma(l, a), ls[i], ls[i+1], xtol=1e-9)
                sd = saddle(lstar, a); eta3 = -25*a*(2*a*a+lstar+2)*(5*a*a*lstar+6*a*a-3*lstar**3-12*lstar**2-15*lstar-6)/64
                print(f"   LOOP CANDIDATE l*={lstar:.9f} m={5*a} b={3*lstar+5:.7f} saddle={np.round(sd[0],5)} trace(saddle)={np.trace(sd[1]):+.5f} eta3={eta3:+.4f}")
