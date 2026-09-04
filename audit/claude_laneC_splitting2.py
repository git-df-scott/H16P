#!/usr/bin/env python3
"""Lane C: robust separatrix-splitting function on the stratum m=5a, b=3l+5.
Saddle S (the one whose stable branch enters from the origin region).
Section Sigma: the line through q = S + 0.15*vs (a point on the local stable
eigendirection, origin side) orthogonal to vs. Integrate the unstable branch
of S that winds around the origin until it crosses Sigma; the splitting is the
signed coordinate of the crossing along the direction orthogonal to vs (with
the true stable branch crossing at ~0 by construction). Zero => loop."""
import numpy as np, sympy as S
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
x, y = S.symbols('x y')
def saddles(l, a):
    m, b = 5*a, 3*l+5
    P = -y+l*x**2+m*x*y+y**2; Q = x+a*x**2+b*x*y
    J = S.Matrix([[S.diff(P,x),S.diff(P,y)],[S.diff(Q,x),S.diff(Q,y)]])
    out = []
    for s in S.solve([P,Q],[x,y],dict=True):
        px, py = complex(s[x]), complex(s[y])
        if abs(px.imag)>1e-9 or abs(py.imag)>1e-9: continue
        Jn = np.array(J.subs(s).evalf(), dtype=float)
        if np.linalg.det(Jn) < 0: out.append((np.array([px.real, py.real]), Jn))
    return out
def field(l, a):
    m, b = 5*a, 3*l+5
    return lambda t, u: [-u[1]+l*u[0]**2+m*u[0]*u[1]+u[1]**2, u[0]+a*u[0]**2+b*u[0]*u[1]]
def split(l, a, verbose=False):
    f = field(l, a)
    for pt, Jn in saddles(l, a):
        if np.linalg.norm(pt) < 1e-9: continue
        w, v = np.linalg.eig(Jn); w = w.real; v = v.real
        vu = v[:, np.argmax(w)]/np.linalg.norm(v[:, np.argmax(w)]); vs = v[:, np.argmin(w)]/np.linalg.norm(v[:, np.argmin(w)])
        # choose the stable direction pointing toward the origin side: integrate backward a bit and see which side winds
        for ss in (+1, -1):
            q = pt+0.15*ss*vs
            n = np.array([-vs[1], vs[0]])  # normal to vs
            for su in (+1, -1):
                u0 = pt+1e-7*su*vu
                def ev(t, u): return np.dot(u-q, vs)   # crossing the line through q orthogonal... we want crossing of Sigma: plane {u: dot(u-q, vs)=0}
                ev.terminal = True; ev.direction = 0
                sol = solve_ivp(f, (0, 100), u0, rtol=1e-11, atol=1e-14, max_step=0.02, events=ev)
                if len(sol.t_events[0]) == 0: continue
                if np.linalg.norm(sol.y[:, -1]) > 50: continue
                hit = sol.y_events[0][0]
                # require the crossing to be near q (within 1.0) and the path to have wound around the origin
                if np.linalg.norm(hit-q) > 1.0: continue
                ang = np.unwrap(np.arctan2(sol.y[1], sol.y[0])); wind = (ang[-1]-ang[0])/(2*np.pi)
                if abs(wind) < 0.6: continue
                d = np.dot(hit-q, n)
                if verbose: print(f"      saddle {np.round(pt,4)} ss={ss} su={su} wind={wind:+.2f} splitting={d:+.5f}")
                return d, pt, wind
    return None
if __name__ == "__main__":
    for a in (1.2, 1.5, 2.0, 3.0):
        lo = -1-np.sqrt(1+3*a*a)
        ls = np.linspace(lo+0.02, -0.3, 30); vals = []
        for l in ls:
            r = split(l, a); vals.append(None if r is None else r[0])
        print(f"a={a}: l from {ls[0]:.3f} to {ls[-1]:.3f}: " + " ".join("  .  " if v is None else f"{v:+.3f}" for v in vals))
        pairs = [(ls[i], ls[i+1]) for i in range(len(ls)-1) if vals[i] is not None and vals[i+1] is not None and vals[i]*vals[i+1] < 0]
        for l1, l2 in pairs:
            lstar = brentq(lambda l: split(l, a)[0], l1, l2, xtol=1e-10)
            d, pt, wind = split(lstar, a, verbose=True)
            eta3 = -25*a*(2*a*a+lstar+2)*(5*a*a*lstar+6*a*a-3*lstar**3-12*lstar**2-15*lstar-6)/64
            m, b = 5*a, 3*lstar+5
            print(f"   LOOP at l*={lstar:.10f} (m={m}, b={b:.8f}); saddle {np.round(pt,6)}; eta3={eta3:+.5f}; (0,1) det={-(1+b):+.4f} trace={m:+.3f}")
