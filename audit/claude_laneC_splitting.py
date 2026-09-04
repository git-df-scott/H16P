#!/usr/bin/env python3
"""Lane C: signed separatrix splitting on the third-order weak-focus stratum
in the finite-saddle region 3a^2 > l^2+2l. For the saddle whose stable
branch comes from the origin nest, follow the unstable branch that winds
around the origin and measure, at its return to the saddle's stable
eigendirection line, the signed offset. A sign change along a path in (l,a)
locates a homoclinic loop surrounding the third-order weak focus."""
import numpy as np, sympy as S
from scipy.integrate import solve_ivp
x, y = S.symbols('x y')
def data(l, a):
    m, b = 5*a, 3*l+5
    P = -y+l*x**2+m*x*y+y**2; Q = x+a*x**2+b*x*y
    J = S.Matrix([[S.diff(P,x),S.diff(P,y)],[S.diff(Q,x),S.diff(Q,y)]])
    f = lambda t, u: [-u[1]+l*u[0]**2+m*u[0]*u[1]+u[1]**2, u[0]+a*u[0]**2+b*u[0]*u[1]]
    sad = []
    for s in S.solve([P,Q],[x,y],dict=True):
        px, py = complex(s[x]), complex(s[y])
        if abs(px.imag)>1e-9 or abs(py.imag)>1e-9: continue
        Jn = np.array(J.subs(s).evalf(), dtype=float)
        if np.linalg.det(Jn) < 0: sad.append((np.array([px.real, py.real]), Jn))
    return f, sad
def splitting(l, a):
    f, sad = data(l, a)
    best = None
    for pt, Jn in sad:
        w, v = np.linalg.eig(Jn); w = w.real; v = v.real
        vu = v[:, np.argmax(w)]; vs = v[:, np.argmin(w)]; vu /= np.linalg.norm(vu); vs /= np.linalg.norm(vs)
        for sgn in (+1, -1):
            u0 = pt+sgn*1e-7*vu
            # stop when trajectory returns near the saddle after winding around the origin
            def ev(t, u): return np.hypot(*(u-pt))-0.3
            ev.direction = -1; ev.terminal = False
            sol = solve_ivp(f, (0, 80), u0, rtol=1e-11, atol=1e-14, max_step=0.02, events=ev)
            if len(sol.t_events[0]) == 0: continue
            # find first return time after leaving; compute winding about origin up to that time
            te = sol.t_events[0][0]; mask = sol.t <= te
            ang = np.unwrap(np.arctan2(sol.y[1, mask], sol.y[0, mask])); wind = (ang[-1]-ang[0])/(2*np.pi)
            if abs(wind) < 0.5: continue
            # signed offset of the returning point from the stable eigenline through the saddle, measured along vu
            ret = sol.y_events[0][0]-pt
            off = np.dot(ret, vu)  # component along unstable direction: 0 means on the stable manifold (to first order)
            best = (pt, sgn, off, wind)
    return best
for a in (0.7, 1.0, 1.5):
    print(f"a={a}: finite-saddle l-interval ({-1-np.sqrt(1+3*a*a):.3f},{-1+np.sqrt(1+3*a*a):.3f})")
    for l in np.linspace(-1-np.sqrt(1+3*a*a)+0.15, -1+np.sqrt(1+3*a*a)-0.15, 9):
        r = splitting(l, a)
        if r is None: print(f"   l={l:+.3f}: no winding return"); continue
        pt, sgn, off, wind = r
        print(f"   l={l:+.3f}: saddle {np.round(pt,3)} branch {sgn:+d} winding {wind:+.2f} signed offset {off:+.4e}")
