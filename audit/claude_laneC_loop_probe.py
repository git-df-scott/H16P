#!/usr/bin/env python3
"""Lane C: at a stratum point with a finite saddle, integrate the saddle's
separatrices and report whether the unstable branch heading toward the
origin nest returns near the saddle (loop candidate) or escapes."""
import numpy as np
from scipy.integrate import solve_ivp
import sympy as S
def system(l, a):
    m, b = 5*a, 3*l+5
    f = lambda t, u: [-u[1]+l*u[0]**2+m*u[0]*u[1]+u[1]**2, u[0]+a*u[0]**2+b*u[0]*u[1]]
    return f, m, b
def equilibria(l, a):
    x, y = S.symbols('x y'); m, b = 5*a, 3*l+5
    P = -y+l*x**2+m*x*y+y**2; Q = x+a*x**2+b*x*y
    J = S.Matrix([[S.diff(P,x),S.diff(P,y)],[S.diff(Q,x),S.diff(Q,y)]])
    out = []
    for s in S.solve([P,Q],[x,y],dict=True):
        px, py = complex(s[x]), complex(s[y])
        if abs(px.imag)>1e-9 or abs(py.imag)>1e-9: continue
        Jn = np.array(J.subs(s).evalf(), dtype=float); w, v = np.linalg.eig(Jn)
        out.append(((px.real, py.real), Jn, w, v))
    return out
for l, a in ((-2.0, 1.0), (-1.0, 1.0), (-3.0, 0.7)):
    f, m, b = system(l, a)
    eqs = equilibria(l, a)
    saddles = [e for e in eqs if np.linalg.det(e[1]) < 0]
    print(f"l={l} a={a} m={m} b={b}: equilibria {[np.round(e[0],4) for e in eqs]}, saddles {[np.round(e[0],4) for e in saddles]}")
    for (pt, Jn, w, v) in saddles:
        for sgn_dir in (+1, -1):
            for which, tdir in (("unstable", 1.0), ("stable", -1.0)):
                idx = int(np.argmax(w.real)) if which == "unstable" else int(np.argmin(w.real))
                vec = v[:, idx].real; vec /= np.linalg.norm(vec)
                u0 = np.array(pt)+sgn_dir*1e-6*vec
                ev = lambda t, u: np.hypot(*(u-np.array(pt)))-1e-6*0.5
                sol = solve_ivp(f, (0, tdir*60), u0, rtol=1e-10, atol=1e-13, max_step=0.05, dense_output=False)
                traj = sol.y; end = traj[:, -1]
                dist_back = np.min(np.hypot(traj[0, 50:]-pt[0], traj[1, 50:]-pt[1])) if traj.shape[1] > 50 else np.inf
                # winding around origin
                ang = np.unwrap(np.arctan2(traj[1], traj[0])); wind = (ang[-1]-ang[0])/(2*np.pi)
                print(f"   saddle {np.round(pt,3)} {which} branch dir {sgn_dir:+d}: end={np.round(end,3)} |end|={np.linalg.norm(end):.2e} closest return to saddle={dist_back:.3e} winding about origin={wind:+.2f}")
