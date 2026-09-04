#!/usr/bin/env python3
"""Fable lane, step 1: original-coordinate Q4 field at rational rho, its finite
singularities, the three infinity directions (node at v=rho, saddles at
v=rho±sqrt(1+rho^2)), the separatrices of the two infinity saddles integrated
into the finite plane, and the outer edge of the period annulus around the
center (0,0) found by integrating from a ray until orbits stop closing."""
import numpy as np, sympy as S
from scipy.integrate import solve_ivp
rho = 1.0
den = 1+rho*rho; b = 2*(1-rho*rho)/den; c = 4*rho/den; kappa = 1+rho*rho
def f(t, u):
    x, y = u
    return [y+(6+b)*x*x+2*c*x*y-(2+b)*y*y, -x+c*x*x+(8-2*b)*x*y-c*y*y]
print(f"rho={rho} kappa={kappa} b={b} c={c}")
X, Y = S.symbols('x y')
P = Y+(6+b)*X**2+2*c*X*Y-(2+b)*Y**2; Q = -X+c*X**2+(8-2*b)*X*Y-c*Y**2
J = S.Matrix([[S.diff(P,X),S.diff(P,Y)],[S.diff(Q,X),S.diff(Q,Y)]])
for s in S.solve([P,Q],[X,Y],dict=True):
    px, py = complex(s[X]), complex(s[Y])
    if abs(px.imag)>1e-9: continue
    Jn = np.array(J.subs(s).evalf(), dtype=float); ev = np.linalg.eigvals(Jn)
    print(f"  finite equilibrium ({px.real:+.6f},{py.real:+.6f}) eig={np.round(ev,4)} det={np.linalg.det(Jn):.3f}")
# infinity directions v=y/x: roots of Q2(1,v)-vP2(1,v)
P2 = lambda v: (6+b)+2*c*v-(2+b)*v*v; Q2 = lambda v: c+(8-2*b)*v-c*v*v
g = np.poly1d(np.polysub(np.poly1d([-c, 8-2*b, c]), np.polymul([1,0], [-(2+b), 2*c, 6+b])))
print("  infinity directions v:", np.round(g.roots, 6), " expected rho, rho±sqrt(1+rho^2) =", rho, rho+np.sqrt(kappa), rho-np.sqrt(kappa))
for v0 in g.roots:
    v0 = v0.real; lam_v = np.polyval(np.polyder(g), v0); lam_z = -P2(v0)
    print(f"    v0={v0:+.6f}: angular eig={lam_v:+.4f} radial eig={lam_z:+.4f} -> {'saddle' if lam_v*lam_z<0 else 'node'}; ratio |in/out|={min(abs(lam_v),abs(lam_z))/max(abs(lam_v),abs(lam_z)):.3f}")
# outer edge of the period annulus: integrate from (x0,0) and measure return
def displacement(x0, T=200):
    ev = lambda t,u: u[1]; ev.direction = -1 if x0>0 else 1; ev.terminal = False
    sol = solve_ivp(f, (0, T), [x0, 0.0], rtol=1e-11, atol=1e-13, max_step=0.02, events=ev)
    hits = [h[0] for h in sol.y_events[0] if abs(h[0]-x0) < 10*abs(x0)+1]
    if len(hits) < 2: return None, sol
    return hits[1]-x0, sol
print("  annulus scan on the positive x-axis (displacement ~0 inside the annulus):")
for x0 in (0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5):
    d, sol = displacement(x0)
    far = np.max(np.hypot(sol.y[0], sol.y[1]))
    print(f"    x0={x0}: displacement={'escape/none' if d is None else f'{d:+.2e}'}  max radius reached={far:.2f}")
print("  annulus scan on the negative x-axis:")
for x0 in (-0.05, -0.1, -0.15, -0.2, -0.25, -0.3, -0.4):
    d, sol = displacement(x0)
    far = np.max(np.hypot(sol.y[0], sol.y[1]))
    print(f"    x0={x0}: displacement={'escape/none' if d is None else f'{d:+.2e}'}  max radius reached={far:.2f}")
