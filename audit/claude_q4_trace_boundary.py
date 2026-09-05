#!/usr/bin/env python3
"""Fable lane, step 1d: trace the annulus-boundary component of the level curve
Hcal = 4/(9 kappa) from (0.2272,0) by following the curve (flow orthogonal to
the gradient), and report its asymptotic slopes and closest approach to the
origin, i.e. the finite orbit of the boundary graphic."""
import numpy as np, sympy as S
from scipy.integrate import solve_ivp
rho = 1.0; den = 1+rho*rho; b = 2*(1-rho*rho)/den; c = 4*rho/den; kappa = 1+rho*rho
x, y = S.symbols('x y')
Yv = c*x-(2+b)*y
Hc = (8*y*(1+Yv)-S.Rational(2,3)*(1+kappa*Yv**3))**2/(1-8*y+kappa*Yv**2)**3
level = S.expand(S.numer(S.together(Hc - S.Rational(4,9)/kappa)))
Lx = S.lambdify((x, y), S.diff(level, x)); Ly = S.lambdify((x, y), S.diff(level, y)); L = S.lambdify((x, y), level)
def curve(t, u):
    gx, gy = Lx(*u), Ly(*u); n = np.hypot(gx, gy)
    return [-gy/n, gx/n]
for direction in (+1, -1):
    sol = solve_ivp(lambda t,u: [direction*v for v in curve(t,u)], (0, 400), [0.2272112871755, 0.0], rtol=1e-10, atol=1e-12, max_step=0.01)
    X, Y = sol.y; R = np.hypot(X, Y)
    slopes = Y[R > 50]/X[R > 50]
    print(f"direction {direction:+d}: arc length reached {sol.t[-1]:.1f}, max radius {R.max():.1f}, min radius {R.min():.4f}, |level| residual max {np.max(np.abs([L(a,b_) for a,b_ in zip(X[::50],Y[::50])])):.1e}")
    print(f"   slopes at radius>50: first {slopes[0] if len(slopes) else None}, last {slopes[-1] if len(slopes) else None}")
    # crossings of the axes
    cx = [(X[i], Y[i]) for i in range(len(X)-1) if Y[i]*Y[i+1] < 0]
    print(f"   x-axis crossings: {np.round([p[0] for p in cx], 4)}")
