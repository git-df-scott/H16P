#!/usr/bin/env python3
"""Fable lane, step 1b: locate the annulus boundary at rho=1 by bisection on
the positive x-axis, then integrate the finite separatrices of the two
infinity saddles (in the chart x=1/z, y=v/z) to confirm the graphic itinerary:
finite orbit from the antipode of the v=rho-sqrt(kappa) saddle to the
v=rho+sqrt(kappa) saddle, closed by the equator arc that avoids the node."""
import numpy as np
from scipy.integrate import solve_ivp
rho = 1.0; den = 1+rho*rho; b = 2*(1-rho*rho)/den; c = 4*rho/den; kappa = 1+rho*rho
def f(t, u):
    x, y = u
    return [y+(6+b)*x*x+2*c*x*y-(2+b)*y*y, -x+c*x*x+(8-2*b)*x*y-c*y*y]
def returns(x0):
    ev = lambda t,u: u[1]; ev.direction = -1; ev.terminal = False
    sol = solve_ivp(f, (0, 60), [x0, 0.0], rtol=1e-11, atol=1e-13, max_step=0.02, events=ev)
    return len(sol.y_events[0]) >= 2 and np.max(np.hypot(sol.y[0], sol.y[1])) < 1e3
lo, hi = 0.2, 0.25
for _ in range(40):
    mid = (lo+hi)/2
    if returns(mid): lo = mid
    else: hi = mid
print(f"annulus boundary on +x axis: x* in ({lo:.10f},{hi:.10f})")
ev = lambda t,u: u[1]; ev.direction = -1; ev.terminal = False
sol = solve_ivp(f, (0, 60), [lo, 0.0], rtol=1e-12, atol=1e-14, max_step=0.01, events=ev)
r = np.hypot(sol.y[0], sol.y[1]); i = np.argmax(r)
print(f"  last closed orbit: max radius {r[i]:.3f} at direction angle {np.degrees(np.arctan2(sol.y[1,i],sol.y[0,i])):.2f} deg (saddle directions: {np.degrees(np.arctan(rho+np.sqrt(kappa))):.2f}, {np.degrees(np.arctan(rho-np.sqrt(kappa))):.2f}, antipodes +180)")
# separatrices of infinity saddles via the chart z=1/x, v=y/x (x>0 chart) and its antipode
P2 = lambda v: (6+b)+2*c*v-(2+b)*v*v; Q2 = lambda v: c+(8-2*b)*v-c*v*v
for v0, sign in ((rho+np.sqrt(kappa), +1), (rho-np.sqrt(kappa), -1)):
    # eigen-direction into the finite plane: perturb z; x = sign/z (sign -1 = antipodal chart, x<0)
    z0 = 1e-4
    for dv in (0.0,):
        x0 = sign/z0; y0 = (v0+dv)*x0
        # integrate forward if the finite direction is unstable in this chart, else backward
        lam_z = -P2(v0)*sign   # time reversal in the antipodal chart for degree 2
        tdir = 1.0 if lam_z > 0 else -1.0
        s = solve_ivp(f, (0, tdir*40), [x0, y0], rtol=1e-11, atol=1e-13, max_step=0.05)
        end = s.y[:, -1]; ang_end = np.degrees(np.arctan2(end[1], end[0]))
        rmin = np.min(np.hypot(s.y[0], s.y[1]))
        print(f"  saddle v0={v0:+.4f} chart sign {sign:+d}: finite separatrix is {'unstable' if tdir>0 else 'stable'}; closest approach to origin {rmin:.3f}; leaves toward |u|={np.linalg.norm(end):.2e} at angle {ang_end:.2f} deg (slope {end[1]/end[0]:+.4f})")
