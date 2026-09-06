#!/usr/bin/env python3
"""Fully independent check of the generating-space reduction: integrate the
unperturbed reversible field for one period along an oval and accumulate
    oint mu (q dx - p dy) = int_0^T mu(y) [ q(x,y) x' - p(x,y) y' ] dt ,
then compare with the reduced expressions in T_{a-2},T_{a-1},T_a,U."""
import numpy as np, mpmath as mp
from scipy.integrate import solve_ivp
from melnikov import turning, coeffs, basis, annulus_h
mp.mp.dps = 30

def field(a, b):
    def F(t, z):
        x, y = z[0], z[1]
        return [ (b-2)/4 + (1-b)*y + a*x*x + b*y*y, -2*x*y ]
    return F

def loop_moments(a, b, h):
    """[A0,A1,A2,A3,B0,B1] by direct flow integration over one period."""
    y1, y2 = turning(a, b, h)
    F0 = field(a, b)
    z0 = [0.0, y2]                                # on the oval, x=0
    def aug(t, z):
        x, y = z[0], z[1]
        xd, yd = F0(t, z)
        mu = y**(a-1)
        return [xd, yd,
                mu*xd, mu*y*xd, mu*y*y*xd, mu*x*x*xd,      # oint mu y^k dx , mu x^2 dx
                mu*x*yd, mu*x*y*yd]                        # oint mu x dy , mu x y dy
    def ev(t, z): return z[0]
    ev.direction = 0
    sol = solve_ivp(aug, (0, 400.0), z0+[0.0]*6, rtol=1e-12, atol=1e-14, events=ev)
    ts = [t for t in sol.t_events[0] if t > 1e-8]
    if len(ts) < 2: return None, None
    T = ts[1]                                     # back to x=0 on the far side twice
    sol = solve_ivp(aug, (0, T), z0+[0.0]*6, rtol=1e-13, atol=1e-15)
    return sol.y[2:, -1], T

def reduced(a, b, h):
    T2, T1, T0, U = basis(a, b, h)
    return np.array([2*(a-1)*T2, 2*a*T1, 2*(a+1)*T0, (2.0/3.0)*(a-1)*U, -2*T1, -2*T0])

for (a, b) in [(-0.5, 1.0), (-0.3, 0.7)]:
    hc = annulus_h(a, b)
    for e in (0.01, 0.1):
        for d in (+1, -1):
            h = hc + d*e
            try:
                if turning(a, b, h) is None: continue
            except Exception:
                continue
            m, T = loop_moments(a, b, h)
            if m is None: continue
            r = reduced(a, b, h)
            rel = np.abs(m - r)/np.maximum(np.abs(r), 1e-300)
            print("a=%-5g b=%-4g h=%+.6g period=%.6f  max rel err = %.3e" % (a, b, h, T, rel.max()))
            print("   flow    :", np.array2string(m, precision=9))
            print("   reduced :", np.array2string(r, precision=9))
            break
