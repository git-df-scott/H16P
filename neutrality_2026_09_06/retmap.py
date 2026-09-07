#!/usr/bin/env python3
"""First-return map on the section {y=0, x>0}, computed on the Poincare sphere
so orbits leaving every finite window are still tracked.  Zeros of the
displacement d(x)=R(x)-x are the limit cycles."""
import numpy as np
from scipy.integrate import solve_ivp
from sphere3 import make

def R(l, a, x0, T=400.0):
    F = make(l, a, sgn=+1)
    s = np.array([x0, 0.0, 1.0]); s = s/np.linalg.norm(s)
    def ev(t, y): return y[1]
    ev.direction = +1
    def esc(t, y): return y[2] - 1e-10
    esc.terminal = True; esc.direction = -1
    sol = solve_ivp(F, (0, T), s, rtol=1e-12, atol=1e-14, events=[ev, esc])
    for p in sol.y_events[0]:
        if p[0] > 0 and p[2] > 1e-12:
            return p[0]/p[2], sol
    return None, sol

def eta3(l, a):
    return -25*a*(2*a**2+l+2)*(5*a**2*l+6*a**2-3*l**3-12*l**2-15*l-6)/64

if __name__ == "__main__":
    for (l, a) in [(-10.0, 1.0)]:
        print("l=%g a=%g eta3=%+.6g" % (l, a, eta3(l, a)))
        for x0 in [1e-3, 3e-3, 1e-2, 3e-2, 0.06, 0.1, 0.15, 0.2, 0.3, 0.5, 0.8, 1.2, 2.0, 4.0]:
            r, sol = R(l, a, x0)
            if r is None:
                print("  x0=%-8g no return  (s3_end=%.3e t=%.4f)" % (x0, sol.y[2, -1], sol.t[-1]))
            else:
                print("  x0=%-8g R=%.12f  d=%+.6e" % (x0, r, r-x0))
