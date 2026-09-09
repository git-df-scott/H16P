#!/usr/bin/env python3
"""Is the graphic connection on the order-2 + neutral family realised ONLY at
eta_2 = 0?  The family meets the order-3 stratum (m=5a, b=3l+5) exactly at
    a_deg = sqrt(-(l+2)/2).
Bisect the splitting and compare a* with a_deg; then test whether the origin
there is a centre (return map = identity)."""
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
from locate import splitting
from splitting_gen import params
from sphere_gen import make

print("   l        a_deg = sqrt(-(l+2)/2)     a* (splitting zero)        a* - a_deg")
for lv in (-20.0, -14.0, -10.0, -7.0, -5.0, -3.0):
    adeg = np.sqrt(-(lv+2)/2)
    lo, hi = adeg*0.75, adeg*1.25
    try:
        flo, fhi = splitting(lv, lo), splitting(lv, hi)
        if flo is None or fhi is None or flo*fhi > 0:
            print("   %-8.2f %-25.12f  no bracket" % (lv, adeg)); continue
        astar = brentq(lambda t: splitting(lv, t), lo, hi, xtol=1e-14, rtol=8.9e-16)
    except Exception as e:
        print("   %-8.2f %-25.12f  %s" % (lv, adeg, e)); continue
    print("   %-8.2f %-25.12f  %-25.12f  %+.3e" % (lv, adeg, astar, astar - adeg), flush=True)

print("\nAt the collision point the origin should be a CENTRE (eta_1=eta_2=eta_3=0).")
print("Return map on {y=0,x>0} at l=-3, a=1/sqrt(2), m=5/sqrt(2), b=-4:")
lv = -3.0; av = 1/np.sqrt(2); mv, bv = params(lv, av)
print("   m=%.15f  5a=%.15f   b=%.15f" % (mv, 5*av, bv))
F = make(0.0, lv, mv, av, bv, sgn=+1)
for x0 in (0.02, 0.05, 0.1, 0.2, 0.3):
    s = np.array([x0, 0.0, 1.0]); s = s/np.linalg.norm(s)
    def ev(t, y): return y[1]
    ev.direction = +1
    sol = solve_ivp(F, (0, 60.0), s, rtol=1e-13, atol=1e-15, events=ev)
    got = None
    for t, p in zip(sol.t_events[0], sol.y_events[0]):
        if t > 1e-7 and p[0] > 0 and p[2] > 1e-12: got = p[0]/p[2]; break
    print("   x0=%-6g R(x0)=%.14f   d=%+.3e" % (x0, got, got - x0) if got else
          "   x0=%-6g no return" % x0)
