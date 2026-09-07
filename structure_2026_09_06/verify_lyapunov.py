#!/usr/bin/env python3
"""INDEPENDENT check of the focal values: compare the symbolic V_4, V_6, V_8
against the measured asymptotics of the actual return map.
For a focus of order k the displacement satisfies d(x) ~ C_k x^(2k+1), and the
sign of C_k must match the sign of the corresponding Lyapunov quantity."""
import numpy as np
from scipy.integrate import solve_ivp
from sphere_gen import make

def V4(l, m, a, b): return (a*b + 2*a*l - l*m - m)/8
def V6(l, m, a, b):
    return -(10*a**3*b + 20*a**3*l + 13*a**2*b*m + 16*a**2*l*m - 30*a**2*m + 23*a*b**3
             + 159*a*b**2*l + 53*a*b**2 + 350*a*b*l**2 + 232*a*b*l + 3*a*b*m**2 + 30*a*b
             + 248*a*l**3 + 252*a*l**2 - 27*a*l*m**2 + 60*a*l - 29*a*m**2 - 27*b**2*l*m
             - 27*b**2*m - 101*b*l**2*m - 138*b*l*m - 37*b*m - 124*l**3*m - 238*l**2*m
             + l*m**3 + m**3 - 124*l*m - 10*m)/192

def disp(l, m, a, b, x0, T=40.0):
    F = make(0.0, l, m, a, b, sgn=+1)
    s = np.array([x0, 0.0, 1.0]); s = s/np.linalg.norm(s)
    def ev(t, y): return y[1]
    ev.direction = +1
    sol = solve_ivp(F, (0, T), s, rtol=3e-14, atol=1e-16, events=ev)
    for t, p in zip(sol.t_events[0], sol.y_events[0]):
        if t > 1e-7 and p[0] > 0 and p[2] > 1e-12: return p[0]/p[2] - x0
    return None

def order_and_sign(l, m, a, b, xs=(0.004, 0.006, 0.009)):
    ds = [disp(l, m, a, b, x) for x in xs]
    if any(d is None for d in ds): return None
    p = np.polyfit(np.log(xs), np.log(np.abs(ds)), 1)
    return p[0], np.sign(ds[-1]), ds[-1]

print("case 1: generic point, eta_1 != 0  -> expect d ~ x^3, sign(d) = sign(V_4)")
for (l, m, a, b) in [(-10.0, 3.0, 1.0, -25.0), (-3.0, 1.0, 0.5, 2.0), (-6.0, 2.0, 0.8, -4.0)]:
    r = order_and_sign(l, m, a, b)
    v = V4(l, m, a, b)
    print("  l=%-5g m=%-4g a=%-4g b=%-6g : exponent=%.4f  sign(d)=%+d  V_4=%+.6g  sign(V_4)=%+d  %s"
          % (l, m, a, b, r[0], r[1], v, np.sign(v), "OK" if r[1] == np.sign(v) else "MISMATCH"))

print("\ncase 2: eta_1 = 0 (solve for b), eta_2 != 0 -> expect d ~ x^5, sign(d)=sign(V_6)")
for (l, m, a) in [(-10.0, 3.0, 1.0), (-3.0, 1.0, 0.5), (-6.0, 2.0, 0.8), (-3.0, 1.704347826, 0.4)]:
    b = (l*m + m - 2*a*l)/a          # eta_1 = 0
    r = order_and_sign(l, m, a, b)
    v = V6(l, m, a, b)
    print("  l=%-5g m=%-10.6g a=%-4g b=%-11.6g: exponent=%.4f  sign(d)=%+d  V_6=%+.6g  %s"
          % (l, m, a, b, r[0], r[1], v, "OK" if r[1] == np.sign(v) else "MISMATCH"))

print("\ncase 3: order-3 stratum m=5a, b=3l+5 -> expect d ~ x^7")
for (l, a) in [(-10.0, 1.0), (-8.0, 1.2), (-12.0, 0.8)]:
    m, b = 5*a, 3*l + 5
    r = order_and_sign(l, m, a, b, xs=(0.02, 0.03, 0.045))
    eta3 = -25*a*(2*a**2+l+2)*(5*a**2*l+6*a**2-3*l**3-12*l**2-15*l-6)/64
    print("  l=%-5g a=%-4g : exponent=%.4f  sign(d)=%+d  eta_3=%+.6g  %s"
          % (l, a, r[0], r[1], eta3, "OK" if r[1] == np.sign(eta3) else "MISMATCH"))
