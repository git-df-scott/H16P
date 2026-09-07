#!/usr/bin/env python3
"""Independent check of the Melnikov reduction: compare
   M1(h) = oint mu (q dx - p dy)      (direct, parametrised by y on both branches)
with the reduced combination in the four generators.
mu = y^(a-1).  Only q00,q01,q02,q20 and p10,p11 survive; the reduction is
   A_k = oint mu y^k dx      = 2(a-1+k) T_{a-2+k}
   A_3 = oint mu x^2 dx      = (2/3)(a-1) U
   B_0 = oint mu x dy        = -2 T_{a-1}
   B_1 = oint mu x y dy      = -2 T_a
"""
import numpy as np, mpmath as mp
from melnikov import turning, coeffs, basis, annulus_h
mp.mp.dps = 30

def direct(a, b, h):
    y1, y2 = turning(a, b, h)
    A, B, C = coeffs(a, b)
    R = lambda t: h*t**(-a) - A - B*t - C*t**2
    dR = lambda t: -a*h*t**(-a-1) - B - 2*C*t
    sq = lambda t: mp.sqrt(max(R(float(t)), 0.0))
    def I(g): return mp.quad(g, [mp.mpf(y1), mp.mpf(0.5), mp.mpf(y2)])
    # oint x^i y^k dx  (i even) = - int R^{(i-1)/2} R' y^k dy
    A0 = -I(lambda t: dR(float(t))/sq(t)*t**(a-1))
    A1 = -I(lambda t: dR(float(t))/sq(t)*t**(a))
    A2 = -I(lambda t: dR(float(t))/sq(t)*t**(a+1))
    A3 = -I(lambda t: sq(t)*dR(float(t))*t**(a-1))
    # oint x^i y^k dy  (i odd) = -2 int sqrt(R) y^k dy
    B0 = -2*I(lambda t: sq(t)*t**(a-1))
    B1 = -2*I(lambda t: sq(t)*t**(a))
    return [float(v) for v in (A0, A1, A2, A3, B0, B1)]

def reduced(a, b, h):
    T2, T1, T0, U = basis(a, b, h)     # T_{a-2}, T_{a-1}, T_a, U
    return [2*(a-1)*T2, 2*a*T1, 2*(a+1)*T0, (2.0/3.0)*(a-1)*U, -2*T1, -2*T0]

for (a, b) in [(-0.5, 1.0), (-0.3, 0.7), (-2.5, 1.4)]:
    hc = annulus_h(a, b)
    for d in (1.0, -1.0):
        try:
            if turning(a, b, hc + d*1e-3) is None: continue
        except Exception:
            continue
        h = hc + d*1e-2
        try:
            dd, rr = direct(a, b, h), reduced(a, b, h)
        except Exception as e:
            print("a=%g b=%g: %s" % (a, b, e)); break
        rel = [abs(x-y)/max(abs(x), 1e-300) for x, y in zip(dd, rr)]
        print("a=%-5g b=%-4g h=%+.6g  max relative discrepancy %.3e" % (a, b, h, max(rel)))
        print("      direct  :", np.array2string(np.array(dd), precision=10))
        print("      reduced :", np.array2string(np.array(rr), precision=10))
        break
