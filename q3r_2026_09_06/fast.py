#!/usr/bin/env python3
"""Fast, accurate evaluation of the four generators.
R(y) has simple zeros at the turning points y1<y2, so write
    sqrt(R) = sqrt((y-y1)(y2-y)) * sqrt(S(y)),   S = R/((y-y1)(y2-y)) smooth>0,
and integrate with scipy's algebraic-endpoint weight (alpha=beta=1/2)."""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

def coeffs(a, b): return ((b-2)/(4*a), (1-b)/(a+1), b/(a+2))
def hcentre(a, b):
    A, B, C = coeffs(a, b)
    return 0.5**a*(A + B*0.5 + C*0.25)

def turning(a, b, h):
    A, B, C = coeffs(a, b)
    R = lambda t: h*t**(-a) - A - B*t - C*t*t
    if not np.isfinite(R(0.5)) or R(0.5) <= 0: return None
    lo = 0.5
    for _ in range(200):
        lo *= 0.7
        if lo < 1e-14: return None
        if R(lo) <= 0: break
    else: return None
    hi = 0.5
    for _ in range(200):
        hi *= 1.4
        if hi > 1e12: return None
        if R(hi) <= 0: break
    else: return None
    return brentq(R, lo, 0.5, xtol=1e-15, rtol=1e-15), brentq(R, 0.5, hi, xtol=1e-15, rtol=1e-15)

def gens(a, b, h):
    t = turning(a, b, h)
    if t is None: return None
    y1, y2 = t
    A, B, C = coeffs(a, b)
    R = lambda y: h*y**(-a) - A - B*y - C*y*y
    def S(y):
        d = (y - y1)*(y2 - y)
        if d <= 0: 
            # endpoint: use the derivative limit
            eps = 1e-9*(y2-y1)
            yy = min(max(y, y1+eps), y2-eps)
            return R(yy)/((yy-y1)*(y2-yy))
        return R(y)/d
    out = []
    for s in (a-2, a-1, a):
        v, e = quad(lambda y: np.sqrt(max(S(y), 0.0))*y**s, y1, y2,
                    weight='alg', wvar=(0.5, 0.5), limit=200)
        out.append(v)
    v, e = quad(lambda y: max(S(y), 0.0)**1.5*y**(a-2), y1, y2,
                weight='alg', wvar=(1.5, 1.5), limit=200)
    out.append(v)
    return out

if __name__ == "__main__":
    import melnikov as M, mpmath as mp
    mp.mp.dps = 30
    for (a, b, f) in [(-0.5, 1.0, 0.01), (-0.5, 1.0, 0.5), (-0.3, 0.7, 0.05)]:
        h = hcentre(a, b) + f
        g1 = np.array(gens(a, b, h)); g2 = np.array(M.basis(a, b, h))
        print("a=%g b=%g h=%+.6g  fast=%s" % (a, b, h, np.array2string(g1, precision=10)))
        print("                  mp  =%s  max rel %.2e"
              % (np.array2string(g2, precision=10), np.abs(g1-g2).max()/np.abs(g2).max()))
