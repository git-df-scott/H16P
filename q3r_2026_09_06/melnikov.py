#!/usr/bin/env python3
"""First-order Melnikov function on the period annulus of the upper centre
(0,1/2) of the reversible two-centre family.

H = y^a ( x^2 + A + B y + C y^2 ),  A=(b-2)/(4a), B=(1-b)/(a+1), C=b/(a+2)
integrating factor mu = y^(a-1).   Ovals: x^2 = R(y) = h y^{-a} - A - B y - C y^2.

For a general quadratic perturbation (p,q), reversibility kills every term
except q00,q01,q02,q20 (even in x, entering via dx) and p10,p11 (odd in x,
entering via dy).  Reducing by parts,

    M1 in span{ T_{a-2}, T_{a-1}, T_a, U },
    T_s(h) = int sqrt(R) y^s dy,     U(h) = int R^{3/2} y^{a-2} dy.

So the generating space is FOUR dimensional, exactly as for Q4.
"""
import numpy as np
import mpmath as mp

def coeffs(a, b):
    return ( (b-2)/(4*a), (1-b)/(a+1), b/(a+2) )

def Rfun(a, b, h):
    A, B, C = coeffs(a, b)
    return lambda yy: h*yy**(-a) - A - B*yy - C*yy**2

def annulus_h(a, b):
    """h at the centre (0,1/2) and at the annulus boundary."""
    A, B, C = coeffs(a, b)
    hc = (0.5)**a*(A + B*0.5 + C*0.25)
    return hc

def turning(a, b, h, ylo=1e-12, yhi=1e6):
    """The two roots of R(y)=0 bracketing y=1/2."""
    R = Rfun(a, b, h)
    if R(0.5) <= 0: return None
    y1 = mp.findroot(lambda t: R(float(t)), 0.5, solver='bisect',
                     tol=1e-40) if False else None
    # bisect outward
    lo = 0.5
    while R(lo) > 0 and lo > ylo: lo *= 0.5
    if R(lo) > 0: return None
    hi = 0.5
    while R(hi) > 0 and hi < yhi: hi *= 2.0
    if R(hi) > 0: return None
    f = lambda t: R(float(t))
    y1 = mp.findroot(f, (mp.mpf(lo), mp.mpf(0.5)), solver='anderson', tol=mp.mpf('1e-30'))
    y2 = mp.findroot(f, (mp.mpf(0.5), mp.mpf(hi)), solver='anderson', tol=mp.mpf('1e-30'))
    return float(y1), float(y2)

def basis(a, b, h):
    """[T_{a-2}, T_{a-1}, T_a, U] at level h."""
    t = turning(a, b, h)
    if t is None: return None
    y1, y2 = t
    A, B, C = coeffs(a, b)
    R = lambda yy: h*yy**(-a) - A - B*yy - C*yy**2
    def integ(g):
        return mp.quad(g, [mp.mpf(y1), mp.mpf(0.5), mp.mpf(y2)])
    out = []
    for s in (a-2, a-1, a):
        out.append(float(integ(lambda yy: mp.sqrt(max(R(float(yy)), 0.0))*yy**s)))
    out.append(float(integ(lambda yy: max(R(float(yy)), 0.0)**mp.mpf('1.5')*yy**(a-2))))
    return out

if __name__ == "__main__":
    mp.mp.dps = 30
    a, b = -0.5, 1.0
    hc = annulus_h(a, b)
    print("a=%g b=%g   h(centre)=%.12f" % (a, b, hc))
    for f in (1.001, 1.01, 1.05, 1.2, 1.5, 2.0, 3.0, 6.0):
        h = hc*f
        t = turning(a, b, h)
        v = basis(a, b, h)
        print("  h/hc=%-6g  turning=%s  basis=%s"
              % (f, None if t is None else np.round(t, 6),
                 None if v is None else np.array2string(np.array(v), precision=8)))
