#!/usr/bin/env python3
"""High-precision Chebyshev test for the 4-dimensional generating space
{T_{a-2}, T_{a-1}, T_a, U} on the period annulus of the upper centre.

If every 4x4 Wronskian-type determinant det[f_i(h_j)] over increasing
4-point subsets has the same sign, the space is an extended Chebyshev system
on the sample and M1 has at most 3 zeros there -> at most 3 cycles in this
annulus at first order -> 3+1 = 4, no five-cycle lead.
Double precision cannot decide this: the determinants are ~1e-14 relative.
"""
import sys, itertools
import numpy as np, mpmath as mp

def coeffs(a, b): return ((b-2)/(4*a), (1-b)/(a+1), b/(a+2))

def hcentre(a, b):
    A, B, C = coeffs(a, b)
    return mp.mpf(0.5)**a*(A + B*mp.mpf('0.5') + C*mp.mpf('0.25'))


def _bisect(f, lo, hi, iters=None):
    """Plain bisection to full working precision on the bracket width."""
    if iters is None:
        iters = int(3.4*mp.mp.dps) + 40
    flo = f(lo)
    for _ in range(iters):
        mid = (lo + hi)/2
        fm = f(mid)
        if fm == 0: return mid
        if (fm > 0) == (flo > 0): lo, flo = mid, fm
        else: hi = mid
    return (lo + hi)/2

def turning(a, b, h):
    A, B, C = coeffs(a, b)
    R = lambda t: h*t**(-a) - A - B*t - C*t*t
    half = mp.mpf('0.5')
    if R(half) <= 0: return None
    lo = half
    for _ in range(400):
        lo *= mp.mpf('0.7')
        if lo < mp.mpf('1e-20'): return None
        if R(lo) <= 0: break
    hi = half
    for _ in range(400):
        hi *= mp.mpf('1.4')
        if hi > mp.mpf('1e14'): return None
        if R(hi) <= 0: break
    y1 = _bisect(R, lo, half)
    y2 = _bisect(R, half, hi)
    return y1, y2

def gens(a, b, h):
    t = turning(a, b, h)
    if t is None: return None
    y1, y2 = t
    A, B, C = coeffs(a, b)
    R = lambda y: h*y**(-a) - A - B*y - C*y*y
    def S(y):
        d = (y - y1)*(y2 - y)
        v = R(y)/d if d > 0 else mp.mpf(0)
        return v if v > 0 else mp.mpf(0)
    # y = (y1+y2)/2 + ((y2-y1)/2) sin(th) removes both sqrt endpoints exactly
    c0 = (y1 + y2)/2; r0 = (y2 - y1)/2
    def I(fun, power):
        def g(th):
            y = c0 + r0*mp.sin(th)
            w = (r0*mp.cos(th))**power           # sqrt((y-y1)(y2-y)) = r0 cos(th)
            return mp.re(fun(y))*w*r0*mp.cos(th)
        return mp.quad(g, [-mp.pi/2, 0, mp.pi/2])
    out = []
    for s in (a-2, a-1, a):
        out.append(mp.re(I(lambda y: mp.sqrt(S(y))*y**s, 1)))
    out.append(mp.re(I(lambda y: S(y)**mp.mpf('1.5')*y**(a-2), 3)))
    return out

def run(a, b, n=22, dps=70):
    mp.mp.dps = dps
    a, b = mp.mpf(a), mp.mpf(b)
    hc = hcentre(a, b)
    d = None
    for dd in (mp.mpf(1), mp.mpf(-1)):
        if turning(a, b, hc + dd*mp.mpf('1e-6')) is not None: d = dd; break
    if d is None: return None
    hi = mp.mpf('1e-6')
    while hi < mp.mpf('1e8') and turning(a, b, hc + d*hi*mp.mpf('1.6')) is not None:
        hi *= mp.mpf('1.6')
    F = []
    for k in range(n):
        e = hi*mp.mpf('1e-6')**(mp.mpf(n-1-k)/(n-1))
        g = gens(a, b, hc + d*e)
        if g is None: continue
        F.append(g)
    F = mp.matrix(F).T            # 4 x N
    N = F.cols
    # normalise rows then columns by positive scalars: sign pattern is unchanged
    for r in range(4):
        m = max(abs(F[r, c]) for c in range(N))
        for c in range(N): F[r, c] /= m
    for c in range(N):
        m = max(abs(F[r, c]) for r in range(4))
        for r in range(4): F[r, c] /= m
    sgns, mags, undecided = [], [], 0
    tol = mp.mpf(10)**(-mp.mp.dps + 15)
    for js in itertools.combinations(range(N), 4):
        M = mp.matrix(4, 4)
        for r in range(4):
            for c in range(4): M[r, c] = F[r, js[c]]
        dd = mp.re(mp.det(M))
        if abs(dd) < tol:
            undecided += 1; continue
        sgns.append(mp.sign(dd)); mags.append(abs(dd))
    return N, sgns, mags, undecided

if __name__ == "__main__":
    for (a, b) in [(-0.5, 1.0), (-0.25, 0.9), (-2.5, 1.2)]:
        r = run(a, b)
        if r is None:
            print("a=%g b=%g : no annulus" % (a, b)); continue
        N, sgns, mags, und = r
        pos = sum(1 for s in sgns if s > 0); neg = sum(1 for s in sgns if s < 0)
        if not sgns:
            print("a=%-6g b=%-5g N=%d  ALL %d subsets undecided at dps=%d" % (a, b, N, und, mp.mp.dps)); continue
        print("a=%-6g b=%-5g N=%d decided=%d undecided=%d  det signs: +%d / -%d  |det| in [%s, %s]"
              % (a, b, N, len(sgns), und, pos, neg, mp.nstr(min(mags), 4), mp.nstr(max(mags), 4)))
        print("   -> %s" % ("CONSTANT SIGN: extended Chebyshev on the sample, at most 3 zeros"
                            if pos == 0 or neg == 0 else
                            "SIGN CHANGE: not Chebyshev on the sample -- 4 zeros may be possible"))
