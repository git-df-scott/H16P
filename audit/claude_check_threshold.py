#!/usr/bin/env python3
"""Claude independent checker for the Strike-3 threshold certificate.
1. Re-derives H at the frozen rational late-root point by direct mpmath
   quadrature of u*F*q (hypergeometric evaluation, no series, no repository
   evaluator) and compares with the JSON enclosures.
2. Re-implements the exact rational series enclosure from scratch with an
   independently derived tail bound and confirms the signs and the 5/11 value.
3. Hostile controls: the analytic path at r=5/11 exactly (H(5/11)=0 by
   construction), r slightly below/above, and the box-radius bound.
"""
import os, sys, json
from fractions import Fraction as Q
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "q4"))
import mpmath as mp
mp.mp.dps = 50
one6, five6 = mp.mpf(1)/6, mp.mpf(5)/6
F  = lambda t: mp.hyp2f1(one6, five6, 1, t)
Fp = lambda t: mp.mpf(5)/36*mp.hyp2f1(one6+1, five6+1, 2, t)
M  = lambda t: 1-6*(1-t)*Fp(t)/F(t)
def H_quad(A, B, eta, t):
    f = lambda u: u*F(u)*(A+B*u-1+(u-eta)*M(u))
    pts = [0, t/2, t] if t < 0.9 else [0, mp.mpf('0.5'), mp.mpf('0.9'), (mp.mpf('0.9')+t)/2, t]
    return mp.quad(f, pts)
cert = json.load(open(os.path.join(os.path.dirname(__file__), "..", "q4", "data", "third_threshold_certificate.json")))
A, B, eta = (Q(cert["parameters"][k]) for k in ("A", "B", "eta"))
Am, Bm, em = (mp.mpf(x.numerator)/x.denominator for x in (A, B, eta))
print("== 1. Direct quadrature versus JSON enclosures ==")
ok = True
for w in cert["witnesses"]:
    t = Q(w["t"]); tm = mp.mpf(t.numerator)/t.denominator
    val = H_quad(Am, Bm, em, tm)
    lo, hi = map(mp.mpf, w["H_enclosure"])
    inside = lo <= val <= hi
    ok &= inside and (val > 0) == (w["certified_sign"] > 0)
    print(f"t={w['t']:6s} quad={mp.nstr(val,25)} in [{mp.nstr(lo,20)},{mp.nstr(hi,20)}] -> {inside}")
tval = H_quad(Am, Bm, em, mp.mpf(5)/11)
print("H(5/11) quad =", mp.nstr(tval, 20), " JSON:", cert["threshold_evaluation"]["H_enclosure"])
ok &= tval > 0
# independent exact rational series (own derivation)
print("== 2. Independent exact rational series with own tail bound ==")
def H_exact(A, B, eta, t, N):
    # H = sum f_n[(A-1)t^{n+2}/(n+2)+B t^{n+3}/(n+3)] + sum_{n>=1} d_n [t^{n+2}/(n+2) - eta t^{n+1}/(n+1)]
    f = Q(1); tn = Q(1); S = Q(0)
    for n in range(N+1):
        S += f*tn*((A-1)*t*t/(n+2)+B*t**3/(n+3))
        if n: S += f*tn*Q(6*n, 6*n-1)*(t*t/(n+2)-eta*t/(n+1))
        f *= Q((6*n+1)*(6*n+5), 36*(n+1)**2); tn *= t
    # tail: for n>N, f_n<=f_{N+1}, d_n/f_n<=6(N+1)/(6N+5), t^n sums to t^{N+1}/(1-t)
    tail = f*tn/(1-t)*(abs(A-1)*t*t/(N+3)+abs(B)*t**3/(N+4)+Q(6*(N+1),6*(N+1)-1)*(t*t/(N+3)+abs(eta)*t/(N+2)))
    return S-tail, S+tail
for w in cert["witnesses"]:
    t = Q(w["t"]); lo, hi = H_exact(A, B, eta, t, 1200)
    sgn = 1 if lo > 0 else (-1 if hi < 0 else 0)
    print(f"t={w['t']:6s} exact sign={sgn:+d} width={float(hi-lo):.1e} margin={float(min(abs(lo),abs(hi))):.3e}")
    ok &= sgn == w["certified_sign"] and min(abs(lo), abs(hi)) > Q(1, 10**5)
lo, hi = H_exact(A, B, eta, Q(5, 11), 600); print("H(5/11) exact lower bound >0:", lo > 0); ok &= lo > 0
# box bound: |dH|<= r * int_0^t uF(u)(2+u) du <= r[t^2/(1-t)+t^3/(3(1-t))]; check the second inequality numerically at t=31/32
t = mp.mpf(31)/32
lhs = mp.quad(lambda u: u*F(u)*(2+u), [0, mp.mpf('0.9'), t]); rhs = t*t/(1-t)+t**3/(3*(1-t))
print("box bound integral", mp.nstr(lhs, 10), "<=", mp.nstr(rhs, 10), lhs <= rhs); ok &= lhs <= rhs
# check that the box radius claim uses the largest witness (31/32), i.e. 122047/3072
assert Q(31,32)**2/(1-Q(31,32))+Q(31,32)**3/(3*(1-Q(31,32))) == Q(122047, 3072)
print("== 3. Hostile controls on the analytic path at r=5/11 and both sides ==")
from q4_threshold_path import coefficients_from_r, threshold_anchors, primitive_value_closed
mp.mp.dps = 40
for r in (mp.mpf(5)/11-mp.mpf('1e-6'), mp.mpf(5)/11, mp.mpf(5)/11+mp.mpf('1e-6')):
    co = coefficients_from_r(r)
    h_at_thr = H_quad(co[0], co[1], co[2], mp.mpf(5)/11)
    # first primitive root by bisection on [0.3, 0.6]
    f = lambda x: H_quad(co[0], co[1], co[2], x)
    lo_, hi_ = mp.mpf('0.3'), mp.mpf('0.6')
    assert f(lo_) > 0 and f(hi_) < 0
    for _ in range(60):
        mid = (lo_+hi_)/2
        if f(mid) > 0: lo_ = mid
        else: hi_ = mid
    print(f"r={mp.nstr(r,12)}  H(5/11)={mp.nstr(h_at_thr,6)}  first root={mp.nstr((lo_+hi_)/2,14)}  eta={mp.nstr(co[2],10)}")
print("THRESHOLD CERTIFICATE INDEPENDENTLY", "VERIFIED" if ok else "FAILED")
