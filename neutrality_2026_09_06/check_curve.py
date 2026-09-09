#!/usr/bin/env python3
"""High-precision test of the conjecture:  r(Gamma)=1  <=>  l = -2-2a^2
(the first factor of eta_3) on the two-non-antipodal-saddle region."""
import mpmath as mp
mp.mp.dps = 50

def rr(l, a):
    l, a = mp.mpf(l), mp.mpf(a)
    rts = mp.polyroots([1, 5*a, -(2*l+5), -a], maxsteps=200, extraprec=200)
    us = [r.real for r in rts if abs(mp.im(r)) < mp.mpf('1e-30')]
    sad = []
    for uu in us:
        le = -3*uu**2 - 10*a*uu + 2*l + 5
        lt = -(uu**2 + 5*a*uu + l)
        if le*lt < 0: sad.append((uu, le, lt))
    if len(sad) != 2: return None
    (u1, e1, t1), (u2, e2, t2) = sad
    return abs(e1/t1)*abs(t2/e2), u1, u2

print("ON the curve l = -2-2a^2:")
for a in ['3', '3.5', '4', '4.5', '5', '7', '10']:
    a = mp.mpf(a); l = -2 - 2*a**2
    v = rr(l, a)
    print("  a=%-5s l=%-9s r-1 = %s" % (a, l, mp.nstr(v[0]-1, 8) if v else "n/a"))
print("\nOFF the curve (l shifted by +-0.5):")
for a in ['4', '5']:
    a = mp.mpf(a)
    for d in ['-0.5', '-0.05', '0.05', '0.5']:
        l = -2 - 2*a**2 + mp.mpf(d)
        v = rr(l, a)
        print("  a=%-3s l=%-10s r-1 = %s" % (a, mp.nstr(l, 8), mp.nstr(v[0]-1, 8) if v else "n/a"))
print("\neta_3 first factor 2a^2+l+2 on the curve is identically 0, so eta_3=0 there.")
