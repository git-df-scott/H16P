#!/usr/bin/env python3
"""Regression tests demanded by the review.

TEST 1 (Bautin baseline).  For the leading Bautin model D(r) = c r (r^2-r0^2)^3,
        D_rrr(r0) = 48 c r0^4,  D_rrrr(r0) = 480 c r0^3,  nu = G/(r0 H) = 1/10.
        So nu = 0.1 is the GENERIC value at a triple cycle, not a small
        "distance to a swallow-tail".

TEST 2 (coordinate rank).  With a00 = a01+a11-a10-a20-a pinning A=(1,-1) and
        I0 = a00/a11, I1 = a11^2/a20, I2 = a01^2/a11, I3 = a10/(a01 a11),
        det d(a,I0,I1,I2,I3)/d(a,a20,a11,a01,a10) = (2a-a01-a10-2a20)/(a11^2 a20^2).
        The numerator is exactly L = det J(A): the chart degenerates precisely
        where the marked point stops being an antisaddle.  The a11, a20 poles are
        coordinate-ratio poles, NOT geometric degeneration."""
import sympy as S

print("TEST 1 -- Bautin baseline")
r, r0, c = S.symbols('r r0 c')
D = c*r*(r**2 - r0**2)**3
G = S.simplify(S.diff(D, r, 3).subs(r, r0))
H = S.simplify(S.diff(D, r, 4).subs(r, r0))
print("   D_rrr(r0)  =", G, "   expected 48 c r0^4   ->", S.simplify(G - 48*c*r0**4) == 0)
print("   D_rrrr(r0) =", H, "   expected 480 c r0^3  ->", S.simplify(H - 480*c*r0**3) == 0)
nu = S.simplify(G/(r0*H))
print("   nu = G/(r0 H) =", nu, "  -> baseline nu = 1/10 = 0.1")
assert S.simplify(G - 48*c*r0**4) == 0 and S.simplify(H - 480*c*r0**3) == 0
assert nu == S.Rational(1, 10)
print("   PASS\n")

print("TEST 2 -- coordinate rank with a00 and I0 included")
a, a20, a11, a01, a10 = S.symbols('a a20 a11 a01 a10')
a00 = a01 + a11 - a10 - a20 - a
I0 = a00/a11; I1 = a11**2/a20; I2 = a01**2/a11; I3 = a10/(a01*a11)
Jm = S.Matrix([a, I0, I1, I2, I3]).jacobian([a, a20, a11, a01, a10])
det = S.simplify(S.factor(Jm.det()))
claim = (2*a - a01 - a10 - 2*a20)/(a11**2*a20**2)
print("   det =", det)
print("   claimed =", S.simplify(claim))
print("   difference =", S.simplify(det - claim))
assert S.simplify(det - claim) == 0
print("   PASS -- numerator is L = det J(A) = 2a - a01 - a10 - 2a20")
print("   => the (a, I1, I2, I3) set used earlier is NOT a complete coordinate")
print("      system; I0 was missing, and the a11/a20 poles are chart poles.")
