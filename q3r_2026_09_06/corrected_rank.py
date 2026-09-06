#!/usr/bin/env python3
"""Corrected first-order structure for Q_3^R after Astra's relation.

Relation (Astra 2026-09-06; = eq. (11) of REVERSIBLE_RESEED_2026_09_05.md in
area-moment notation):
    (a+2) U + 3[ (b-2)/4 T_{a-2} + (1-b) T_{a-1} + b T_a ] = 0.

1. Is the remaining triple {T_{a-2}, T_{a-1}, T_a} independent (rank exactly 3)?
2. The corrected kernel of the coefficient -> M_1 map.
3. Regression test: the affine shear y -> y + eps x must have M_1 == 0.
"""
import itertools, numpy as np, mpmath as mp, sympy as S
import cheb_hp as H
mp.mp.dps = 60

print("=== 1. rank of {T_{a-2}, T_{a-1}, T_a} ===")
for (a, b) in [(-0.5, 1.0), (-0.25, 0.9), (-2.5, 1.2)]:
    a = mp.mpf(a); b = mp.mpf(b)
    hc = H.hcentre(a, b)
    d = None
    for dd in (mp.mpf(1), mp.mpf(-1)):
        if H.turning(a, b, hc + dd*mp.mpf('1e-6')) is not None: d = dd; break
    hi = mp.mpf('1e-6')
    while hi < mp.mpf('1e6') and H.turning(a, b, hc + d*hi*mp.mpf('1.6')) is not None:
        hi *= mp.mpf('1.6')
    rows = []
    for k in range(12):
        e = hi*mp.mpf('1e-5')**(mp.mpf(11-k)/11)
        g = H.gens(a, b, hc + d*e)
        if g is None: continue
        rows.append(g[:3])
    Fm = mp.matrix(rows).T                    # 3 x N
    N = Fm.cols
    for r in range(3):
        s = max(abs(Fm[r, c]) for c in range(N))
        for c in range(N): Fm[r, c] /= s
    dets = []
    for js in itertools.combinations(range(N), 3):
        Mt = mp.matrix(3, 3)
        for r in range(3):
            for c in range(3): Mt[r, c] = Fm[r, js[c]]
        dets.append(mp.re(mp.det(Mt)))
    mags = [abs(x) for x in dets]
    sg = [mp.sign(x) for x in dets if abs(x) > mp.mpf(10)**(-mp.mp.dps+12)]
    pos = sum(1 for s in sg if s > 0); neg = sum(1 for s in sg if s < 0)
    print("  a=%-6s b=%-5s N=%d  |det| in [%s, %s]  signs +%d/-%d  -> %s"
          % (mp.nstr(a,4), mp.nstr(b,4), N, mp.nstr(min(mags), 3), mp.nstr(max(mags), 3), pos, neg,
             "RANK 3, and Chebyshev on this sample" if (pos == 0) != (neg == 0)
             else "rank 3 but not Chebyshev on this sample"))

print("\n=== 2. corrected kernel conditions ===")
a, b, q00, q01, q02, q20, p10, p11 = S.symbols('a b q00 q01 q02 q20 p10 p11')
# M_1 = 2(a-1)q00 T2 + 2(a q01+p10) T1 + 2((a+1)q02+p11) T0 + (2/3)(a-1)q20 U
# eliminate U = -3/(a+2) [ (b-2)/4 T2 + (1-b) T1 + b T0 ]
U_sub = {'T2': -3*(b-2)/(4*(a+2)), 'T1': -3*(1-b)/(a+2), 'T0': -3*b/(a+2)}
cT2 = S.simplify(2*(a-1)*q00 + S.Rational(2,3)*(a-1)*q20*U_sub['T2'])
cT1 = S.simplify(2*(a*q01 + p10) + S.Rational(2,3)*(a-1)*q20*U_sub['T1'])
cT0 = S.simplify(2*((a+1)*q02 + p11) + S.Rational(2,3)*(a-1)*q20*U_sub['T0'])
print("  coeff of T_{a-2}:", S.factor(cT2))
print("  coeff of T_{a-1}:", S.factor(cT1))
print("  coeff of T_a    :", S.factor(cT0))
sol = S.solve([cT2, cT1, cT0], [q00, p10, p11], dict=True)
print("\n  M_1 == 0  <=>", sol)
print("  -> three conditions on the six visible coefficients, so the kernel is")
print("     3 (visible) + 6 (invisible) = NINE dimensional, not eight.")

print("\n=== 3. regression: the affine shear y -> y + eps x ===")
# dP = (b-1)x - 2bxy ; dQ = (a+2)x^2 + (b-2)/4 + (1-b)y + by^2
sub = {q00: (b-2)/4, q01: 1-b, q02: b, q20: a+2, p10: b-1, p11: -2*b}
print("  shear coefficients: q00=(b-2)/4, q01=1-b, q02=b, q20=a+2, p10=b-1, p11=-2b")
for nm, c in (("T_{a-2}", cT2), ("T_{a-1}", cT1), ("T_a", cT0)):
    print("    coeff of %-8s -> %s" % (nm, S.simplify(c.subs(sub))))
print("  the shear has q00 != 0 and q20 != 0, so it VIOLATES the old claimed")
print("  kernel conditions q00 = q20 = 0, yet its M_1 vanishes identically.")
