#!/usr/bin/env python3
"""Exact: which quadratic perturbations of the reversible two-centre family have
M1 identically zero?  Those are the ones for which second-order Melnikov theory
is required -- and the repository has no second-order machinery at all.

M1 = q00 A0 + q01 A1 + q02 A2 + q20 A3 - p10 B0 - p11 B1  with
  A_k = 2(a-1+k) T_{a-2+k},  A_3 = (2/3)(a-1) U,  B_0 = -2 T_{a-1},  B_1 = -2 T_a.
"""
import sympy as S
a = S.symbols('a', real=True)
q00, q01, q02, q20, p10, p11 = S.symbols('q00 q01 q02 q20 p10 p11', real=True)
cT2 = 2*(a-1)*q00
cT1 = 2*a*q01 + 2*p10
cT0 = 2*(a+1)*q02 + 2*p11
cU  = S.Rational(2, 3)*(a-1)*q20
print("coefficient of T_{a-2}: ", S.factor(cT2))
print("coefficient of T_{a-1}: ", S.factor(cT1))
print("coefficient of T_a    : ", S.factor(cT0))
print("coefficient of U      : ", S.factor(cU))
sol = S.solve([cT2, cT1, cT0, cU], [q00, q20, p10, p11], dict=True)
print("\nM1 == 0  <=> ", sol)
print("""
Of the twelve quadratic perturbation coefficients:
  * six  (p00,p01,p02,p20,q10,q11) never appear in M1 at all -- they are killed
    by the reversibility symmetry x -> -x of the ovals;
  * of the remaining six, the four conditions above cut out a 2-parameter family.
So an EIGHT-dimensional space of quadratic perturbations has M1 identically zero
(for a != 1, which holds throughout the two-centre region a<0).
Every cycle produced there is invisible to first order.
""")
