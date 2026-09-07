#!/usr/bin/env python3
"""Symbolic proof that the neutrality set of the two-saddle infinite graphic on
the order-3 weak-focus stratum is exactly  2a^2 + l + 2 = 0.

lam_eq(u) = A u^2 + B u + C,  A=-3, B=-10a, C=2l+5
lam_tr(u) = D u^2 + E u + F,  D=-1, E=-5a,  F=-l
Neutrality  lam_eq(u1)lam_tr(u2) = lam_tr(u1)lam_eq(u2)  factors as
  (u1-u2) * [ (AE-BD) e2 + (AF-CD) e1 + (BF-CE) ] = 0,  e1=u1+u2, e2=u1u2.
u1,u2,u3 are the roots of  u^3 + 5a u^2 - (2l+5) u - a.
"""
import sympy as S
l, a, e1, e2, u3 = S.symbols('l a e1 e2 u3', real=True)

A, B, C = -3, -10*a, 2*l + 5
D, E, F = -1, -5*a, -l
cof = S.expand((A*E - B*D)*e2 + (A*F - C*D)*e1 + (B*F - C*E))
print("cofactor =", S.factor(cof))

# Vieta for u^3 + 5a u^2 - (2l+5) u - a with roots u1,u2,u3
rel = [S.Eq(e1 + u3, -5*a),
       S.Eq(e2 + u3*e1, -(2*l + 5)),
       S.Eq(e2*u3, a)]
sol = S.solve(rel[:2], [e1, e2], dict=True)[0]
cof_u3 = S.simplify(cof.subs(sol))
last = S.simplify(rel[2].lhs.subs(sol) - rel[2].rhs)
print("\nafter using two Vieta relations:")
print("  cofactor(u3) =", S.factor(S.numer(S.together(cof_u3))))
print("  cubic(u3)    =", S.factor(S.numer(S.together(last))))

res = S.factor(S.resultant(S.numer(S.together(cof_u3)), S.numer(S.together(last)), u3))
print("\nresultant eliminating u3:")
print(" ", res)
print("\neta_3 =", S.factor(-25*a*(2*a**2+l+2)*(5*a**2*l+6*a**2-3*l**3-12*l**2-15*l-6)/64))
