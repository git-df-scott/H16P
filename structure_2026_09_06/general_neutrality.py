#!/usr/bin/env python3
"""Neutrality of the two-saddle infinite graphic in the FULL Shi chart, with no
stratum imposed.

  xdot = lam x - y + l x^2 + m x y + y^2,   ydot = x + a x^2 + b x y

Infinite directions: roots of  G(u) = a + (b-l)u - m u^2 - u^3.
  lam_eq(u) = G'(u) = -3u^2 - 2m u + (b-l)
  lam_tr(u) = -(u^2 + m u + l)
Neither depends on lam: the neutrality condition lives on (l,m,a,b) alone.
"""
import sympy as S
u, e1, e2, u3, l, m, a, b, lam = S.symbols('u e1 e2 u3 l m a b lam', real=True)

A, B, C = -3, -2*m, b - l
D, E, F = -1, -m, -l
cof = S.expand((A*E - B*D)*e2 + (A*F - C*D)*e1 + (B*F - C*E))
print("cofactor  =", S.factor(cof))

# u^3 + m u^2 - (b-l) u - a = 0 with roots u1,u2,u3 ; e1=u1+u2, e2=u1u2
sol = S.solve([S.Eq(e1 + u3, -m), S.Eq(e2 + u3*e1, -(b - l))], [e1, e2], dict=True)[0]
cof_u3 = S.numer(S.together(S.simplify(cof.subs(sol))))
cubic_u3 = S.numer(S.together(S.simplify((e2*u3).subs(sol) - a)))
print("cofactor(u3) =", S.factor(cof_u3))
print("cubic(u3)    =", S.factor(cubic_u3))

N = S.factor(S.resultant(cof_u3, cubic_u3, u3))
print("\nN(l,m,a,b) = resultant =")
print(" ", N)
print("\nexpanded:", S.expand(N))

print("\n--- restricted to the order-3 stratum m=5a, b=3l+5 ---")
Nstrat = S.factor(N.subs({m: 5*a, b: 3*l + 5}))
print(" ", Nstrat)
eta3 = S.factor(-25*a*(2*a**2+l+2)*(5*a**2*l+6*a**2-3*l**3-12*l**2-15*l-6)/64)
print(" eta3 =", eta3)
print(" ratio N/eta3 on the stratum:", S.simplify(Nstrat/eta3))
