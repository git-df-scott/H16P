#!/usr/bin/env python3
"""First-order Melnikov generating space for the reversible two-centre family
  P0 = (b-2)/4 + (1-b) y + a x^2 + b y^2 ,   Q0 = -2 x y      (REVERSIBLE_RESEED (1))
Integrating factor and first integral, then the reduction of a general quadratic
perturbation's Melnikov integral to a finite generating set."""
import sympy as S
x, y, a, b, h = S.symbols('x y a b h', real=True)
P0 = (b-2)/4 + (1-b)*y + a*x**2 + b*y**2
Q0 = -2*x*y
s = S.symbols('s', real=True)
mu = y**s
print("div(mu*F) =", S.factor(S.simplify(S.diff(mu*P0, x) + S.diff(mu*Q0, y))))
sol = S.solve(S.Eq(S.simplify(S.diff(mu*P0, x) + S.diff(mu*Q0, y)), 0), s)
print("integrating factor exponent s =", sol)

A = (b-2)/(4*a); B = (1-b)/(a+1); C = b/(a+2)
H = y**a*(x**2 + A + B*y + C*y**2)
mu = y**(a-1)
print("\nH =", H)
print("check  dH/dx + mu*Q0 =", S.simplify(S.diff(H, x) + mu*Q0))
print("check  dH/dy - mu*P0 =", S.simplify(S.diff(H, y) - mu*P0))
print("\ncentres:", S.solve([P0, Q0], [x, y], dict=True))
print("H at the upper centre (0,1/2):", S.simplify(H.subs({x: 0, y: S.Rational(1,2)})))
