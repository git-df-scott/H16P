#!/usr/bin/env python3
"""Step 1. Infinite singularities of the Shi chart on the order-3 weak-focus
stratum (lambda=0, m=5a, b=3l+5), via Poincare compactification.

Chart U1: x=1/z, y=u/z, time rescaled by z (orientation preserved for z>0):
    u' = G(1,u) + O(z),   z' = -z*P2(1,u) + O(z^2)
with G(1,u) = Q2(1,u) - u*P2(1,u).
Eigenvalues at (u*,0):  lam1 = dG/du(u*),  lam2 = -P2(1,u*).
"""
import sympy as S

x, y, u, z, l, a = S.symbols('x y u z l a', real=True)
m, b = 5*a, 3*l + 5
P = -y + l*x**2 + m*x*y + y**2
Q = x + a*x**2 + b*x*y
P2 = l*x**2 + m*x*y + y**2
Q2 = a*x**2 + b*x*y

G = S.expand(Q2.subs({x: 1, y: u}) - u*P2.subs({x: 1, y: u}))
print("G(1,u) =", S.factor(S.collect(G, u)))
print("P2(1,u) =", S.collect(S.expand(P2.subs({x: 1, y: u})), u))

# direction x=0 is an infinite singularity iff coeff of y^3 in x*Q2-y*P2 vanishes
cubic = S.expand(x*Q2 - y*P2)
print("coeff of y^3 in x*Q2-y*P2 :", S.Poly(cubic, y).coeff_monomial(y**3))

lam1 = S.expand(S.diff(G, u))
lam2 = S.expand(-P2.subs({x: 1, y: u}))
print("lam1(u) =", S.collect(lam1, u))
print("lam2(u) =", S.collect(lam2, u))

# discriminant of the direction cubic: where the number of real directions changes
disc = S.factor(S.discriminant(S.Poly(-G, u)))
print("discriminant of direction cubic (in l,a):")
print(S.simplify(disc))

# eta_3 (inherited from audit/claude_laneC_splitting4.py)
eta3 = -25*a*(2*a**2 + l + 2)*(5*a**2*l + 6*a**2 - 3*l**3 - 12*l**2 - 15*l - 6)/64
print("eta3 =", eta3)
