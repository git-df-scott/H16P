#!/usr/bin/env python3
"""Is the graphic-neutrality condition N independent of the Lyapunov quantities?

eta_1 ~ V_4,  eta_2 ~ V_6 mod eta_1,  eta_3 ~ V_8 mod (eta_1,eta_2).
N(l,m,a,b) is the neutrality resultant of the two-saddle infinite graphic.
On the order-3 stratum N = 640 eta_3.  The question is what happens one step
lower, on the ORDER-2 stratum {eta_1 = 0, eta_2 != 0}."""
import sympy as S
l, m, a, b = S.symbols('l m a b', real=True)

eta1 = a*b + 2*a*l - l*m - m
eta2 = (6*a**2*b + 12*a**2*l + 3*a*b*m - a*m + b**3 - b**2*l - 4*b**2
        - 6*b*l**2 - 11*b*l - b*m**2 - 5*b - 6*l**2 - 10*l + m**2)
N = (a**2*m**3 - a*b**3 - 6*a*b**2*l - 12*a*b*l**2 + 3*a*b*l*m**2 - a*b*m**4
     - 8*a*l**3 + 6*a*l**2*m**2 - a*l*m**4 + b**4*m + 3*b**3*l*m + b**2*l*m**3
     - 4*b*l**3*m + b*l**2*m**3)

print("=== sanity: the order-3 stratum is m=5a, b=3l+5 ===")
sub3 = {m: 5*a, b: 3*l + 5}
print("eta_1 there:", S.simplify(eta1.subs(sub3)))
print("eta_2 there:", S.simplify(eta2.subs(sub3)))
print("N     there:", S.factor(N.subs(sub3)))

print("\n=== is N a multiple of eta_1? ===")
q, r = S.reduced(S.expand(N), [S.expand(eta1)], l, m, a, b)
print("N mod eta_1 =", S.factor(r))
print("  -> N is", "IN" if r == 0 else "NOT in", "the ideal (eta_1)")

print("\n=== the order-2 stratum: solve eta_1 = 0 for b ===")
bsol = S.solve(eta1, b)
print("b =", bsol)
bs = bsol[0]
N2 = S.factor(S.simplify(S.numer(S.together(N.subs(b, bs)))))
e2s = S.factor(S.simplify(S.numer(S.together(eta2.subs(b, bs)))))
print("\nN  on {eta_1=0} (numerator) =", N2)
print("\neta_2 on {eta_1=0} (numerator) =", e2s)
print("\ngcd(N, eta_2) on the stratum =", S.factor(S.gcd(S.expand(N2), S.expand(e2s))))
