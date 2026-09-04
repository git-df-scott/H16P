#!/usr/bin/env python3
"""Lane C, step 1: exact Lyapunov (focus) quantities of the Shi chart
  x' = lam x - y + l x^2 + m x y + y^2,   y' = x + a x^2 + b x y,
computed symbolically at lam=0 by the standard Lyapunov-function method,
to verify the third-order weak-focus stratum m=5a, b=3l+5 recorded in
ATTACK_MATRIX.md and to obtain the exact sign of the third quantity
at the seed (l,a)=(-10,1)."""
import sympy as S
x, y, l, m, a, b = S.symbols('x y l m a b')
P = -y + l*x**2 + m*x*y + y**2
Q = x + a*x**2 + b*x*y
# Lyapunov function V = (x^2+y^2)/2 + sum V_k homogeneous; dV/dt = sum eta_k (x^2+y^2)^{k+1}
N = 8
V = (x**2+y**2)/2
etas = []
def homog(deg):
    cs = S.symbols(f'c{deg}_0:{deg+1}')
    return sum(c*x**(deg-i)*y**i for i, c in enumerate(cs)), list(cs)
for k in range(3, N+1):
    Vk, cs = homog(k)
    Vt = V+Vk
    dV = S.expand(S.diff(Vt, x)*P + S.diff(Vt, y)*Q)
    terms = S.Poly(dV, x, y)
    # keep only degree-k terms of dV; require them to equal eta*(x^2+y^2)^{k/2} if k even, else 0
    deg_terms = sum(c*x**i*y**j for (i, j), c in terms.terms() if i+j == k)
    if k % 2 == 0:
        eta = S.symbols(f'eta{k//2-1}')
        eqs = S.Poly(S.expand(deg_terms - eta*(x**2+y**2)**(k//2)), x, y).coeffs()
        sol = S.solve(eqs, cs+[eta], dict=True)[0]
        etas.append(S.factor(sol[eta]))
        V = Vt.subs(sol)
        # any free coefficient: set to zero
        V = V.subs({c: 0 for c in cs})
    else:
        eqs = S.Poly(deg_terms, x, y).coeffs()
        sol = S.solve(eqs, cs, dict=True)[0]
        V = Vt.subs(sol).subs({c: 0 for c in cs})
for i, e in enumerate(etas, 1):
    print(f"eta_{i} =", e)
e1, e2, e3 = etas[:3]
print("eta_1 = 0 solves m =", S.solve(e1, m))
sub1 = {m: S.solve(e1, m)[0]}
print("eta_2 | eta_1=0 :", S.factor(e2.subs(sub1)))
print("  b =", S.solve(S.factor(e2.subs(sub1)), b))
sub2 = {b: 3*l+5}  # the weak-focus stratum; b=-2l and the quadratic factor are center strata
e3r = S.factor(e3.subs(sub1).subs(sub2))
print("eta_3 | eta_1=eta_2=0 :", e3r)
print("eta_3 at seed (l,a)=(-10,1):", e3r.subs({l: -10, a: 1}))
print("eta_3 on the stratum vanishes where:", S.solve(S.numer(S.together(e3r)), a))
