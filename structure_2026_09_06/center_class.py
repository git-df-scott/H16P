#!/usr/bin/env python3
"""Which quadratic centre class is the collision curve?
Curve:  lam=0, a = sqrt(-(l+2)/2), m = 5a, b = 3l+5   (so it lies on the
order-3 stratum with 2a^2+l+2 = 0, i.e. the first factor of eta_3)."""
import sympy as S
x, y, l = S.symbols('x y l', real=True)
a = S.sqrt(-(l+2)/2)
m, b = 5*a, 3*l+5
P = -y + l*x**2 + m*x*y + y**2
Q = x + a*x**2 + b*x*y
print("field on the curve (parameter l < -2):")
print("  xdot =", S.simplify(P)); print("  ydot =", S.simplify(Q))
div = S.simplify(S.diff(P, x) + S.diff(Q, y))
print("\ndivergence =", div, "  -> Hamiltonian?", div == 0)

# invariant straight line  f = c0 + c1 x + c2 y  with  df/dt = (k0+k1 x+k2 y) f
c0, c1, c2, k0, k1, k2 = S.symbols('c0 c1 c2 k0 k1 k2')
f = c0 + c1*x + c2*y
cond = S.expand(S.diff(f, x)*P + S.diff(f, y)*Q - (k0 + k1*x + k2*y)*f)
eqs = [S.expand(c) for c in S.Poly(cond, x, y).coeffs()]
sol = S.solve(eqs + [S.Eq(c1**2 + c2**2, 1)], [c0, c1, c2, k0, k1, k2], dict=True)
real = [s for s in sol if all(S.simplify(S.im(v)) == 0 for v in s.values() if v.free_symbols == set())]
print("\ninvariant straight lines found:", len(sol), "(any with real coefficients ->  Lotka-Volterra class)")
for s in sol[:4]:
    print("   ", {k: S.simplify(v) for k, v in s.items()})

# reversibility test: is the field invariant under (x,y)->(-x,y), t->-t  or a rotation of it?
Pr = S.simplify(-P.subs(x, -x)); Qr = S.simplify(Q.subs(x, -x))
print("\nreversible about x=0 ?", S.simplify(Pr - P) == 0 and S.simplify(Qr + Q) == 0)

for lv in (-3, -5, -10):
    av = S.sqrt(S.Rational(-(lv+2), 2))
    print("\nl=%d: a=%s  m=%s  b=%d" % (lv, av, 5*av, 3*lv+5))
    Pn = P.subs(l, lv); Qn = Q.subs(l, lv)
    sing = S.solve([Pn, Qn], [x, y], dict=True)
    print("   finite singular points:", [(S.nsimplify(s[x]), S.nsimplify(s[y])) for s in sing
                                          if S.im(S.N(s[x])) == 0 and S.im(S.N(s[y])) == 0])
