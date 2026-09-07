#!/usr/bin/env python3
"""Lyapunov quantities of the origin in the Shi chart with lam=0:
   xdot = -y + l x^2 + m x y + y^2,  ydot = x + a x^2 + b x y.
Formal first integral F = (x^2+y^2)/2 + sum F_k with
   Fdot = V_4 (x^2+y^2)^2 + V_6 (x^2+y^2)^3 + ...
L(F_k) = -y F_k,x + x F_k,y is invertible on odd degrees; on even degree 2j its
cokernel is spanned by (x^2+y^2)^j, which is where the obstruction V sits."""
import sympy as S

x, y, l, m, a, b = S.symbols('x y l m a b', real=True)
P = -y + l*x**2 + m*x*y + y**2
Q = x + a*x**2 + b*x*y
P2 = l*x**2 + m*x*y + y**2
Q2 = a*x**2 + b*x*y

def homs(k):
    return [x**(k-i)*y**i for i in range(k+1)]

def L(f):
    return S.expand(-y*S.diff(f, x) + x*S.diff(f, y))

KMAX = 8
F = {2: (x**2 + y**2)/2}
V = {}
for k in range(3, KMAX+1):
    src = S.expand(-(S.diff(F[k-1], x)*P2 + S.diff(F[k-1], y)*Q2))
    basis = homs(k)
    cs = S.symbols('c0:%d' % len(basis))
    Fk = sum(c*mn for c, mn in zip(cs, basis))
    if k % 2 == 0:
        v = S.Symbol('v%d' % k)
        target = src - v*(x**2 + y**2)**(k//2)
        unk = list(cs) + [v]
    else:
        target = src
        unk = list(cs)
    eq = S.expand(L(Fk) - target)
    poly = S.Poly(eq, x, y)
    eqs = [S.expand(c) for c in poly.coeffs()]
    sol = S.solve(eqs, unk, dict=True)
    assert sol, "degree %d unsolved" % k
    sol = sol[0]
    free = [s for s in unk if s not in sol]
    sub = {f: 0 for f in free}
    F[k] = S.expand(Fk.subs(sol).subs(sub))
    if k % 2 == 0:
        V[k] = S.factor(S.simplify(S.Symbol('v%d' % k).subs(sol).subs(sub)))
        print("V_%d =" % k, V[k])

print()
eta1 = S.factor(V[4])
print("eta_1 ~ V_4 =", eta1)
print()
V6r = S.factor(S.simplify(S.reduced(S.expand(S.numer(S.together(V[6]))),
                                    [S.expand(S.numer(S.together(eta1)))], x, y, l, m, a, b)[1]))
print("V_6 reduced modulo eta_1 =", V6r)
