#!/usr/bin/env python3
"""Exact identification of the neutrality set.

Saddle roots u1,u2 of  C(u)=u^3+5a u^2-(2l+5)u-a  carry
   lam_eq(u) = -3u^2-10a u+2l+5,   lam_tr(u) = -(u^2+5a u+l).
Neutrality of the two-saddle graphic:  lam_eq(u1) lam_tr(u2) = +/- lam_tr(u1) lam_eq(u2).
We eliminate u1,u2 and compare the result with eta_3."""
import sympy as S

u, u1, u2, l, a = S.symbols('u u1 u2 l a', real=True)
C = u**3 + 5*a*u**2 - (2*l+5)*u - a
le = lambda t: -3*t**2 - 10*a*t + 2*l + 5
lt = lambda t: -(t**2 + 5*a*t + l)

E = S.expand(le(u1)*lt(u2) - lt(u1)*le(u2))     # antisymmetric branch
F = S.expand(le(u1)*lt(u2) + lt(u1)*le(u2))     # symmetric branch
print("E factors as:", S.factor(E))

# E is divisible by (u1-u2); the graphic needs u1 != u2, so use the cofactor
Ec = S.simplify(S.cancel(E/(u1-u2)))
print("E/(u1-u2) =", S.expand(Ec))

# eliminate the two roots: express in elementary symmetric functions of u1,u2
e1, e2 = S.symbols('e1 e2', real=True)
def to_sym(expr):
    p = S.Poly(S.expand(expr), u1, u2)
    out = 0
    for (i, j), c in p.terms():
        # symmetrise
        out += c*S.symmetrize_helper if False else 0
    return None
# simpler: both roots of a common quadratic factor u^2 - e1 u + e2 of C
q = u**2 - e1*u + e2
rem = S.rem(S.Poly(C, u), S.Poly(q, u))
r0 = S.expand(rem.as_expr().coeff(u, 0)); r1 = S.expand(rem.as_expr().coeff(u, 1))
print("\nC = (u - u3)(u^2 - e1 u + e2) requires remainder 0:")
print("  coeff u^1:", S.factor(r1))
print("  coeff u^0:", S.factor(r0))

# rewrite Ec and F in e1,e2
def sym_reduce(expr):
    ex = S.expand(expr)
    ex = ex.subs({u1**2: e1*u1 - e2, u2**2: e1*u2 - e2})
    ex = S.expand(ex)
    ex = ex.subs({u1*u2: e2})
    ex = S.expand(ex)
    ex = ex.subs({u1 + u2: e1})
    # any leftover linear-in-u1/u2 asymmetry must cancel
    return S.simplify(S.expand(ex))

Ecs = sym_reduce(Ec); Fs = sym_reduce(F)
print("\nEc in (e1,e2):", S.factor(Ecs))
print("F  in (e1,e2):", S.factor(Fs))

sol = S.solve([S.Eq(r1, 0), S.Eq(r0, 0)], [e1, e2], dict=True)
print("\nnumber of (e1,e2) branches from the remainder conditions:", len(sol))
for k, s in enumerate(sol):
    ev_E = S.simplify(Ecs.subs(s)); ev_F = S.simplify(Fs.subs(s))
    print(f"  branch {k}: E -> {S.factor(ev_E)}")
    print(f"             F -> {S.factor(ev_F)}")
