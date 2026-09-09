#!/usr/bin/env python3
"""Focal values of the Cherkas family at a=-2 on V1=0, as exact functions of
(a11, a10, a20); then test the Newton limit against the centre conditions."""
import sympy as S
U, V, a11, a10, a20 = S.symbols('U V a11 a10 a20')
u, v = S.symbols('u v')

a = S.Integer(-2); a01 = -a11 - 3
a00 = a01 + a11 - a10 - a20 - a
Pt = S.expand((1 + (u+1)*(v-1)))
Qt = S.expand(a00 + a10*(u+1) + a20*(u+1)**2 + a01*(v-1) + a11*(u+1)*(v-1) + a*(v-1)**2)
A11 = S.diff(Pt,u).subs({u:0,v:0}); A12 = S.diff(Pt,v).subs({u:0,v:0})
w2 = S.simplify(-a10 + a11 - 2*a20 - 1); w = S.sqrt(w2)
sub = {u: U, v: (-A11*U + w*V)/A12}
Ud = S.expand(Pt.subs(sub, simultaneous=True))
Vd = S.expand(S.simplify((A12*Qt.subs(sub, simultaneous=True) + A11*Ud)/w))
# rescale time by w so linear part is exactly (V, -U)
Ud = S.expand(S.cancel(Ud/w)); Vd = S.expand(S.cancel(Vd/w))
Ud2 = S.expand(Ud - V); Vd2 = S.expand(Vd + U)      # nonlinear parts
def homs(k): return [U**(k-i)*V**i for i in range(k+1)]
def L(f): return S.expand(V*S.diff(f,U) - U*S.diff(f,V))
F = {2: (U**2+V**2)/2}; Vs = {}
for k in range(3, 7):
    src = S.expand(-(S.diff(F[k-1],U)*Ud2 + S.diff(F[k-1],V)*Vd2))
    src = S.expand(src)
    basis = homs(k); cs = S.symbols('c0:%d' % len(basis))
    Fk = sum(c*m for c, m in zip(cs, basis))
    if k % 2 == 0:
        vv = S.Symbol('v%d' % k); target = src - vv*(U**2+V**2)**(k//2); unk = list(cs)+[vv]
    else:
        target = src; unk = list(cs)
    eq = S.expand(L(Fk) - target)
    sol = S.solve([S.expand(c) for c in S.Poly(eq, U, V).coeffs()], unk, dict=True)[0]
    free = [s for s in unk if s not in sol]
    F[k] = S.expand(Fk.subs(sol).subs({f:0 for f in free}))
    if k % 2 == 0:
        Vs[k] = S.simplify(S.Symbol('v%d'%k).subs(sol).subs({f:0 for f in free}))
V3 = S.simplify(S.factor(S.numer(S.together(Vs[4]))))
print("V3 numerator (a=-2, V1=0):"); print("  ", V3)
print()
sols = S.solve(V3, a10)
print("V3 = 0  =>  a10 =", [S.simplify(s) for s in sols])
import mpmath as mp
mp.mp.dps = 40
vals = {a11: S.Rational(0), a20: S.Rational(-3254, 675)}
a11n = mp.mpf("8.067217612039481856440047"); a10n = mp.mpf("16.53363941409600843877444")
for k, s in enumerate(sols):
    f = S.lambdify((a11, a20), s, 'mpmath')
    got = f(a11n, mp.mpf(-3254)/675)
    print("   branch %d at a11=%s : a10 = %s   (Newton gave %s)  diff %s"
          % (k, mp.nstr(a11n, 20), mp.nstr(got, 22), mp.nstr(a10n, 22), mp.nstr(got - a10n, 6)))
