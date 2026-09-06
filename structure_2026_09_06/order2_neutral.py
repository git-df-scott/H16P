#!/usr/bin/env python3
"""The order-2 weak focus with a NEUTRAL boundary graphic: an explicit
two-parameter family, and the non-degeneracy checks.

On {eta_1 = 0} (which fixes b), the neutrality resultant factors as
    N = m^3 * S * W,        S = a^3 - a^2 m + 3 a l + a - l m - m
    eta_2 = m * (5a - m) * W
so the branch S = 0 gives NEUTRALITY WITHOUT eta_2 = 0.  Solving S for m:
    m = a (a^2 + 3l + 1) / (a^2 + l + 1),      b = -2l + m (l+1)/a .
"""
import sympy as S
import numpy as np

l, m, a, b = S.symbols('l m a b', real=True)
eta1 = a*b + 2*a*l - l*m - m
eta2e = (6*a**2*b + 12*a**2*l + 3*a*b*m - a*m + b**3 - b**2*l - 4*b**2
         - 6*b*l**2 - 11*b*l - b*m**2 - 5*b - 6*l**2 - 10*l + m**2)
Se = a**3 - a**2*m + 3*a*l + a - l*m - m
msol = S.solve(Se, m)[0]
bsol = S.solve(eta1, b)[0].subs(m, msol)
print("m(l,a) =", S.factor(S.simplify(msol)))
print("b(l,a) =", S.factor(S.simplify(bsol)))
sub = {m: msol, b: bsol}
print("\ncheck eta_1 on the family:", S.simplify(eta1.subs(sub)))
print("check S     on the family:", S.simplify(Se.subs(sub)))
eta2f = S.factor(S.simplify(eta2e.subs(sub)))
print("\neta_2 on the family =", eta2f)
print("  -> generically NONZERO: the focus stays of order exactly two")

# numeric non-degeneracy sweep
fm = S.lambdify((l, a), msol, 'numpy'); fb = S.lambdify((l, a), bsol, 'numpy')
fe2 = S.lambdify((l, a), eta2f, 'numpy')

def infinity(lv, mv, av, bv):
    """roots of G(u)=a+(b-l)u-m u^2-u^3; classify; return saddle ratio product."""
    r = np.roots([-1.0, -mv, bv-lv, av])
    us = sorted(v.real for v in r if abs(v.imag) < 1e-9)
    sad = []
    for u in us:
        le = -3*u*u - 2*mv*u + (bv - lv)
        lt = -(u*u + mv*u + lv)
        if le*lt < 0: sad.append((u, le, lt))
    return len(us), sad

def finite_kinds(lv, mv, av, bv):
    x, y = S.symbols('x y', real=True)
    P = -y + lv*x**2 + mv*x*y + y**2; Q = x + av*x**2 + bv*x*y
    sols = S.solve([P, Q], [x, y], dict=True)
    J = S.Matrix([[S.diff(P, x), S.diff(P, y)], [S.diff(Q, x), S.diff(Q, y)]])
    out = []
    for s in sols:
        try:
            xv, yv = complex(s[x]), complex(s[y])
        except Exception:
            continue
        if abs(xv.imag) > 1e-9 or abs(yv.imag) > 1e-9: continue
        Jn = np.array(J.subs({x: xv.real, y: yv.real}).evalf(), dtype=float)
        out.append("saddle" if np.linalg.det(Jn) < 0 else "antisaddle")
    return out

print("\n l       a      m          b          eta_2      #inf dirs  #inf saddles  r(saddle pair)   finite")
hits = []
for lv in (-12.0, -6.0, -3.0, -1.5, -0.5, 0.5, 2.0, 5.0):
    for av in (-2.0, -0.8, 0.4, 1.0, 2.5):
        try:
            mv = float(fm(lv, av)); bv = float(fb(lv, av)); e2 = float(fe2(lv, av))
        except Exception:
            continue
        if not np.isfinite(mv) or not np.isfinite(bv): continue
        nd, sad = infinity(lv, mv, av, bv)
        rr = np.nan
        if len(sad) == 2:
            (u1, e1a, t1), (u2, e2a, t2) = sad
            rr = abs(e1a/t1)*abs(t2/e2a)
        fk = finite_kinds(lv, mv, av, bv)
        ok = (len(sad) == 2 and abs(rr-1) < 1e-9 and abs(e2) > 1e-9
              and fk.count("antisaddle") >= 2)
        print(" %-7.2f %-6.2f %-10.5f %-10.5f %-10.3g %-10d %-13d %-16.12f %s%s"
              % (lv, av, mv, bv, e2, nd, len(sad), rr, fk, "   <== CANDIDATE" if ok else ""))
        if ok: hits.append((lv, av, mv, bv))
print("\npoints passing every non-degeneracy check:", len(hits))
