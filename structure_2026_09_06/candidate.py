#!/usr/bin/env python3
"""Structure of a point on the order-2 + neutral-graphic family.
Shi chart:  xdot = lam x - y + l x^2 + m x y + y^2,  ydot = x + a x^2 + b x y."""
import numpy as np, sympy as S

def params(lv, av):
    mv = av*(av**2 + 3*lv + 1)/(av**2 + lv + 1)
    bv = -2*lv + mv*(lv + 1)/av
    return mv, bv

def report(lv, av, lam=0.0):
    mv, bv = params(lv, av)
    x, y = S.symbols('x y', real=True)
    P = lam*x - y + lv*x**2 + mv*x*y + y**2
    Q = x + av*x**2 + bv*x*y
    J = S.Matrix([[S.diff(P, x), S.diff(P, y)], [S.diff(Q, x), S.diff(Q, y)]])
    print("l=%g a=%g -> m=%.10g b=%.10g   (lam=%g)" % (lv, av, mv, bv, lam))
    print("  xdot = %s\n  ydot = %s" % (S.nsimplify(P, rational=True), S.nsimplify(Q, rational=True)))
    for s in S.solve([P, Q], [x, y], dict=True):
        try: xv, yv = complex(s[x]), complex(s[y])
        except Exception: continue
        if abs(xv.imag) > 1e-9 or abs(yv.imag) > 1e-9: continue
        Jn = np.array(J.subs({x: xv.real, y: yv.real}).evalf(), dtype=float)
        det, tr = np.linalg.det(Jn), np.trace(Jn)
        disc = tr*tr - 4*det
        kind = ("saddle" if det < 0 else
                ("focus" if disc < 0 else "node") + ("(weak)" if abs(tr) < 1e-12 else
                 ("(unstable)" if tr > 0 else "(stable)")))
        print("    finite (%+.6f, %+.6f): det=%+.4g tr=%+.4g  %s" % (xv.real, yv.real, det, tr, kind))
    r = np.roots([-1.0, -mv, bv-lv, av])
    for v in sorted(r, key=lambda z: z.real):
        if abs(v.imag) > 1e-9: continue
        u = v.real
        le = -3*u*u - 2*mv*u + (bv-lv); lt = -(u*u + mv*u + lv)
        k = "SADDLE" if le*lt < 0 else "node"
        print("    infinite u=%+.8f: lam_eq=%+.6g lam_tr=%+.6g  %s  ratio=%.10f"
              % (u, le, lt, k, abs(le/lt)))

for (lv, av) in [(-6.0, 1.0), (-12.0, 1.0), (-3.0, 0.4), (-6.0, 0.4)]:
    report(lv, av); print()
