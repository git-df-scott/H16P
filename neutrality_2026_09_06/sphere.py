#!/usr/bin/env python3
"""Poincare-sphere compactification of the Shi chart on the order-3 stratum.

s=(s1,s2,s3) on S^2, x=s1/s3, y=s2/s3.  With
    Pbar = P0 s3^2 + P1(s1,s2) s3 + P2(s1,s2),   Qbar likewise,
the field (time-rescaled by s3, orientation preserved on s3>0) is
    s1' = Pbar - s1 W,  s2' = Qbar - s2 W,  s3' = -s3 W,   W = s1 Pbar + s2 Qbar.
It is polynomial and smooth across the equator s3=0.
"""
import numpy as np

def make(l, a):
    m, b = 5.0*a, 3.0*l + 5.0
    def F(t, s):
        s1, s2, s3 = s
        Pb = -s2*s3 + l*s1*s1 + m*s1*s2 + s2*s2
        Qb = s1*s3 + a*s1*s1 + b*s1*s2
        W = s1*Pb + s2*Qb
        return np.array([Pb - s1*W, Qb - s2*W, -s3*W])
    return F

def jac(l, a, s):
    h = 1e-7; F = make(l, a); J = np.zeros((3, 3))
    for j in range(3):
        e = np.zeros(3); e[j] = h
        J[:, j] = (F(0, s+e) - F(0, s-e))/(2*h)
    return J

def inf_points(l, a):
    """Infinite singularities as unit vectors on the equator (both antipodes)."""
    r = np.roots([1.0, 5*a, -(2*l+5), -a])
    pts = []
    for v in r:
        if abs(v.imag) < 1e-9:
            uu = v.real; n = np.hypot(1.0, uu)
            pts.append((uu, np.array([1.0/n, uu/n, 0.0])))
            pts.append((uu, np.array([-1.0/n, -uu/n, 0.0])))
    return pts

def finite_points(l, a):
    """Finite equilibria lifted to the sphere."""
    import sympy as S
    x, y = S.symbols('x y', real=True)
    m, b = 5*a, 3*l+5
    P = -y + l*x**2 + m*x*y + y**2; Q = x + a*x**2 + b*x*y
    sols = S.solve([P, Q], [x, y], dict=True)
    out = []
    for s in sols:
        xv, yv = complex(s[x]), complex(s[y])
        if abs(xv.imag) < 1e-9 and abs(yv.imag) < 1e-9:
            xv, yv = xv.real, yv.real
            n = np.sqrt(1+xv*xv+yv*yv)
            out.append(((xv, yv), np.array([xv/n, yv/n, 1.0/n])))
    return out

if __name__ == "__main__":
    for (l, a) in [(-10.0, 1.0), (-12.0, 0.8), (-8.0, 1.2)]:
        print(f"=== l={l} a={a} ===")
        for xy, s in finite_points(l, a):
            J = jac(l, a, s); w = np.linalg.eigvals(J)
            w = sorted(w, key=lambda t: -abs(t.real))[:2]
            print(f"  finite {xy}: sphere eigs {np.round(w,6)}")
        for uu, s in inf_points(l, a):
            J = jac(l, a, s); w = np.linalg.eigvals(J)
            w = np.array(sorted(w, key=lambda t: t.real))
            print(f"  infinite u={uu:+.6f} s={np.round(s,4)}: eigs {np.round(w.real,6)}")
