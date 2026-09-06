#!/usr/bin/env python3
"""Poincare-sphere flow with a drift-correcting term that is attracting in the
integration direction (so backward orbits are obtained by negating the field,
never by negative time)."""
import numpy as np

def raw(l, a):
    m, b = 5.0*a, 3.0*l + 5.0
    def F(s):
        s1, s2, s3 = s
        Pb = -s2*s3 + l*s1*s1 + m*s1*s2 + s2*s2
        Qb = s1*s3 + a*s1*s1 + b*s1*s2
        W = s1*Pb + s2*Qb
        return np.array([Pb - s1*W, Qb - s2*W, -s3*W])
    return F

def make(l, a, sgn=+1, k=20.0):
    R = raw(l, a)
    def F(t, s):
        f = sgn*R(s); n2 = float(s @ s)
        return f - s*((s @ f)/n2) - k*s*(n2 - 1.0)
    return F

def jac(l, a, s, h=1e-7):
    R = raw(l, a); J = np.zeros((3, 3))
    for j in range(3):
        e = np.zeros(3); e[j] = h
        J[:, j] = (R(s+e) - R(s-e))/(2*h)
    return J

def inf_points(l, a):
    r = np.roots([1.0, 5*a, -(2*l+5), -a]); pts = []
    for v in r:
        if abs(v.imag) < 1e-9:
            uu = v.real; n = np.hypot(1.0, uu)
            pts.append((uu, np.array([1.0/n, uu/n, 0.0])))
            pts.append((uu, np.array([-1.0/n, -uu/n, 0.0])))
    return pts

def tangent_eigs(l, a, s0):
    J = jac(l, a, s0); w, v = np.linalg.eig(J); w = w.real; v = v.real
    rad = max(range(3), key=lambda i: abs(np.dot(v[:, i]/np.linalg.norm(v[:, i]), s0)))
    t = [i for i in range(3) if i != rad]
    return w[t], v[:, t]
