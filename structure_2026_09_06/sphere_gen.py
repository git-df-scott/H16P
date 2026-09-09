#!/usr/bin/env python3
"""Poincare-sphere flow for the general Shi chart, drift-stabilised."""
import numpy as np

def raw(lam, l, m, a, b):
    def F(s):
        s1, s2, s3 = s
        Pb = (lam*s1 - s2)*s3 + l*s1*s1 + m*s1*s2 + s2*s2
        Qb = s1*s3 + a*s1*s1 + b*s1*s2
        W = s1*Pb + s2*Qb
        return np.array([Pb - s1*W, Qb - s2*W, -s3*W])
    return F

def make(lam, l, m, a, b, sgn=+1, k=20.0):
    R = raw(lam, l, m, a, b)
    def F(t, s):
        f = sgn*R(s); n2 = float(s @ s)
        return f - s*((s @ f)/n2) - k*s*(n2 - 1.0)
    return F

def jac(lam, l, m, a, b, s, h=1e-7):
    R = raw(lam, l, m, a, b); J = np.zeros((3, 3))
    for j in range(3):
        e = np.zeros(3); e[j] = h
        J[:, j] = (R(s+e) - R(s-e))/(2*h)
    return J

def inf_points(l, m, a, b):
    """(u, unit vector, type) for every infinite singularity, both antipodes."""
    r = np.roots([-1.0, -m, b-l, a]); out = []
    for v in r:
        if abs(v.imag) > 1e-9: continue
        u = v.real; n = np.hypot(1.0, u)
        le = -3*u*u - 2*m*u + (b-l); lt = -(u*u + m*u + l)
        kind = "saddle" if le*lt < 0 else "node"
        out.append((u, np.array([1.0/n, u/n, 0.0]), kind, le, lt))
        out.append((u, np.array([-1.0/n, -u/n, 0.0]), kind, -le, -lt))
    return out

def tangent_eigs(lam, l, m, a, b, s0):
    J = jac(lam, l, m, a, b, s0); w, v = np.linalg.eig(J); w = w.real; v = v.real
    rad = max(range(3), key=lambda i: abs(np.dot(v[:, i]/np.linalg.norm(v[:, i]), s0)))
    t = [i for i in range(3) if i != rad]
    return w[t], v[:, t]
