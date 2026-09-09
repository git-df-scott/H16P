#!/usr/bin/env python3
"""Drift-stabilised Poincare-sphere flow.  s.F = W(1-|s|^2) vanishes on S^2,
so subtracting the radial part changes nothing on the sphere but removes the
numerical drift off it."""
import numpy as np

def make(l, a, stabilise=True):
    m, b = 5.0*a, 3.0*l + 5.0
    def F(t, s):
        s1, s2, s3 = s
        Pb = -s2*s3 + l*s1*s1 + m*s1*s2 + s2*s2
        Qb = s1*s3 + a*s1*s1 + b*s1*s2
        W = s1*Pb + s2*Qb
        f = np.array([Pb - s1*W, Qb - s2*W, -s3*W])
        if stabilise:
            n2 = s1*s1 + s2*s2 + s3*s3
            f = f - s*(np.dot(s, f)/n2) - 10.0*s*(n2 - 1.0)
        return f
    return F

def jac(l, a, s, h=1e-7):
    F = make(l, a, stabilise=False); J = np.zeros((3, 3))
    for j in range(3):
        e = np.zeros(3); e[j] = h
        J[:, j] = (F(0, s+e) - F(0, s-e))/(2*h)
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
