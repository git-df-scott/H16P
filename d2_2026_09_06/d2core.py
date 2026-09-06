#!/usr/bin/env python3
"""Core for Conjecture D2 on the order-two stratum of the Shi chart.

  xdot = -y + l x^2 + m x y + y^2 ,   ydot = x (1 + a x + b y),
  lambda = 0,   eta_1 = 0  <=>  m = a(b+2l)/(l+1).

eta_2 = a (b+2l) (b-3l-5) [ a^2(b+2l+1) - (b+1)(l+1)^2 ] / (48 (l+1)^2)
      (Astra/Fable exact; independently reproduced this session up to sign)

Non-origin equilibria: (0,1), and the points of the line 1+ax+by=0 with P=0.
sigma = div X(S) = (b+2l) x_S + m y_S.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def mval(a, b, l):
    return a*(b + 2*l)/(l + 1)

def eta2(a, b, l):
    C3 = a**2*(b + 2*l + 1) - (b + 1)*(l + 1)**2
    return a*(b + 2*l)*(b - 3*l - 5)*C3/(48*(l + 1)**2)

def C3(a, b, l):
    return a**2*(b + 2*l + 1) - (b + 1)*(l + 1)**2

def field(a, b, l):
    m = mval(a, b, l)
    def F(t, z):
        x, y = z
        return [-y + l*x*x + m*x*y + y*y, x*(1 + a*x + b*y)]
    return F

def jac(a, b, l, x, y):
    m = mval(a, b, l)
    return np.array([[2*l*x + m*y, -1 + m*x + 2*y],
                     [1 + 2*a*x + b*y, b*x]])

def div(a, b, l, x, y):
    return (b + 2*l)*x + mval(a, b, l)*y

def saddles(a, b, l):
    """Non-origin equilibria that are saddles, with (x,y,sigma)."""
    m = mval(a, b, l)
    out = []
    # (0,1)
    for (x, y) in [(0.0, 1.0)]:
        J = jac(a, b, l, x, y)
        if np.linalg.det(J) < 0: out.append((x, y, div(a, b, l, x, y)))
    # line 1+ax+by=0 -> y = -(1+ax)/b  (b != 0); substitute into P=0
    if abs(b) > 1e-14:
        # P = -y + l x^2 + m x y + y^2 with y = -(1+a x)/b
        # multiply by b^2:  b(1+ax) + l b^2 x^2 - m b x (1+a x) + (1+a x)^2 = 0
        A2 = l*b*b - m*b*a + a*a
        A1 = a*b - m*b + 2*a
        A0 = b + 1.0
        if abs(A2) > 1e-14:
            disc = A1*A1 - 4*A2*A0
            if disc >= 0:
                for s in (+1, -1):
                    x = (-A1 + s*np.sqrt(disc))/(2*A2)
                    y = -(1 + a*x)/b
                    J = jac(a, b, l, x, y)
                    if np.linalg.det(J) < 0:
                        out.append((x, y, div(a, b, l, x, y)))
        elif abs(A1) > 1e-14:
            x = -A0/A1; y = -(1 + a*x)/b
            J = jac(a, b, l, x, y)
            if np.linalg.det(J) < 0: out.append((x, y, div(a, b, l, x, y)))
    return out

def splitting(a, b, l, S, T=60.0):
    """Signed separation of the saddle's two separatrices on the ray from the
    origin through the point antipodal to S.  Zero <=> homoclinic loop."""
    F = field(a, b, l)
    xs, ys = S[0], S[1]
    J = jac(a, b, l, xs, ys)
    w, v = np.linalg.eig(J); w = w.real; v = v.real
    iu, is_ = int(np.argmax(w)), int(np.argmin(w))
    vu = v[:, iu]/np.linalg.norm(v[:, iu]); vs = v[:, is_]/np.linalg.norm(v[:, is_])
    th0 = np.arctan2(ys, xs) + np.pi                      # ray opposite the saddle
    c, s = np.cos(th0), np.sin(th0)
    def ev(t, z): return z[0]*s - z[1]*c                  # crosses the line through origin
    ev.terminal = True
    res = {}
    for tag, vec, sgn in (("u", vu, +1.0), ("s", vs, -1.0)):
        best = None
        for eps in (+1e-7, -1e-7):
            z0 = np.array([xs, ys]) + eps*vec
            sol = solve_ivp(F, (0, sgn*T), z0, rtol=1e-11, atol=1e-13, events=ev,
                            dense_output=False)
            if len(sol.t_events[0]) == 0: continue
            p = sol.y_events[0][0]
            if p[0]*c + p[1]*s <= 0: continue             # must be on the correct half-ray
            r = np.hypot(p[0], p[1])
            if best is None or r < best: best = r
        res[tag] = best
    if res["u"] is None or res["s"] is None: return None
    return res["u"] - res["s"]
