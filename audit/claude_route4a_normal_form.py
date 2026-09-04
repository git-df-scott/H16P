#!/usr/bin/env python3
"""Route 4a: rotate the reversible center on the stratum into the Q3^R normal
form X' = -Y(1+kX), Y' = X + pX^2 + qY^2 and list equilibria."""
import numpy as np, sympy as S
X, Y = S.symbols('X Y')
def center_l(a):
    r = np.roots([-3, -12, 5*a*a-15, 6*a*a-6]); r = r[np.isreal(r)].real
    return min(r, key=lambda v: abs(v+1.19))
for a in (1.0, 2.0, 0.6):
    l = center_l(a); m = 5*a; b = 3*l+5
    xs = np.roots([l - m*a/b + a*a/(b*b), a/b - m/b + 2*a/(b*b), 1.0/b + 1.0/(b*b)])
    reals = [(0.0,0.0),(0.0,1.0)] + [(float(r.real), float(-(1+a*r.real)/b)) for r in xs if abs(r.imag)<1e-9]
    sad = [p for p in reals if p[0] < -1e-6][0]
    th = np.arctan2(sad[1], sad[0]); c, sn = np.cos(th), np.sin(th)
    x = c*X - sn*Y; y = sn*X + c*Y
    P = -y + l*x**2 + m*x*y + y**2; Q = x + a*x**2 + b*x*y
    Pr = S.expand(c*P + sn*Q); Qr = S.expand(-sn*P + c*Q)
    def coeffs(e):
        pol = S.Poly(e, X, Y); return {str(k): round(float(v), 8) for k, v in pol.terms() if abs(float(v)) > 1e-8}
    print(f"a={a}, l_c={l:.12f}, m={m}, b={b:.8f}, axis angle={th:.6f}, saddle {np.round(sad,6)}")
    print("   X' =", coeffs(Pr)); print("   Y' =", coeffs(Qr))
    print("   equilibria (X,Y):", [(round(c*p+sn*q,5), round(-sn*p+c*q,5)) for p, q in reals])
