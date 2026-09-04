#!/usr/bin/env python3
"""Lane C: on the stratum m=5a, b=3l+5, the finite equilibria besides (0,0)
and the one on x=0 satisfy the quadratic in x from the resultant. Determine
where it has real roots (=> two more finite equilibria, generically saddle
plus antisaddle), so a finite homoclinic loop around the weak focus is
topologically possible. Then sample equilibrium types there."""
import numpy as np, sympy as S
x, y, l, a = S.symbols('x y l a', real=True)
m, b = 5*a, 3*l+5
P = -y + l*x**2 + m*x*y + y**2; Q = x + a*x**2 + b*x*y
res = S.factor(S.resultant(P, Q, y))
quad = [f for f in S.Mul.make_args(res) if S.degree(f, x) == 2][0]
c2, c1, c0 = S.Poly(quad, x).all_coeffs()
disc = S.factor(c1**2-4*c2*c0)
print("quadratic in x:", S.expand(quad)); print("discriminant:", disc)
print("leading coefficient:", S.factor(c2))
# region scan
import itertools
J = S.Matrix([[S.diff(P,x), S.diff(P,y)],[S.diff(Q,x), S.diff(Q,y)]])
print("sign of discriminant on a grid (l in [-20,5], a in (0,3]):")
for lv in np.linspace(-20, 5, 11):
    row = []
    for av in np.linspace(0.25, 3, 12):
        dv = float(disc.subs({l: lv, a: av}))
        row.append('+' if dv > 0 else '-')
    print(f"  l={lv:6.1f}: {''.join(row)}")
# sample a point with 4 finite equilibria and classify
for lv, av in ((-2.0, 1.0), (2.0, 1.0), (-10.0, 3.0), (-5.0, 0.5)):
    sols = S.solve([P.subs({l:lv,a:av}), Q.subs({l:lv,a:av})], [x,y], dict=True)
    out = []
    for s in sols:
        px, py = complex(s[x]), complex(s[y])
        if abs(px.imag) > 1e-9 or abs(py.imag) > 1e-9: continue
        Jn = np.array(J.subs({l:lv,a:av}).subs(s).evalf(), dtype=float)
        out.append((round(px.real,4), round(py.real,4), "saddle" if np.linalg.det(Jn) < 0 else "antisaddle", round(float(np.trace(Jn)),3)))
    print(f"l={lv} a={av}: eta3={float((-25*av*(2*av**2+lv+2)*(5*av**2*lv+6*av**2-3*lv**3-12*lv**2-15*lv-6))/64):.4g}  equilibria: {out}")
