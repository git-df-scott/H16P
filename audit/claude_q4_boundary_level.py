#!/usr/bin/env python3
"""Fable lane, step 1c: the annulus boundary as the level set of the rational
first integral Hcal = [8y(1+Yv) - (2/3)(1+kappa Yv^3)]^2 / [1-8y+kappa Yv^2]^3,
Yv = c x - (2+b) y, at the loop value 4/(9 kappa). Check: Hcal(0,0)=4/9 (center),
trace the level curve around the origin, and find its asymptotic directions."""
import numpy as np, sympy as S
rho = 1.0; den = 1+rho*rho; b = 2*(1-rho*rho)/den; c = 4*rho/den; kappa = 1+rho*rho
x, y = S.symbols('x y')
Yv = c*x-(2+b)*y
Hc = (8*y*(1+Yv)-S.Rational(2,3)*(1+kappa*Yv**3))**2/(1-8*y+kappa*Yv**2)**3
print("Hcal(0,0) =", Hc.subs({x:0,y:0}), " loop value 4/(9 kappa) =", S.Rational(4,9)/kappa)
# invariance check: dH/dt = Hx*P + Hy*Q should vanish identically
P = y+(6+b)*x**2+2*c*x*y-(2+b)*y**2; Q = -x+c*x**2+(8-2*b)*x*y-c*y**2
dH = S.simplify(S.diff(Hc,x)*P+S.diff(Hc,y)*Q)
print("dHcal/dt simplified:", dH)
# level curve at loop value: numerator - value*denominator = 0 (polynomial)
level = S.expand(S.numer(S.together(Hc - S.Rational(4,9)/kappa)))
print("degree of level polynomial:", S.Poly(level, x, y).total_degree())
# asymptotic directions: highest-degree homogeneous part
poly = S.Poly(level, x, y); deg = poly.total_degree()
top = sum(cf*x**i*y**j for (i,j), cf in poly.terms() if i+j == deg)
v = S.symbols('v')
dirs = S.Poly(S.expand(top.subs({x:1, y:v})), v)
print("top-degree part in v=y/x:", S.factor(dirs.as_expr()))
print("real asymptotic slopes:", [complex(r) for r in dirs.nroots() if abs(complex(r).imag) < 1e-9], " saddle slopes rho±sqrt(kappa):", rho+np.sqrt(kappa), rho-np.sqrt(kappa))
# numerical trace: intersections of the level curve with the x-axis and the y-axis
f_x = S.lambdify(x, level.subs(y, 0)); f_y = S.lambdify(y, level.subs(x, 0))
print("level curve on x-axis (roots):", [complex(r) for r in S.Poly(level.subs(y,0), x).nroots() if abs(complex(r).imag)<1e-9])
print("level curve on y-axis (roots):", [complex(r) for r in S.Poly(level.subs(x,0), y).nroots() if abs(complex(r).imag)<1e-9])
