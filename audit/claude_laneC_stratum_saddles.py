#!/usr/bin/env python3
"""Lane C: over the Attack-2 box on the third-order weak-focus stratum
(m=5a, b=3l+5, lambda=0), count finite equilibria and finite saddles.
A finite saddle is needed for a finite homoclinic loop around the origin
nest; without one, a fourth cycle in that nest could only come from a
separatrix cycle through infinity."""
import numpy as np, sympy as S
x, y, l, a = S.symbols('x y l a', real=True)
m, b = 5*a, 3*l+5
P = -y + l*x**2 + m*x*y + y**2; Q = x + a*x**2 + b*x*y
# eliminate: Q=0 -> x=0 or 1 + a x + b y = 0
res = S.resultant(P, Q, y)
print("resultant in x (factored):", S.factor(res))
J = S.Matrix([[S.diff(P,x), S.diff(P,y)],[S.diff(Q,x), S.diff(Q,y)]])
for lv in (-12, -10, -8):
    for av in (0.8, 1.0, 1.2):
        sols = S.solve([P.subs({l:lv,a:av}), Q.subs({l:lv,a:av})], [x,y], dict=True)
        reals = [s for s in sols if abs(complex(s[x]).imag)<1e-10 and abs(complex(s[y]).imag)<1e-10]
        kinds = []
        for s in reals:
            Jn = np.array(J.subs({l:lv,a:av}).subs(s).evalf(), dtype=float)
            kinds.append("saddle" if np.linalg.det(Jn) < 0 else "antisaddle")
        print(f"l={lv} a={av}: {len(reals)} real equilibria: {kinds}")
