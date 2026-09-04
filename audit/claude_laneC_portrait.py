#!/usr/bin/env python3
"""Lane C: equilibria and saddle structure of the Shi seed system
x'=-y-10x^2+5xy+y^2, y'=x+x^2-25xy, and of the Songling four-cycle point
(lambda,delta,eps small). Identifies which saddle separatrices could bound
the nest around the origin (candidate for a loop-born fourth cycle)."""
import numpy as np, sympy as S
x, y = S.symbols('x y')
def analyze(l, m, a, b, lam=0):
    P = lam*x - y + l*x**2 + m*x*y + y**2; Q = x + a*x**2 + b*x*y
    sols = S.solve([P, Q], [x, y], dict=True)
    J = S.Matrix([[S.diff(P, x), S.diff(P, y)], [S.diff(Q, x), S.diff(Q, y)]])
    for s in sols:
        pt = (complex(s[x]), complex(s[y]))
        if abs(pt[0].imag) > 1e-12 or abs(pt[1].imag) > 1e-12: continue
        Jn = np.array(J.subs(s).evalf(), dtype=float); ev = np.linalg.eigvals(Jn)
        kind = "saddle" if ev[0].real*ev[1].real < 0 else ("focus/center" if abs(ev[0].imag) > 1e-9 else "node")
        print(f"  equilibrium ({pt[0].real:+.6f},{pt[1].real:+.6f}) eig={np.round(ev,5)} {kind}  det={np.linalg.det(Jn):.4f} tr={np.trace(Jn):.4f}")
print("Shi seed (l,m,a,b)=(-10,5,1,-25):"); analyze(-10, 5, 1, -25)
print("Degree-2 part at infinity: directions where x*Q2-y*P2=0 (invariant directions):")
t = S.symbols('t')
P2 = -10 + 5*t + t**2; Q2 = 1 - 25*t   # x=1,y=t
print("  roots t of Q2 - t*P2 =", [complex(r) for r in S.Poly(S.expand(Q2 - t*P2), t).nroots()])
