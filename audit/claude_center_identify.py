#!/usr/bin/env python3
"""Identify the center at the loop point on the stratum (a=1, l=-1.1835...):
reversibility test across the line through the origin and the saddle, and
finite singularities. Also: finite singular points of the Q4 center family."""
import numpy as np, sympy as S
x, y = S.symbols('x y')
a = 1.0; l = -1.183503419072; m = 5*a; b = 3*l+5
P = -y+l*x**2+m*x*y+y**2; Q = x+a*x**2+b*x*y
sols = [(complex(s[x]), complex(s[y])) for s in S.solve([P,Q],[x,y],dict=True)]
print("center-curve point a=1: equilibria", [(round(p.real,5), round(q.real,5)) for p,q in sols if abs(p.imag)<1e-9 and abs(q.imag)<1e-9])
# reversibility: rotate so the saddle lies on the x-axis; test P(x,-y)=-P(x,y)?? For reversible systems w.r.t. y->-y with time reversal: dx/dt odd in y, dy/dt even in y.
sad = [(p.real,q.real) for p,q in sols if abs(p.imag)<1e-9 and abs(q.imag)<1e-9 and abs(p)>1e-6][0]
th = np.arctan2(sad[1], sad[0])
c, s_ = np.cos(th), np.sin(th)
X, Y = S.symbols('X Y')
xr = c*X - s_*Y; yr = s_*X + c*Y   # (x,y) = R(theta)(X,Y)
Pr = S.expand(c*P.subs({x:xr,y:yr}) + s_*Q.subs({x:xr,y:yr}))   # dX/dt
Qr = S.expand(-s_*P.subs({x:xr,y:yr}) + c*Q.subs({x:xr,y:yr}))  # dY/dt
def parts(expr):
    pol = S.Poly(expr, X, Y); even = odd = 0
    for (i,j), cf in pol.terms():
        if j % 2 == 0: even += abs(float(cf))
        else: odd += abs(float(cf))
    return even, odd
print("rotated field: dX/dt (even_in_Y, odd_in_Y) =", parts(Pr), "  dY/dt (even,odd) =", parts(Qr))
print("reversible w.r.t. Y->-Y, t->-t requires dX/dt odd in Y and dY/dt even in Y")
# Q4 family other singular points
print("--- Q4 centers: finite singular points for rho in (0.5,1,2) ---")
for rho in (0.5, 1.0, 2.0):
    den = 1+rho*rho; bq = 2*(1-rho*rho)/den; cq = 4*rho/den
    Pq = y+(6+bq)*x**2+2*cq*x*y-(2+bq)*y**2; Qq = -x+cq*x**2+(8-2*bq)*x*y-cq*y**2
    J = S.Matrix([[S.diff(Pq,x),S.diff(Pq,y)],[S.diff(Qq,x),S.diff(Qq,y)]])
    out = []
    for s in S.solve([Pq,Qq],[x,y],dict=True):
        px, py = complex(s[x]), complex(s[y])
        if abs(px.imag)>1e-9 or abs(py.imag)>1e-9: continue
        Jn = np.array(J.subs(s).evalf(), dtype=float); ev = np.linalg.eigvals(Jn)
        kind = "saddle" if np.linalg.det(Jn)<0 else ("focus" if abs(ev[0].imag)>1e-9 else "node")
        out.append((round(px.real,4), round(py.real,4), kind, round(float(np.trace(Jn)),3)))
    print(f"rho={rho} kappa={1+rho*rho}: {out}")
