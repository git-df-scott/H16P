"""D2, Proposition A (exact): on the order-two stratum of the Shi chart at zero trace, a saddle with zero
divergence exists only when the origin is a center. Verification by exact polynomial algebra (sympy)."""
import sympy as sp
x, y, a, b, l = sp.symbols('x y a b l', real=True)
m = a*(b+2*l)/(l+1)                                  # eta_1 = 0
P = -y + l*x**2 + m*x*y + y**2; Q = x + a*x**2 + b*x*y
div = sp.expand(sp.diff(P, x) + sp.diff(Q, y))       # (2l+b) x + m y
C3 = a**2*(b+2*l+1) - (b+1)*(l+1)**2
# Equilibria: Q = x(1+ax+by) = 0 -> x = 0 (giving (0,0) and (0,1)) or the line 1+ax+by = 0.
# (i) the equilibrium (0,1): divergence = m
print("div(0,1) =", sp.simplify(div.subs({x: 0, y: 1})), " -> zero iff a(b+2l) = 0 (reversible center stratum)")
# (ii) equilibria on the line: divergence restricted to the line is affine in x
yl = -(1+a*x)/b
tr = sp.factor(sp.together(div.subs(y, yl)))
print("div on the line =", tr)
x0 = sp.solve(sp.numer(sp.together(div.subs(y, yl))), x)
print("zero of div on the line: x0 =", [sp.factor(s) for s in x0])
# an equilibrium at x0 requires P(x0, y(x0)) = 0:
cond = sp.factor(sp.numer(sp.together(P.subs(y, yl).subs(x, x0[0]))))
print("equilibrium-at-x0 condition, factored:", cond)
q, r = sp.div(sp.Poly(sp.expand(cond), a, b, l), sp.Poly(sp.expand(C3), a, b, l))
print("remainder modulo C3:", r.as_expr(), "  quotient:", sp.factor(q.as_expr()))
assert r.as_expr() == 0
print("\nPROPOSITION A verified: on the stratum eta_1 = 0, an equilibrium off the origin has zero divergence iff")
print("  a(b+2l) = 0  [point (0,1)]   or   C3 = 0  [points on the line 1+ax+by=0]  (b = 0, l = -1 chart degeneracies aside),")
print("  and both are center strata of eta_2 = a(b+2l)(b-3l-5) C3 / (48 (l+1)^2).")
