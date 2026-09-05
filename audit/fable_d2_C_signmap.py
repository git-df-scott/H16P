"""(C): on the order-three stratum (Shi line b = 3l+5, m = 5a, lam = 0) map sign(sigma * eta_3) over all saddles
(on the line 1+ax+by = 0 and at (0,1)). Li-Cherkas + parity: a loop around the order-three focus needs sigma*eta_3 < 0."""
import numpy as np, sympy as sp
x, y, a, l = sp.symbols('x y a l', real=True)
b = 3*l+5; m = 5*a
P = -y + l*x**2 + m*x*y + y**2; Q = x + a*x**2 + b*x*y
eta3 = -25*a*(2*a**2+l+2)*(5*a**2*l+6*a**2-3*l**3-12*l**2-15*l-6)/64
J = sp.Matrix([[sp.diff(P,x), sp.diff(P,y)],[sp.diff(Q,x), sp.diff(Q,y)]])
yl = -(1+a*x)/b
Pl = sp.numer(sp.together(P.subs(y, yl)))            # quadratic in x
coeffs = [sp.lambdify((a,l), c) for c in sp.Poly(Pl, x).all_coeffs()]
trf = sp.lambdify((x,y,a,l), J.trace()); detf = sp.lambdify((x,y,a,l), J.det()); eta3f = sp.lambdify((a,l), eta3)
ylf = sp.lambdify((x,a,l), yl)
rows = []; viol = []
for av in np.linspace(-4, 4, 161):
    for lv in np.linspace(-4, 3, 141):
        if abs(lv+1) < 1e-9 or abs(3*lv+5) < 1e-9 or abs(av) < 1e-9: continue
        c2, c1, c0 = [float(c(av, lv)) for c in coeffs]
        disc = c1*c1-4*c2*c0
        e3 = eta3f(av, lv)
        cands = []
        if abs(c2) > 1e-12 and disc >= 0:
            for xs in ((-c1+np.sqrt(disc))/(2*c2), (-c1-np.sqrt(disc))/(2*c2)):
                cands.append((xs, ylf(xs, av, lv)))
        cands.append((0.0, 1.0))
        for (xs, ys) in cands:
            d = detf(xs, ys, av, lv); s = trf(xs, ys, av, lv)
            if d < 0:   # saddle
                rows.append((av, lv, xs, ys, s, e3, np.sign(s*e3)))
                if s*e3 < 0: viol.append((av, lv, xs, ys, s, e3))
rows = np.array(rows, dtype=float)
print("saddles examined:", len(rows), " sign(sigma*eta3): +", int((rows[:,6] > 0).sum()), " -", int((rows[:,6] < 0).sum()), " 0", int((rows[:,6] == 0).sum()))
print("fraction with sigma*eta3 < 0 (loop-compatible):", (rows[:,6] < 0).mean())
if viol:
    v = np.array(viol); print("loop-compatible saddles: a range", v[:,0].min(), v[:,0].max(), " l range", v[:,1].min(), v[:,1].max())
    on_line = v[np.abs(v[:,2]) > 1e-9]; at01 = v[np.abs(v[:,2]) <= 1e-9]
    print("  on the line:", len(on_line), " at (0,1):", len(at01))
np.save('audit/fable_engine/data/D2_C_signmap.npy', rows)
