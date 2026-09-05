"""F16: numerical Dulac-coefficient rank at the neutral hemicycle of the two-center system
   x' = (b-2)/4 + (1-b) y + a x^2 + b y^2,   y' = -2xy,   a = -1 (neutral: infinity eigenvalues -a, -(a+2) = 1, -1).
Upper center (0, 1/2); its annulus is the upper half plane, bounded by the hemicycle (line y=0 + arc at infinity).
For each unfolding direction p_i in {da, db, e0 (y' += e0), e1 (x' += e1 x), e2 (x' += e2 xy)} we compute the
first-order log-displacement G_i(u) = dDelta u/dp_i along the ray from the center at angle pi/2 for large log-radius u,
by central finite differences, and fit it to the neutral Dulac series in w = e^{-u}:
   G(u) ~ c0 + c1 (u w) + c2 w + c3 (u w^2) + c4 w^2 + ...
The rank of the matrix (c_k, i) tells how many leading coefficients the five parameters control independently,
i.e. how many cycles can be born from the hemicycle at first order (if the next coefficient is nonzero).
"""
import sys, numpy as np, retmap as rm
b0 = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
def field(a, b, e0, e1, e2):
    return np.array([(b-2)/4, e1, 1-b, a, e2, b, e0, 0, 0, 0, -2.0, 0], float)
def center(a, b, e0, e1, e2):
    c = field(a, b, e0, e1, e2); eq = rm.equilibria(c)
    pts = [p for p in eq if p[1] > 0]
    return c, min(pts, key=lambda p: abs(p[1]-0.5))
Y = np.geomspace(1e-7, 0.3, 97)        # section: height y above the invariant line, on the ray pointing down from the center
def disp(params):
    c, pt = center(*params)
    U0 = np.log(pt[1]-Y)                     # log distance from the center along the downward ray
    u1, S, st = rm.returns_log(c[None], np.array([pt]), U0[None], th0=-np.pi/2, rtol=1e-13, umax=60, Smax=5000, maxsteps=3_000_000)
    ok = st[0] == 0
    yret = pt[1]-np.exp(u1[0])
    return np.where(ok, yret-Y, np.nan), ok
U = Y
base = np.array([-1.0, b0, 0.0, 0.0, 0.0])
D0, ok0 = disp(base); print("base (center) max |D| =", np.nanmax(np.abs(D0)), "valid down to y =", Y[ok0].min())
names = ['da', 'db', 'e0', 'e1', 'e2']; h = 1e-6; G = []
for i in range(5):
    p = base.copy(); p[i] += h; Dp, okp = disp(p)
    p = base.copy(); p[i] -= h; Dm, okm = disp(p)
    g = (Dp-Dm)/(2*h); G.append(g)
    ok = okp & okm
    print(f"{names[i]}: valid y range {Y[ok].min():.1e}..{Y[ok].max():.1e}, G at y=1e-2: {g[np.argmin(abs(Y-1e-2))]:+.4e}, y=1e-4: {g[np.argmin(abs(Y-1e-4))]:+.4e}, y=1e-6: {g[np.argmin(abs(Y-1e-6))]:+.4e}")
G = np.array(G)
np.save('data/f16_G_b%.2f.npy' % b0, np.vstack([Y, G]))
# fit on the asymptotic window
for yhi in (1e-2, 1e-3, 1e-4):
    m = (Y <= yhi) & np.all(np.isfinite(G), axis=0)
    w = Y[m]; u = -np.log(w); ulo = yhi
    A = np.column_stack([np.ones_like(u), u*w, w, u*w**2, w**2, u*w**3, w**3])
    C = np.linalg.lstsq(A, G[:, m].T, rcond=None)[0]      # 7 x 5 coefficients
    resid = np.linalg.norm(A @ C - G[:, m].T, axis=0)/np.linalg.norm(G[:, m].T, axis=0)
    print(f"\nfit window y <= {yhi:g} ({m.sum()} points): relative residuals per direction {np.round(resid, 5)}")
    lab = ['1', 'w ln', 'w', 'w2 ln', 'w^2', 'w3 ln', 'w^3']
    for k in range(5):
        print(f"  {lab[k]:6s}: " + " ".join(f"{names[i]}={C[k,i]:+.3e}" for i in range(5)))
    for kk in (2, 3, 4):
        s = np.linalg.svd(C[:kk], compute_uv=False)
        print(f"  leading {kk} coefficients x 5 directions: singular values {np.round(s/s[0], 6)} -> rank {int(np.sum(s > 1e-6*s[0]))}")
