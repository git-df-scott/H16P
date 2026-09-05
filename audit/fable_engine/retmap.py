"""Python driver for the Fable return-map engine (libretmap.so).

Coefficient vector c (12): P = c0+c1x+c2y+c3x^2+c4xy+c5y^2, Q = c6..c11.
"""
import ctypes, os, numpy as np
_lib = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), "libretmap.so"))
_lib.returns.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                         ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.POINTER(ctypes.c_double),
                         ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int),
                         ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_long]

def _p(a, t=ctypes.c_double):
    return a.ctypes.data_as(ctypes.POINTER(t))

def shift(coef, foc):
    """Translate the quadratic field so that the point foc becomes the origin (exactly an equilibrium)."""
    c = np.array(coef, float, copy=True); x0, y0 = foc
    for o in (0, 6):
        c0, c1, c2, c3, c4, c5 = c[o:o+6]
        c[o+1] = c1 + 2*c3*x0 + c4*y0
        c[o+2] = c2 + c4*x0 + 2*c5*y0
        c[o] = 0.0
    return c

def returns(coef, foc, dir_, radii, rtol=1e-10, Rmax=1e8, Tmax=1e5, maxsteps=5_000_000):
    """coef (n,12), foc (n,2), dir (n,2), radii (n,nr). Returns R, T, status arrays (n,nr).
    The field is re-expanded about each focus before integration (tolerances then scale with r)."""
    coef = np.array([shift(c, f) for c, f in zip(coef, foc)], float)
    foc = np.zeros_like(np.asarray(foc, float))
    coef = np.ascontiguousarray(coef, float); foc = np.ascontiguousarray(foc, float)
    dir_ = np.ascontiguousarray(dir_, float); radii = np.ascontiguousarray(radii, float)
    n, nr = radii.shape
    R = np.empty((n, nr)); T = np.empty((n, nr)); st = np.empty((n, nr), dtype=np.int32)
    _lib.returns(n, _p(coef), _p(foc), _p(dir_), nr, _p(radii), _p(R), _p(T), _p(st, ctypes.c_int),
                 rtol, Rmax, Tmax, maxsteps)
    return R, T, st

def shi_coef(lam, l, m, a, b, n=1.0):
    """Shi chart: x' = lam x - y + l x^2 + m x y + n y^2, y' = x + a x^2 + b x y."""
    return np.array([0, lam, -1, l, m, n, 0, 1, 0, a, b, 0], float)

def equilibria(c):
    """All real finite equilibria of the quadratic field via resultant in y."""
    c = np.asarray(c, float)
    if abs(c[5]) < 1e-14 and abs(c[11]) < 1e-14:
        # resultant in y degenerates: swap x<->y (P: c0 + c2 y + c1 x + c5 y^2 + c4 xy + c3 x^2)
        sw = c[[0, 2, 1, 5, 4, 3, 6, 8, 7, 11, 10, 9]]
        if abs(sw[5]) < 1e-14 and abs(sw[11]) < 1e-14:
            return np.zeros((0, 2))
        return equilibria(sw)[:, ::-1]
    # P as polynomial in y: A1 y^2 + B1(x) y + C1(x)
    A1 = np.poly1d([c[5]]); B1 = np.poly1d([c[4], c[2]]); C1 = np.poly1d([c[3], c[1], c[0]])
    A2 = np.poly1d([c[11]]); B2 = np.poly1d([c[10], c[8]]); C2 = np.poly1d([c[9], c[7], c[6]])
    res = (A1*C2 - A2*C1)**2 - (A1*B2 - A2*B1)*(B1*C2 - B2*C1)
    co = res.coeffs
    if np.all(np.abs(co) < 1e-14):
        return np.zeros((0, 2))
    xs = np.roots(co)
    pts = []
    for x in xs:
        if abs(x.imag) > 1e-7*(1+abs(x.real)):
            continue
        x = x.real
        # solve for y: common root of two quadratics in y; use the linear combination
        a1, b1, c1 = A1(x), B1(x), C1(x); a2, b2, c2 = A2(x), B2(x), C2(x)
        den = a1*b2 - a2*b1
        num = a2*c1 - a1*c2
        if abs(den) > 1e-12*(abs(a1*b2)+abs(a2*b1)+1e-300):
            ys = [num/den]
        else:
            ys = [r.real for r in np.roots([a1, b1, c1]) if abs(r.imag) < 1e-7] if abs(a1) > 1e-14 else ([-c1/b1] if abs(b1) > 1e-14 else [])
        for y in ys:
            P = c[0]+c[1]*x+c[2]*y+c[3]*x*x+c[4]*x*y+c[5]*y*y
            Q = c[6]+c[7]*x+c[8]*y+c[9]*x*x+c[10]*x*y+c[11]*y*y
            if abs(P)+abs(Q) < 1e-6*(1+abs(x)+abs(y))**2:
                # Newton polish
                for _ in range(5):
                    J = jac(c, x, y)
                    try:
                        d = np.linalg.solve(J, -np.array([P, Q]))
                    except np.linalg.LinAlgError:
                        break
                    x += d[0]; y += d[1]
                    P = c[0]+c[1]*x+c[2]*y+c[3]*x*x+c[4]*x*y+c[5]*y*y
                    Q = c[6]+c[7]*x+c[8]*y+c[9]*x*x+c[10]*x*y+c[11]*y*y
                pts.append((x, y))
    # dedup
    out = []
    for p in pts:
        if not any(abs(p[0]-q[0])+abs(p[1]-q[1]) < 1e-8*(1+abs(p[0])+abs(p[1])) for q in out):
            out.append(p)
    return np.array(out).reshape(-1, 2)

def jac(c, x, y):
    return np.array([[c[1]+2*c[3]*x+c[4]*y, c[2]+c[4]*x+2*c[5]*y],
                     [c[7]+2*c[9]*x+c[10]*y, c[8]+c[10]*x+2*c[11]*y]])

def antisaddles(c):
    """Equilibria with det J > 0 (foci/nodes/centers), with trace and discriminant."""
    out = []
    for (x, y) in equilibria(c):
        J = jac(c, x, y); det = np.linalg.det(J); tr = np.trace(J)
        if det > 0:
            out.append(dict(pt=(x, y), tr=tr, det=det, focus=(tr*tr < 4*det)))
    return out

def count_nest(c, focus, rmin, rmax, nr=64, dir_=(1.0, 0.0), rtol=1e-10, Rmax=1e8, Tmax=1e5, refine=True):
    """Count sign changes of D(r) = R(r) - r on a geometric grid; refine roots by bisection.
    Returns dict with roots (r, stability), grid data, near-miss info."""
    c = np.asarray(c, float)
    rad = np.geomspace(rmin, rmax, nr)
    R, T, st = returns(c[None], np.array([focus]), np.array([dir_]), rad[None], rtol, Rmax, Tmax)
    R, T, st = R[0], T[0], st[0]
    ok = st == 0
    # truncate at first failure (nest boundary)
    if not ok.all():
        k = int(np.argmin(ok)); rad, R, T, st = rad[:k], R[:k], T[:k], st[:k]
    D = R - rad
    roots = []
    for i in range(len(D)-1):
        if D[i] == 0 or D[i]*D[i+1] < 0:
            lo, hi = rad[i], rad[i+1]; Dlo, Dhi = D[i], D[i+1]
            if refine:
                for _ in range(50):
                    mid = 0.5*(lo+hi)
                    Rm, _, sm = returns(c[None], np.array([focus]), np.array([dir_]), np.array([[mid]]), rtol, Rmax, Tmax)
                    if sm[0, 0] != 0:
                        break
                    Dm = Rm[0, 0]-mid
                    if Dm*Dlo < 0: hi, Dhi = mid, Dm
                    else: lo, Dlo = mid, Dm
                    if hi-lo < 1e-12*mid: break
            root = 0.5*(lo+hi)
            stab = 'S' if Dlo > 0 else 'U'   # D>0 inside, <0 outside: orbits move toward cycle
            roots.append((root, stab))
    # near-miss: local extrema of D/r that do not cross zero
    q = D/np.maximum(rad, 1e-300)
    near = []
    for i in range(1, len(q)-1):
        if (q[i]-q[i-1])*(q[i+1]-q[i]) < 0 and np.sign(q[i-1]) == np.sign(q[i+1]) == np.sign(q[i]):
            near.append((rad[i], q[i]))
    return dict(roots=roots, rad=rad, D=D, T=T, status=st, near=near, boundary=(len(rad) < nr))

def cycles_all_nests(c, nr=64, rmin_frac=1e-4, rmax=1e6, **kw):
    """For every antisaddle, count the cycles around it. Returns list of (antisaddle, nest)."""
    res = []
    eq = equilibria(c)
    for a in antisaddles(c):
        x, y = a['pt']
        others = [np.hypot(x-p[0], y-p[1]) for p in eq if np.hypot(x-p[0], y-p[1]) > 1e-9]
        scale = min(others) if others else 1.0
        nest = count_nest(c, (x, y), rmin_frac*scale, rmax, nr=nr, **kw)
        res.append((a, nest))
    return res

if __name__ == "__main__":
    # KKL incumbent: x' = y + x^2 + x y ; y' = -10x^2 + 11/5 xy + c y^2 + alpha x + beta y
    cc, al, be = 0.7, -363889/5000, 3/2000
    coef = np.array([0, 0, 1, 1, 1, 0, 0, al, be, -10, 11/5, cc], float)
    import time; t0 = time.time()
    for a, nest in cycles_all_nests(coef, nr=80, rmax=1e7, Tmax=1e6):
        print("antisaddle", a['pt'], "tr", a['tr'], "focus", a['focus'])
        print("  roots:", [(round(r, 6), s) for r, s in nest['roots']], "boundary hit:", nest['boundary'],
              "grid", len(nest['rad']), "status last", nest['status'][-1] if len(nest['status']) else None)
        print("  near-misses:", nest['near'][:5])
    print("time", time.time()-t0)


NOISE = 5e-12   # relative displacement below which a sign change is treated as integration noise

def count_signs(rad, D, noise=NOISE):
    """Indices i where D changes sign between rad[i] and rad[i+1], both values above the noise floor."""
    idx = []
    for i in range(len(D)-1):
        if D[i]*D[i+1] < 0 and min(abs(D[i]), abs(D[i+1])) > noise*rad[i]:
            idx.append(i)
    return idx

def edge_refine(coef, foc, dir_, lo, hi, rtol, Rmax, Tmax, maxsteps, rounds=8):
    """Batched bisection between the last valid radius lo and the first failing radius hi (arrays over sets).
    Returns (r_edge, D_edge) at the innermost valid point found."""
    coef = np.asarray(coef, float); foc = np.asarray(foc, float); dir_ = np.asarray(dir_, float)
    lo = np.array(lo, float); hi = np.array(hi, float)
    r_ok = lo.copy(); D_ok = np.full(len(lo), np.nan)
    for _ in range(rounds):
        mid = np.sqrt(lo*hi)
        R, T, st = returns(coef, foc, dir_, mid[:, None], rtol, Rmax, Tmax, maxsteps)
        ok = st[:, 0] == 0
        lo = np.where(ok, mid, lo); hi = np.where(ok, hi, mid)
        r_ok = np.where(ok, mid, r_ok); D_ok = np.where(ok, R[:, 0]-mid, D_ok)
    return r_ok, D_ok


def away_dir(pt, eq):
    """Unit ray direction from pt pointing away from the nearest other equilibrium (default +x)."""
    x, y = pt; best = None
    for q in eq:
        d = np.hypot(x-q[0], y-q[1])
        if d > 1e-9 and (best is None or d < best[0]): best = (d, q)
    if best is None: return (1.0, 0.0)
    v = np.array([x-best[1][0], y-best[1][1]]); v /= np.linalg.norm(v)
    return (float(v[0]), float(v[1]))

# ---------------- compactified (log-polar) return map ----------------
_libl = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), "libretmap_log.so"))
_libl.returns_log.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.c_double,
                              ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int),
                              ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_long]

def returns_log(coef, foc, u0, th0=0.0, rtol=1e-12, umax=40.0, Smax=1e6, maxsteps=5_000_000):
    """Compactified return: coef (n,12), foc (n,2), u0 (n,nr) = log radii along the ray at angle th0 from the
    focus. Returns (u1, S1, status); displacement in log radius is u1 - u0 (D/r ~ exp(u1-u0)-1)."""
    coef = np.array([shift(c, f) for c, f in zip(coef, foc)], float)
    coef = np.ascontiguousarray(coef); u0 = np.ascontiguousarray(u0, float)
    n, nr = u0.shape
    u1 = np.empty((n, nr)); S1 = np.empty((n, nr)); st = np.empty((n, nr), dtype=np.int32)
    _libl.returns_log(n, _p(coef), nr, _p(u0), th0, _p(u1), _p(S1), _p(st, ctypes.c_int), rtol, umax, Smax, maxsteps)
    return u1, S1, st
