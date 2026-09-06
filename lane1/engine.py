"""Lane 1 engine driver: Andronov-Hopf curve beta*(s) of a uniformly rotated
quadratic family.

Coefficient convention (vec12), matching SEEDS.json:
    P = c0 + c1 x + c2 y + c3 x^2 + c4 x y + c5 y^2
    Q = c6 + c7 x + c8 y + c9 x^2 + c10 x y + c11 y^2
    xdot = P, ydot = Q.

Rotated family  X_b = (P cos b - Q sin b, P sin b + Q cos b).
"""
import ctypes, hashlib, os, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "retmap1.c")
LIB = os.path.join(HERE, "libretmap1.so")


def _build():
    src_hash = hashlib.sha256(open(SRC, "rb").read()).hexdigest()
    stamp = LIB + ".sha256"
    if os.path.exists(LIB) and os.path.exists(stamp):
        if open(stamp).read().strip() == src_hash:
            return src_hash
    cmd = ["cc", "-O3", "-march=native", "-fno-fast-math", "-fopenmp",
           "-shared", "-fPIC", SRC, "-o", LIB, "-lm"]
    subprocess.run(cmd, check=True)
    open(stamp, "w").write(src_hash)
    return src_hash


ENGINE_HASH = _build()
ENGINE_NAME = "lane1/retmap1.c"

_lib = ctypes.CDLL(LIB)
_dbl = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
_int = np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS")

_lib.d_curve.restype = None
_lib.d_curve.argtypes = [_dbl, ctypes.c_double, _dbl, ctypes.c_int,
                         ctypes.c_double, ctypes.c_double, ctypes.c_double,
                         ctypes.c_double, ctypes.c_long,
                         _dbl, _int, _dbl]
_lib.betastar.restype = None
_lib.betastar.argtypes = [_dbl, ctypes.c_double, _dbl, ctypes.c_int,
                          ctypes.c_double, ctypes.c_double, ctypes.c_double,
                          ctypes.c_long, ctypes.c_double, ctypes.c_double,
                          ctypes.c_int, _dbl, _int, _dbl, _dbl]

DEFAULTS = dict(rtol=1e-12, Tmax=400.0, Rmax=1.0e4, nstep=4_000_000)


# ----------------------------------------------------------------- fields
def local10(vec12, focus):
    """Re-expand the field about `focus`; returns the local 10-vector.
    (bug A2: never integrate in a chart whose origin is far from the orbit)"""
    c = np.asarray(vec12, dtype=float)
    fx, fy = float(focus[0]), float(focus[1])
    out = np.empty(10)
    for k, off in ((0, 0), (5, 6)):
        c0, c1, c2, c3, c4, c5 = c[off:off + 6]
        out[k + 0] = c1 + 2 * c3 * fx + c4 * fy          # u
        out[k + 1] = c2 + c4 * fx + 2 * c5 * fy          # v
        out[k + 2] = c3                                  # u^2
        out[k + 3] = c4                                  # uv
        out[k + 4] = c5                                  # v^2
    return np.ascontiguousarray(out)


def residual_at(vec12, pt):
    c = np.asarray(vec12, float)
    x, y = pt
    m = np.array([1, x, y, x * x, x * y, y * y])
    return float(c[0:6] @ m), float(c[6:12] @ m)


def equilibria(vec12, tol=1e-9):
    """All real equilibria of the quadratic field, by resultant + polish."""
    import numpy.polynomial.polynomial as npoly
    c = np.asarray(vec12, float)
    # P and Q as polynomials in y with x-dependent coefficients
    # P = (c5) y^2 + (c2 + c4 x) y + (c0 + c1 x + c3 x^2)
    # Q = (c11) y^2 + (c8 + c10 x) y + (c6 + c7 x + c9 x^2)
    def coeffs(o):
        return (np.array([c[o + 5]]),
                np.array([c[o + 2], c[o + 4]]),
                np.array([c[o + 0], c[o + 1], c[o + 3]]))
    a2, a1, a0 = coeffs(0)
    b2, b1, b0 = coeffs(6)
    mul = npoly.polymul
    sub = lambda p, q: npoly.polysub(p, q)
    # Sylvester resultant of two quadratics in y
    res = sub(mul(sub(mul(a2, b0), mul(b2, a0)), sub(mul(a2, b0), mul(b2, a0))),
              mul(sub(mul(a2, b1), mul(b2, a1)), sub(mul(a1, b0), mul(b1, a0))))
    res = np.trim_zeros(res, 'b')
    pts = []
    if len(res) > 1:
        for xr in np.roots(res[::-1]):
            if abs(xr.imag) > 1e-8 * max(1.0, abs(xr.real)):
                continue
            x = xr.real
            # solve P=0 for y, keep roots that also kill Q
            pa = np.array([c[5], c[2] + c[4] * x, c[0] + c[1] * x + c[3] * x * x])
            ys = np.roots(pa) if abs(pa[0]) > 1e-14 else (
                np.array([-pa[2] / pa[1]]) if abs(pa[1]) > 1e-14 else np.array([]))
            for yr in np.atleast_1d(ys):
                if abs(np.imag(yr)) > 1e-8 * max(1.0, abs(np.real(yr))):
                    continue
                y = float(np.real(yr))
                p, q = residual_at(c, (x, y))
                sc = max(1.0, abs(x), abs(y)) ** 2
                if abs(p) < 1e-6 * sc and abs(q) < 1e-6 * sc:
                    pts.append(newton_eq(c, x, y))
    # dedupe
    keep = []
    for p in pts:
        if p is None:
            continue
        if not any(np.hypot(p[0] - k[0], p[1] - k[1]) < 1e-7 * max(1, abs(p[0]), abs(p[1]))
                   for k in keep):
            keep.append(p)
    return keep


def jac(vec12, pt):
    c = np.asarray(vec12, float)
    x, y = pt
    return np.array([[c[1] + 2 * c[3] * x + c[4] * y, c[2] + c[4] * x + 2 * c[5] * y],
                     [c[7] + 2 * c[9] * x + c[10] * y, c[8] + c[10] * x + 2 * c[11] * y]])


def newton_eq(vec12, x, y, iters=60):
    for _ in range(iters):
        p, q = residual_at(vec12, (x, y))
        J = jac(vec12, (x, y))
        d = np.linalg.det(J)
        if abs(d) < 1e-300:
            return None
        dx = (-p * J[1, 1] + q * J[0, 1]) / d
        dy = (-q * J[0, 0] + p * J[1, 0]) / d
        x += dx; y += dy
        if abs(dx) + abs(dy) < 1e-15 * (1 + abs(x) + abs(y)):
            break
    return (x, y)


def foci(vec12):
    """Equilibria that are foci or centres (complex eigenvalues)."""
    out = []
    for p in equilibria(vec12):
        J = jac(vec12, p)
        tr, dt = J[0, 0] + J[1, 1], np.linalg.det(J)
        if dt > 0 and tr * tr - 4 * dt < 0:
            out.append((p, tr, dt))
    return out


# ----------------------------------------------------------------- calls
def d_curve(loc10, phi, s, b=0.0, rtol=None, **kw):
    o = dict(DEFAULTS); o.update(kw)
    if rtol is not None:
        o["rtol"] = rtol
    s = np.ascontiguousarray(np.asarray(s, float))
    n = s.size
    D = np.empty(n); st = np.empty(n, np.int32); T = np.empty(n)
    _lib.d_curve(np.ascontiguousarray(loc10), float(phi), s, n, float(b),
                 o["rtol"], o["Tmax"], o["Rmax"], o["nstep"], D, st, T)
    return D, st, T


def d_curve_noisy(loc10, phi, s, b=0.0, **kw):
    """Two-tolerance displacement (PROTOCOL rule 1).
    noise = 10*|D(rtol) - D(10*rtol)| + 5e-12*s."""
    o = dict(DEFAULTS); o.update(kw)
    D1, st1, T = d_curve(loc10, phi, s, b, rtol=o["rtol"], **{k: v for k, v in o.items() if k != "rtol"})
    D2, st2, _ = d_curve(loc10, phi, s, b, rtol=o["rtol"] * 10, **{k: v for k, v in o.items() if k != "rtol"})
    s = np.asarray(s, float)
    noise = np.where((st1 == 0) & (st2 == 0), 10 * np.abs(D1 - D2) + 5e-12 * s, np.inf)
    return D1, st1, noise, T


def count_sign_changes(s, D, st, noise):
    """PROTOCOL rule 1: a bracket counts only if both endpoints are above noise."""
    br = []
    ok = (st == 0) & np.isfinite(D)
    idx = np.where(ok)[0]
    for a, b in zip(idx[:-1], idx[1:]):
        if b != a + 1:
            continue                      # never bracket across a failure (bug B2)
        if D[a] * D[b] < 0 and min(abs(D[a]), abs(D[b])) > max(noise[a], noise[b]):
            br.append((s[a], s[b], D[a], D[b]))
    return br


def betastar(loc10, phi, s, dirhint=0, bmax=1.5, btol=1e-10, **kw):
    o = dict(DEFAULTS); o.update(kw)
    s = np.ascontiguousarray(np.asarray(s, float))
    n = s.size
    b = np.empty(n); st = np.empty(n, np.int32); d0 = np.empty(n); nf = np.empty(n)
    _lib.betastar(np.ascontiguousarray(loc10), float(phi), s, n,
                  o["rtol"], o["Tmax"], o["Rmax"], o["nstep"],
                  float(bmax), float(btol), int(dirhint), b, st, d0, nf)
    return b, st, d0, nf


def rotation_direction(loc10, phi, s_probe, **kw):
    """Global sign of dD/db (constant in s for a rotated family)."""
    eps = 1e-4
    sp = np.array([s_probe], float)
    Dp, stp, _ = d_curve(loc10, phi, sp, +eps, **kw)
    Dm, stm, _ = d_curve(loc10, phi, sp, -eps, **kw)
    if stp[0] or stm[0]:
        return 0
    return 1 if Dp[0] > Dm[0] else -1


# ------------------------------------------------------- extrema of beta*
def interior_extrema(s, b, st, min_prom_rel=1e-9):
    """Interior local extrema of beta*(s) on the resolved domain.

    Returns (list of (index, kind, prominence), height_range).  Prominence is
    measured as the smaller of the two drops to the neighbouring extrema /
    endpoints, normalised by the total height range of the curve.
    """
    ok = (st == 0) & np.isfinite(b)
    idx = np.where(ok)[0]
    if idx.size < 5:
        return [], 0.0
    # longest contiguous run
    runs, cur = [], [idx[0]]
    for i in idx[1:]:
        if i == cur[-1] + 1:
            cur.append(i)
        else:
            runs.append(cur); cur = [i]
    runs.append(cur)
    run = max(runs, key=len)
    if len(run) < 5:
        return [], 0.0
    bb = b[run]
    rng = bb.max() - bb.min()
    if rng <= 0:
        return [], 0.0
    ext = []
    for j in range(1, len(bb) - 1):
        if (bb[j] - bb[j - 1]) * (bb[j + 1] - bb[j]) < 0:
            kind = "max" if bb[j] > bb[j - 1] else "min"
            ext.append([j, kind, 0.0])
    # prominence: distance in b to the nearest neighbouring turning point
    marks = [0] + [e[0] for e in ext] + [len(bb) - 1]
    for k, e in enumerate(ext):
        j = e[0]
        left, right = marks[k], marks[k + 2]
        e[2] = min(abs(bb[j] - bb[left]), abs(bb[j] - bb[right])) / rng
    ext = [(run[e[0]], e[1], e[2]) for e in ext if e[2] > min_prom_rel]
    return ext, rng
