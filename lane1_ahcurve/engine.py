"""lane1_ahcurve/engine.py -- driver for the Andronov-Hopf curve engine.

Coefficient convention (12-vector), identical to the repository's:
    P = c0 + c1 x + c2 y + c3 x^2 + c4 x y + c5 y^2
    Q = c6 + c7 x + c8 y + c9 x^2 + c10 x y + c11 y^2
    xdot = P, ydot = Q.

Everything the C engine sees is LOCAL: the field re-expanded about the focus,
so c0 = c6 = 0 and the orbit radius is the coordinate magnitude.
"""
import ctypes, hashlib, os, math
from fractions import Fraction
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.join(_HERE, "ahcurve.c")
_LIB = os.path.join(_HERE, "libahcurve.so")

def engine_hash():
    with open(_SRC, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()[:16]

ENGINE = "lane1_ahcurve/ahcurve.c@" + engine_hash()

_l = ctypes.CDLL(_LIB)
_dp = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
_ip = np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS")

_l.returns_rot.argtypes = [ctypes.c_int,_dp,_dp,_dp,ctypes.c_int,_dp,_dp,_dp,_ip,
                           ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_long]
_l.returns_lin.argtypes = [ctypes.c_int,_dp,_dp,_dp,_dp,ctypes.c_int,_dp,_dp,_dp,_ip,
                           ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_long]
_l.ahcurve.argtypes = [ctypes.c_int,_dp,ctypes.c_void_p,_dp,_dp,ctypes.c_int,
                       ctypes.c_int,_dp,_dp,_ip,_ip,
                       ctypes.c_double,ctypes.c_double,
                       ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_long]

DEFAULTS = dict(rtol=1e-12, Rmax=1e6, Tmax=2.0e5, maxsteps=4_000_000)

# ----------------------------------------------------------------- algebra --

def local_expand(c, x0, y0):
    """Re-expand P,Q about (x0,y0).  Exact when the inputs are Fractions."""
    out = []
    for base in (0, 6):
        c0, c1, c2, c3, c4, c5 = c[base:base+6]
        C0 = c0 + c1*x0 + c2*y0 + c3*x0*x0 + c4*x0*y0 + c5*y0*y0
        C1 = c1 + 2*c3*x0 + c4*y0
        C2 = c2 + c4*x0 + 2*c5*y0
        out += [C0, C1, C2, c3, c4, c5]
    return out

def jacobian_at_origin(cl):
    """Linear part of a LOCAL field: [[P_u,P_v],[Q_u,Q_v]]."""
    return ((cl[1], cl[2]), (cl[7], cl[8]))

def hopf_angle(cl):
    """b with trace(R_b L) = 0, i.e. the angle at which the focus is weak.
    trace = (a+d) cos b + (b'-c) sin b  ->  b = atan2(-(a+d), (b'-c))."""
    (a, bb), (cc, d) = jacobian_at_origin(cl)
    b = math.atan2(-(a + d), (bb - cc))
    # trace(b) is pi-periodic in sign structure; take the representative in
    # (-pi/2, pi/2] so that b = 0 is the identity member of the family.
    while b > math.pi/2:
        b -= math.pi
    while b <= -math.pi/2:
        b += math.pi
    return b

def is_focus(cl):
    (a, bb), (cc, d) = jacobian_at_origin(cl)
    tr, det = a + d, a*d - bb*cc
    return det > 0 and tr*tr - 4*det < 0

def equilibria(c, tol=1e-9):
    """Real solutions of P=Q=0 by resultant-free numeric means (companion of a
    2x2 polynomial system via numpy roots on the resultant in x)."""
    import numpy.polynomial.polynomial as _P
    P = [float(v) for v in c[0:6]]
    Q = [float(v) for v in c[6:12]]
    # P = (p5) y^2 + (p2 + p4 x) y + (p0 + p1 x + p3 x^2)
    # Q = (q5) y^2 + (q2 + q4 x) y + (q0 + q1 x + q3 x^2)
    out = []
    xs = np.linspace(-1e3, 1e3, 1)  # placeholder, replaced below
    # Sylvester resultant in y of two quadratics in y with coefficients poly in x
    def coeffs(vec, x):
        return (vec[5], vec[2] + vec[4]*x, vec[0] + vec[1]*x + vec[3]*x*x)
    # build resultant as a polynomial in x by sampling + fit (degree <= 4)
    deg = 4
    xs = np.cos(np.pi*(np.arange(deg+1)+0.5)/(deg+1))*10.0
    vals = []
    for x in xs:
        a2, a1, a0 = coeffs(P, x)
        b2, b1, b0 = coeffs(Q, x)
        M = np.array([[a2, a1, a0, 0.0],
                      [0.0, a2, a1, a0],
                      [b2, b1, b0, 0.0],
                      [0.0, b2, b1, b0]])
        vals.append(np.linalg.det(M))
    co = np.polyfit(xs, vals, deg)
    if np.max(np.abs(co)) < 1e-14:
        return []
    rts = np.roots(co)
    for r in rts:
        if abs(r.imag) > 1e-7*max(1.0, abs(r.real)):
            continue
        x = r.real
        a2, a1, a0 = coeffs(P, x)
        cand = []
        if abs(a2) > 1e-12:
            d = a1*a1 - 4*a2*a0
            if d >= 0:
                sq = math.sqrt(d)
                cand += [(-a1+sq)/(2*a2), (-a1-sq)/(2*a2)]
        elif abs(a1) > 1e-12:
            cand.append(-a0/a1)
        for y in cand:
            px = P[0]+P[1]*x+P[2]*y+P[3]*x*x+P[4]*x*y+P[5]*y*y
            qx = Q[0]+Q[1]*x+Q[2]*y+Q[3]*x*x+Q[4]*x*y+Q[5]*y*y
            sc = max(1.0, abs(x), abs(y))**2
            if abs(px) < 1e-6*sc and abs(qx) < 1e-6*sc:
                if not any(abs(x-u) < 1e-7*max(1,abs(x)) and abs(y-v) < 1e-7*max(1,abs(y))
                           for u, v in out):
                    out.append((x, y))
    return out

# --------------------------------------------------------------- the engine --

def _arr(a):
    return np.ascontiguousarray(np.asarray(a, dtype=np.float64))

def returns(cloc, direction, svals, b=0.0, evec=None, t=None, **kw):
    """Return radii R(s) on one field member.  Rotation family unless evec given."""
    o = dict(DEFAULTS); o.update(kw)
    s = _arr(svals); ns = s.size
    C = _arr(cloc).reshape(1, 12)
    D = _arr(direction).reshape(1, 2)
    R = np.empty(ns); T = np.empty(ns); st = np.empty(ns, dtype=np.int32)
    if evec is None:
        _l.returns_rot(1, C, D, _arr([b]), ns, s, R, T, st,
                       o["rtol"], o["Rmax"], o["Tmax"], o["maxsteps"])
    else:
        E = _arr(evec).reshape(1, 12)
        _l.returns_lin(1, C, E, D, _arr([t]), ns, s, R, T, st,
                       o["rtol"], o["Rmax"], o["Tmax"], o["maxsteps"])
    return R, T, st

def displacement(cloc, direction, svals, **kw):
    R, T, st = returns(cloc, direction, svals, **kw)
    D = np.where(st == 0, R - np.asarray(svals, dtype=float), np.nan)
    return D, st

def ah_curve(cloc, direction, svals, mode="rot", evec=None, p0=None,
             ptol=1e-10, span=1.2, **kw):
    """beta*(s): the family parameter that closes the orbit through each s."""
    o = dict(DEFAULTS); o.update(kw)
    s = _arr(svals); ns = s.size
    C = _arr(cloc).reshape(1, 12)
    Dv = _arr(direction).reshape(1, 2)
    if p0 is None:
        p0 = hopf_angle(cloc) if mode == "rot" else 0.0
    P0 = _arr([p0])
    out = np.empty(ns); st = np.empty(ns, dtype=np.int32); nev = np.empty(ns, dtype=np.int32)
    if mode == "rot":
        ep = None
        _l.ahcurve(1, C, None, Dv, P0, 0, ns, s, out, st, nev,
                   ptol, span, o["rtol"], o["Rmax"], o["Tmax"], o["maxsteps"])
    else:
        E = _arr(evec).reshape(1, 12)
        _l.ahcurve(1, C, E.ctypes.data_as(ctypes.c_void_p), Dv, P0, 1, ns, s, out, st, nev,
                   ptol, span, o["rtol"], o["Rmax"], o["Tmax"], o["maxsteps"])
    return out, st, nev

def ah_curve_batch(clocs, dirs, svals, p0s, mode="rot", evecs=None,
                   ptol=1e-10, span=1.2, **kw):
    """Many fields at once; svals is (nsets, ns)."""
    o = dict(DEFAULTS); o.update(kw)
    n = len(clocs); s = _arr(svals).reshape(n, -1); ns = s.shape[1]
    C = _arr(clocs).reshape(n, 12)
    Dv = _arr(dirs).reshape(n, 2)
    P0 = _arr(p0s).reshape(n)
    out = np.empty(n*ns); st = np.empty(n*ns, dtype=np.int32); nev = np.empty(n*ns, dtype=np.int32)
    if mode == "rot":
        _l.ahcurve(n, C, None, Dv, P0, 0, ns, s.ravel(), out, st, nev,
                   ptol, span, o["rtol"], o["Rmax"], o["Tmax"], o["maxsteps"])
    else:
        E = _arr(evecs).reshape(n, 12)
        _l.ahcurve(n, C, E.ctypes.data_as(ctypes.c_void_p), Dv, P0, 1, ns, s.ravel(), out, st, nev,
                   ptol, span, o["rtol"], o["Rmax"], o["Tmax"], o["maxsteps"])
    return out.reshape(n, ns), st.reshape(n, ns), nev.reshape(n, ns)

# ---------------------------------------------------- noise, extrema, counts --

def s_max(cloc, direction, s_lo, s_hi, b=0.0, n=60, **kw):
    """Largest s on a geometric probe for which the return succeeds, and the
    smallest failing s.  PROTOCOL rule 4: never count beyond this."""
    s = np.geomspace(s_lo, s_hi, n)
    R, T, st = returns(cloc, direction, s, b=b, **kw)
    ok = st == 0
    if not ok.any():
        return None, s[0]
    last = np.max(np.nonzero(ok)[0])
    fail = s[last+1] if last+1 < n else None
    return s[last], fail

def curve_with_noise(cloc, direction, svals, mode="rot", evec=None, p0=None,
                     ptol=1e-10, span=1.2, rtol=1e-12, rtol_loose=1e-10, floor=None, **kw):
    """PROTOCOL rule 1 applied to beta*: two-tolerance differential noise."""
    if floor is None:
        floor = 2*ptol
    b1, st1, ne1 = ah_curve(cloc, direction, svals, mode, evec, p0, ptol, span,
                            rtol=rtol, **kw)
    b2, st2, ne2 = ah_curve(cloc, direction, svals, mode, evec, p0, ptol, span,
                            rtol=rtol_loose, **kw)
    ok = (st1 == 0) & (st2 == 0)
    noise = np.where(ok, 10*np.abs(b1-b2) + floor, np.inf)
    return b1, st1, noise, int(ne1.sum()+ne2.sum())

def turning_points(y, tau):
    """Interior extrema of a sampled curve, with hysteresis tau.
    Returns [(index, +1 for max / -1 for min), ...]."""
    y = np.asarray(y, dtype=float)
    n = y.size
    if n < 3:
        return []
    ext = []
    d = 0
    ci = 0
    lo_i = hi_i = 0
    for i in range(1, n):
        v = y[i]
        if not np.isfinite(v):
            continue
        if d == 0:
            if v < y[lo_i]:
                lo_i = i
            if v > y[hi_i]:
                hi_i = i
            if v > y[lo_i] + tau:
                d = 1; ci = lo_i
                # we were descending into lo_i then rose: lo_i is interior min
                if lo_i > 0:
                    ext.append((lo_i, -1))
                ci = i
            elif v < y[hi_i] - tau:
                d = -1
                if hi_i > 0:
                    ext.append((hi_i, +1))
                ci = i
        elif d == 1:
            if v >= y[ci]:
                ci = i
            elif v < y[ci] - tau:
                ext.append((ci, +1)); d = -1; ci = i
        else:
            if v <= y[ci]:
                ci = i
            elif v > y[ci] + tau:
                ext.append((ci, -1)); d = 1; ci = i
    return [(i, sg) for (i, sg) in ext if 0 < i < n-1]

def count_sign_changes(s, D, noise):
    """PROTOCOL rule 1: a cycle only from a sign change of D whose two endpoint
    magnitudes both clear the two-tolerance noise estimate."""
    s = np.asarray(s, float); D = np.asarray(D, float); noise = np.asarray(noise, float)
    br = []
    for i in range(len(D)-1):
        if not (np.isfinite(D[i]) and np.isfinite(D[i+1])):
            continue
        if D[i]*D[i+1] < 0 and min(abs(D[i]), abs(D[i+1])) > max(noise[i], noise[i+1]):
            br.append((s[i], s[i+1]))
    return br

def displacement_with_noise(cloc, direction, svals, b=0.0, rtol=1e-12, rtol_loose=1e-10, **kw):
    D1, st1 = displacement(cloc, direction, svals, b=b, rtol=rtol, **kw)
    D2, st2 = displacement(cloc, direction, svals, b=b, rtol=rtol_loose, **kw)
    ok = (st1 == 0) & (st2 == 0)
    s = np.asarray(svals, float)
    noise = np.where(ok, 10*np.abs(D1-D2) + 5e-12*s, np.inf)
    return D1, st1, noise
