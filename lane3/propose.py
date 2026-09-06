"""Floating point PROPOSAL layer for the Lane 3 certifier.

Nothing in this file is part of any proof.  It locates limit cycles
numerically and proposes rational polygons; every proposal is then decided in
exact rational arithmetic by certify.py / verify.py, which never call anything
here.

Construction of a transversal polygon
-------------------------------------
Work in polar coordinates about a rational centre O inside the cycle.  Because
limit cycles of a quadratic field are convex and enclose exactly one
antisaddle, the polar angle is strictly monotone along the orbit near the
cycle, so an orbit is a graph r(theta).

Let the cycle be at radius r* on the ray theta = theta0 and let the direction
of increasing theta be the contracting one.  Integrate one orbit starting at
r0 = r* + delta for two full turns, giving

    r_a(tau) = r(theta0 + 2 pi tau),      r_b(tau) = r(theta0 + 2 pi (1 + tau))

with r_b strictly between r_a and the cycle.  The closed curve

    R(tau) = (1 - tau) r_b(tau) + tau r_a(tau),   tau in [0,1]

satisfies R(0) = R(1) = r(theta0 + 2 pi), so it is closed, and it drifts from
the inner turn to the outer turn, i.e. AWAY from the cycle, at rate
|r_a - r_b| per turn.  The flow is tangent to each turn, so relative to R the
flow has normal speed about |r_a - r_b| / T pointing TOWARDS the cycle, and the
same construction started at r0 = r* - delta produces the companion curve on
the other side.  Both boundaries are then crossed towards the cycle, which is
exactly the Poincare-Bendixson annulus condition.

The transversality margin is proportional to delta * |1 - multiplier|, so
cycles whose multiplier is very close to 1 need many, short edges: the chord of
a curve of curvature kappa deviates in direction from the curve by about
h*kappa/2, which must stay below the crossing angle.
"""

import math
from fractions import Fraction as F

import numpy as np
from scipy.integrate import solve_ivp


# ------------------------------------------------------------------ field

def as_float(vec12):
    return [float(c) for c in vec12]


def PQ(v, x, y):
    p = v[0] + v[1] * x + v[2] * y + v[3] * x * x + v[4] * x * y + v[5] * y * y
    q = v[6] + v[7] * x + v[8] * y + v[9] * x * x + v[10] * x * y + v[11] * y * y
    return p, q


def jac(v, x, y):
    return np.array([[v[1] + 2 * v[3] * x + v[4] * y, v[2] + v[4] * x + 2 * v[5] * y],
                     [v[7] + 2 * v[9] * x + v[10] * y, v[8] + v[10] * x + 2 * v[11] * y]])


def equilibria(v):
    """All finite equilibria of the quadratic field, numerically (via the
    resultant in y of the two conics)."""
    import numpy.polynomial.polynomial as npp
    # P = (v3) x^2 + (v1 + v4 y) x + (v0 + v2 y + v5 y^2)  -- as poly in x
    def coeffs_in_x(c):
        return [c[3], None, None]
    out = []
    # Sylvester resultant in x of P and Q, both quadratic in x
    ys = None
    A2, A1, A0 = v[3], np.poly1d([v[4], v[1]]), np.poly1d([v[5], v[2], v[0]])
    B2, B1, B0 = v[9], np.poly1d([v[10], v[7]]), np.poly1d([v[11], v[8], v[6]])
    A2 = np.poly1d([v[3]])
    B2 = np.poly1d([v[9]])
    # Sylvester 4x4 determinant with polynomial entries in y
    M = [[A2, A1, A0, np.poly1d([0])],
         [np.poly1d([0]), A2, A1, A0],
         [B2, B1, B0, np.poly1d([0])],
         [np.poly1d([0]), B2, B1, B0]]
    det = _poly_det4(M)
    roots = det.r if det.order > 0 else np.array([])
    for yr in roots:
        if abs(yr.imag) > 1e-8:
            continue
        yv = yr.real
        pa, pb, pc = v[3], v[1] + v[4] * yv, v[0] + v[2] * yv + v[5] * yv * yv
        cands = np.roots([pa, pb, pc]) if abs(pa) > 1e-14 else (
            np.array([-pc / pb]) if abs(pb) > 1e-14 else np.array([]))
        for xr in np.atleast_1d(cands):
            if abs(np.imag(xr)) > 1e-8:
                continue
            xv = float(np.real(xr))
            p, q = PQ(v, xv, yv)
            if abs(p) < 1e-6 and abs(q) < 1e-6:
                # polish
                z = np.array([xv, yv])
                for _ in range(60):
                    p, q = PQ(v, z[0], z[1])
                    J = jac(v, z[0], z[1])
                    try:
                        z = z - np.linalg.solve(J, np.array([p, q]))
                    except np.linalg.LinAlgError:
                        break
                if not any(abs(z[0] - e[0]) + abs(z[1] - e[1]) < 1e-9 for e in out):
                    out.append((z[0], z[1]))
    return out


def _poly_det4(M):
    n = 4
    def det(rows, cols):
        if len(cols) == 1:
            return M[rows[0]][cols[0]]
        tot = np.poly1d([0.0])
        for k, c in enumerate(cols):
            sub = det(rows[1:], [cc for cc in cols if cc != c])
            term = M[rows[0]][c] * sub
            tot = tot + (term if k % 2 == 0 else -term)
        return tot
    return det(list(range(n)), list(range(n)))


def classify(v, e):
    J = jac(v, e[0], e[1])
    tr, dt = J[0, 0] + J[1, 1], np.linalg.det(J)
    disc = tr * tr - 4 * dt
    if dt < 0:
        return "saddle", tr, dt
    if disc < 0:
        return "focus", tr, dt
    return "node", tr, dt


# --------------------------------------------------- polar orbit machinery

def make_rhs(v, O):
    ox, oy = O
    def rhs(th, r):
        rr = r[0]
        c, s = math.cos(th), math.sin(th)
        x, y = ox + rr * c, oy + rr * s
        p, q = PQ(v, x, y)
        num = p * c + q * s
        den = -p * s + q * c
        return [rr * num / den]
    return rhs


def integrate_turns(v, O, th0, r0, nturns, dsign=1, rtol=1e-12, atol=1e-14):
    rhs = make_rhs(v, O)
    th1 = th0 + dsign * 2.0 * math.pi * nturns
    sol = solve_ivp(rhs, (th0, th1), [r0], method="DOP853", dense_output=True,
                    rtol=rtol, atol=atol, max_step=abs(th1 - th0) / 200.0)
    if not sol.success:
        return None
    return sol


def return_map(v, O, th0, r0, dsign=1, rtol=1e-12):
    sol = integrate_turns(v, O, th0, r0, 1, dsign, rtol)
    if sol is None:
        return None
    return float(sol.y[0, -1])


def displacement(v, O, th0, r0, dsign=1, rtol=1e-12):
    g = return_map(v, O, th0, r0, dsign, rtol)
    return None if g is None else g - r0


def scan_cycles(v, O, th0, rlo, rhi, n=200, dsign=1, log=True):
    """Bracket sign changes of the displacement on the ray theta = th0."""
    rs = np.geomspace(rlo, rhi, n) if log else np.linspace(rlo, rhi, n)
    ds = []
    for r in rs:
        try:
            ds.append(displacement(v, O, th0, float(r), dsign))
        except Exception:
            ds.append(None)
    brackets = []
    for i in range(len(rs) - 1):
        a, b = ds[i], ds[i + 1]
        if a is None or b is None:
            continue
        if a == 0 or (a < 0) != (b < 0):
            brackets.append((float(rs[i]), float(rs[i + 1]), a, b))
    return rs, ds, brackets


def refine_cycle(v, O, th0, ra, rb, dsign=1, iters=90):
    fa = displacement(v, O, th0, ra, dsign)
    for _ in range(iters):
        rm = 0.5 * (ra + rb)
        fm = displacement(v, O, th0, rm, dsign)
        if fm is None:
            break
        if (fa < 0) == (fm < 0):
            ra, fa = rm, fm
        else:
            rb = rm
        if abs(rb - ra) <= 1e-15 * max(1.0, abs(ra)):
            break
    return 0.5 * (ra + rb)


def cycle_multiplier(v, O, th0, rstar, dsign=1, h=None):
    """dG/dr at the cycle (the Floquet multiplier) by a centred difference."""
    if h is None:
        h = max(1e-7, 1e-7 * abs(rstar))
    gp = return_map(v, O, th0, rstar + h, dsign)
    gm = return_map(v, O, th0, rstar - h, dsign)
    if gp is None or gm is None:
        return None
    return (gp - gm) / (2 * h)


# ------------------------------------------------------- curve proposal

def blend_curve(v, O, th0, rstar, delta, dsign, ntau, rtol=1e-13):
    """Sample the closed blended curve at ntau angles.  Returns (thetas, radii)
    or None if the two-turn integration failed."""
    sol = integrate_turns(v, O, th0, rstar + delta, 2, dsign, rtol)
    if sol is None:
        return None
    taus = np.linspace(0.0, 1.0, ntau + 1)[:-1]
    th_a = th0 + dsign * 2 * math.pi * taus
    th_b = th0 + dsign * 2 * math.pi * (1.0 + taus)
    ra = np.array([float(sol.sol(t)[0]) for t in th_a])
    rb = np.array([float(sol.sol(t)[0]) for t in th_b])
    R = (1.0 - taus) * rb + taus * ra
    return th_a, R, ra, rb


def curve_eval(sol, th0, dsign, tau):
    """Blended radius at a single tau, from an existing two-turn solution."""
    th_a = th0 + dsign * 2 * math.pi * tau
    th_b = th0 + dsign * 2 * math.pi * (1.0 + tau)
    ra = float(sol.sol(th_a)[0])
    rb = float(sol.sol(th_b)[0])
    return th_a, (1.0 - tau) * rb + tau * ra


def dyadic(x, bits):
    """Nearest rational with denominator 2**bits."""
    d = 1 << bits
    return F(int(round(x * d)), d)


# ----------------------------------------------- phase-parametrised blending
#
# The angle-parametrised blend of the two turns has a gap vector that is radial;
# where the ray from O is oblique to the orbit its normal component collapses
# and the transversality margin with it.  Blending at matched PHASE instead
# keeps the gap vector along the normal deviation between successive turns,
# which is where the whole margin lives.  Closure of the blended curve is exact
# by construction: both ends of the blend are the same trajectory point Z(T1).

def period_estimate(v, O, th0, r0, dsign=1):
    """Time for one turn about O, from the polar system augmented with dt/dtheta."""
    ox, oy = O

    def rhs(th, s):
        rr = s[0]
        c, sn = math.cos(th), math.sin(th)
        x, y = ox + rr * c, oy + rr * sn
        p, q = PQ(v, x, y)
        den = -p * sn + q * c
        return [rr * (p * c + q * sn) / den, rr / den]

    th1 = th0 + dsign * 2 * math.pi
    sol = solve_ivp(rhs, (th0, th1), [r0, 0.0], method="DOP853",
                    rtol=1e-11, atol=1e-13)
    if not sol.success:
        return None
    return abs(float(sol.y[1, -1]))


def two_turn_orbit(v, O, th0, r0, dsign=1, rtol=1e-13, atol=1e-15):
    """Integrate the orbit through O + r0*(cos th0, sin th0) for two returns to
    that ray, in the time direction for which `dsign` is the direction of
    increasing polar angle.  Returns (sol, T1, T2, time_dir) with sol a dense
    solution in TIME, T1 the first return time and T2 the second."""
    ox, oy = O
    c0, s0 = math.cos(th0), math.sin(th0)
    z0 = [ox + r0 * c0, oy + r0 * s0]
    p, q = PQ(v, z0[0], z0[1])
    thdot = (-p * s0 + q * c0)
    if thdot == 0:
        return None
    time_dir = dsign * (1.0 if thdot > 0 else -1.0)
    Test = period_estimate(v, O, th0, r0, dsign)
    if Test is None or not math.isfinite(Test):
        return None
    Tmax = 2.6 * Test

    def rhs(t, z):
        pp, qq = PQ(v, z[0], z[1])
        return [time_dir * pp, time_dir * qq]

    def ev(t, z):
        return c0 * (z[1] - oy) - s0 * (z[0] - ox)
    ev.terminal = False
    ev.direction = 0

    sol = solve_ivp(rhs, (0.0, Tmax), z0, method="DOP853", dense_output=True,
                    events=[ev], rtol=rtol, atol=atol,
                    max_step=Tmax / 400.0)
    if not sol.success:
        return None
    hits = []
    for t in sol.t_events[0]:
        if t <= 1e-9 * Tmax:
            continue
        z = sol.sol(t)
        if c0 * (z[0] - ox) + s0 * (z[1] - oy) > 0:
            hits.append(float(t))
    if len(hits) < 2:
        return None
    return sol, hits[0], hits[1], time_dir


def blend_points(sol, T1, T2, taus):
    """Points of the closed blended curve C(tau) = (1-tau) Z(T1 + tau (T2-T1))
    + tau Z(tau T1) at the given phases.  C(0) = C(1) = Z(T1)."""
    out = []
    for tau in taus:
        za = sol.sol(tau * T1)
        zb = sol.sol(T1 + tau * (T2 - T1))
        out.append(((1.0 - tau) * zb[0] + tau * za[0],
                    (1.0 - tau) * zb[1] + tau * za[1]))
    return out
