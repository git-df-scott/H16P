"""Independent second integrator (PROTOCOL rule 2).

Different method (SciPy DOP853 / Radau, dense event location) and a different
crossing algorithm (SciPy's Brentq on the interpolant) from the C engine.
Same local-chart re-expansion, so it tests the integrator and the crossing
logic, not the algebra.
"""
import numpy as np
from scipy.integrate import solve_ivp

TWO_PI = 2 * np.pi


def make_rhs(loc10, b):
    p1, p2, p3, p4, p5, q1, q2, q3, q4, q5 = loc10
    cb, sb = np.cos(b), np.sin(b)

    def f(t, y):
        u, v = y[0], y[1]
        P = p1 * u + p2 * v + p3 * u * u + p4 * u * v + p5 * v * v
        Q = q1 * u + q2 * v + q3 * u * u + q4 * u * v + q5 * v * v
        fx = P * cb - Q * sb
        fy = P * sb + Q * cb
        return [fx, fy, (u * fy - v * fx) / (u * u + v * v)]
    return f


def ret_once(loc10, phi, s, b=0.0, rtol=1e-12, atol_rel=1e-16,
             Tmax=400.0, Rmax=1e4, method="DOP853"):
    """Returns (R, status).  status 0 ok."""
    ce, se = np.cos(phi), np.sin(phi)
    f = make_rhs(loc10, b)

    def ev(t, y):
        return -se * y[0] + ce * y[1]
    ev.terminal = False
    ev.direction = 0

    def esc(t, y):
        return np.hypot(y[0], y[1]) - Rmax
    esc.terminal = True

    def wind(t, y):
        return abs(y[2]) - 2.4 * TWO_PI
    wind.terminal = True

    y0 = [s * ce, s * se, 0.0]
    sol = solve_ivp(f, (0.0, Tmax), y0, method=method, rtol=rtol,
                    atol=atol_rel * s + 1e-300, events=[ev, esc, wind],
                    dense_output=True)
    if not sol.success:
        return np.nan, 5
    # a valid full-turn crossing takes precedence over the terminal guards,
    # which only stop the integration after the return has already happened
    for te in sol.t_events[0]:
        y = sol.sol(te)
        proj = y[0] * ce + y[1] * se
        th = y[2]
        if proj > 0 and abs(th) > 3.0:
            if abs(abs(th) - TWO_PI) > 0.7:
                return np.nan, 4
            return float(proj), 0
    if sol.t_events[1].size:
        return np.nan, 1
    if sol.t_events[2].size:
        return np.nan, 4
    return np.nan, 2


def d_curve(loc10, phi, sarr, b=0.0, **kw):
    D = np.empty(len(sarr)); st = np.zeros(len(sarr), int)
    for i, s in enumerate(sarr):
        R, k = ret_once(loc10, phi, s, b, **kw)
        st[i] = k
        D[i] = R - s if k == 0 else np.nan
    return D, st
