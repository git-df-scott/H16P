"""Trigger-precision recheck (PROTOCOL rules 3 and 6).

Every bracket endpoint that would decide a trigger is recomputed in binary128
at two tolerances, and once more with mpmath at dps >= 40 through an
independent Taylor-series integrator, so that the sign of D at the endpoint
rests on three separate computations.
"""
import numpy as np
import engine as E


def quad_pair(L, phi, s, b):
    """binary128 D at two tolerances -> (D, noise estimate)."""
    D1, st1, _ = E.d_curve_quad(L, phi, np.array([s]), b, rtol=1e-18)
    D2, st2, _ = E.d_curve_quad(L, phi, np.array([s]), b, rtol=1e-16)
    if st1[0] or st2[0]:
        return None, None, int(st1[0]), int(st2[0])
    return float(D1[0]), 10 * abs(float(D1[0] - D2[0])), 0, 0


def recheck(L, phi, b, brackets):
    out = []
    for (s1, s2, D1, D2) in brackets:
        row = {}
        for tag, s, Dd in (("lo", s1, D1), ("hi", s2, D2)):
            v, nz, k1, k2 = quad_pair(L, phi, float(s), float(b))
            row[tag] = dict(s=float(s), D_double=float(Dd), D_binary128=v,
                            noise_binary128=nz, status=(k1, k2))
        row["sign_change"] = (
            row["lo"]["D_binary128"] is not None
            and row["hi"]["D_binary128"] is not None
            and row["lo"]["D_binary128"] * row["hi"]["D_binary128"] < 0
            and abs(row["lo"]["D_binary128"]) > (row["lo"]["noise_binary128"] or 0)
            and abs(row["hi"]["D_binary128"]) > (row["hi"]["noise_binary128"] or 0))
        out.append(row)
    return out


# ------------------------------------------------------------------ mpmath
def mp_return(L, phi, s, b, dps=45, Tmax=400.0):
    """Independent arbitrary-precision return via mpmath's Taylor ODE solver."""
    from mpmath import mp, mpf, cos, sin, sqrt, atan2, odefun, findroot
    mp.dps = dps
    p1, p2, p3, p4, p5, q1, q2, q3, q4, q5 = [mpf(repr(float(x))) for x in L]
    bb = mpf(repr(float(b)))
    cb, sb = cos(bb), sin(bb)
    ce, se = cos(mpf(repr(float(phi)))), sin(mpf(repr(float(phi))))

    def f(t, y):
        u, v = y[0], y[1]
        P = p1 * u + p2 * v + p3 * u * u + p4 * u * v + p5 * v * v
        Q = q1 * u + q2 * v + q3 * u * u + q4 * u * v + q5 * v * v
        fx = P * cb - Q * sb
        fy = P * sb + Q * cb
        return [fx, fy, (u * fy - v * fx) / (u * u + v * v)]

    ss = mpf(repr(float(s)))
    y0 = [ss * ce, ss * se, mpf(0)]
    g = odefun(f, mpf(0), y0, tol=mpf(10) ** (-dps + 6), degree=None)
    h = lambda t: -se * g(t)[0] + ce * g(t)[1]
    th = lambda t: g(t)[2]

    # march on a coarse grid to bracket the return crossing
    t, dt = mpf(0), mpf("0.02")
    prev = h(t)
    have = False
    while t < Tmax:
        t2 = t + dt
        cur = h(t2)
        if not have and abs(cur) > mpf(10) ** (-dps + 10):
            have = True
        elif have and prev * cur < 0:
            tc = findroot(h, (t, t2), solver="anderson", tol=mpf(10) ** (-2 * dps + 20))
            y = g(tc)
            proj = y[0] * ce + y[1] * se
            if proj > 0 and abs(y[2]) > 3:
                return proj - ss, tc, y[2]
        t, prev = t2, cur
    return None, None, None
