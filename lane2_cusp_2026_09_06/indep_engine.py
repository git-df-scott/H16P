"""INDEPENDENT second engine (PROTOCOL rule 2): mpmath Taylor integrator for the
Cherkas normal form, written from the ODE directly with a different step-size rule,
a different event solver, and no shared code with cusp_engine.cpp.

Only D and D_x (by the standard variational equation) are provided; this engine
exists to cross-check the C++ engine's D, its root positions and its D_x, not to
carry the continuation.
"""
import mpmath as mp


def field(p, x, y):
    a, a20, a11, a01, a10, a00 = p
    return (1 + x * y, a00 + a10 * x + a20 * x * x + a01 * y + a11 * x * y + a * y * y)


def params(a, a20, a11, a01, a10):
    a, a20, a11, a01, a10 = map(mp.mpf, (a, a20, a11, a01, a10))
    a00 = a01 + a11 - a10 - a20 - a
    return (a, a20, a11, a01, a10, a00)


def taylor_coeffs(p, x0, y0, N):
    """Time-Taylor coefficients of the scalar (non-jet) solution."""
    a, a20, a11, a01, a10, a00 = p
    X = [mp.mpf(0)] * (N + 1)
    Y = [mp.mpf(0)] * (N + 1)
    X[0], Y[0] = x0, y0
    for k in range(N):
        XY = sum(X[j] * Y[k - j] for j in range(k + 1))
        X2 = sum(X[j] * X[k - j] for j in range(k + 1))
        Y2 = sum(Y[j] * Y[k - j] for j in range(k + 1))
        dx = XY + (1 if k == 0 else 0)
        dy = a10 * X[k] + a20 * X2 + a01 * Y[k] + a11 * XY + a * Y2 + (a00 if k == 0 else 0)
        X[k + 1] = dx / (k + 1)
        Y[k + 1] = dy / (k + 1)
    return X, Y


def poly(C, t):
    s = C[-1]
    for k in range(len(C) - 2, -1, -1):
        s = s * t + C[k]
    return s


def dpoly(C, t):
    s = mp.mpf(0)
    for k in range(len(C) - 1, 0, -1):
        s = s * t + k * C[k]
    return s


def ret(a, a20, a11, a01, a10, x0, side="right", N=24, tol=None, tmax=400, maxsteps=200000):
    """First return to {y=-1, side of x=1}.  Returns (R, T) or (None, reason)."""
    if tol is None:
        tol = mp.mpf(10) ** (-(mp.mp.dps - 6))
    p = params(a, a20, a11, a01, a10)
    x, y = mp.mpf(x0), mp.mpf(-1)
    t = mp.mpf(0)
    armed = False
    prev = None
    for step in range(maxsteps):
        X, Y = taylor_coeffs(p, x, y, N)
        mN = max(abs(X[N]), abs(Y[N]))
        mN1 = max(abs(X[N - 1]), abs(Y[N - 1]))
        h = mp.mpf("0.25")
        if mN > 0:
            h = min(h, (tol / mN) ** (mp.mpf(1) / N))
        if mN1 > 0:
            h = min(h, (tol / mN1) ** (mp.mpf(1) / (N - 1)))
        h *= mp.mpf("0.8")
        if h < mp.mpf("1e-13"):
            return None, "stepsize"
        lo = mp.mpf(0)
        if not armed:
            lo = h * mp.mpf("1e-9")
            prev = poly(Y, lo) + 1
            armed = True
        NSUB = 6
        for i in range(1, NSUB + 1):
            ta = lo + (h - lo) * mp.mpf(i - 1) / NSUB
            tb = lo + (h - lo) * mp.mpf(i) / NSUB
            va = prev if i == 1 else poly(Y, ta) + 1
            vb = poly(Y, tb) + 1
            if (va > 0) != (vb > 0):
                # regula-falsi + Newton
                A, B, fa = ta, tb, va
                for _ in range(200):
                    m = (A + B) / 2
                    fm = poly(Y, m) + 1
                    if (fm > 0) == (fa > 0):
                        A, fa = m, fm
                    else:
                        B = m
                    if B - A < abs(B) * mp.mpf(10) ** (-(mp.mp.dps - 4)):
                        break
                tr = (A + B) / 2
                for _ in range(60):
                    f = poly(Y, tr) + 1
                    d = dpoly(Y, tr)
                    if d == 0:
                        break
                    dt = f / d
                    tr -= dt
                    if abs(dt) < abs(tr) * mp.mpf(10) ** (-(mp.mp.dps - 4)):
                        break
                xr = poly(X, tr)
                ok = (xr > 1) if side == "right" else (xr < 1)
                if ok and (t + tr) > mp.mpf("1e-9"):
                    return xr, t + tr
            prev = vb
        x = poly(X, h)
        y = poly(Y, h)
        t += h
        if abs(x) > mp.mpf("1e7") or abs(y) > mp.mpf("1e7"):
            return None, "bound"
        if t > tmax:
            return None, "time"
    return None, "steps"


def D(a, a20, a11, a01, a10, x0, side="right", **kw):
    r = ret(a, a20, a11, a01, a10, x0, side=side, **kw)
    if r[0] is None:
        return None, r[1]
    return r[0] - mp.mpf(x0), r[1]
