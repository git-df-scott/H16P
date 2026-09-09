"""Validated Taylor integration for planar quadratic vector fields.

    x' = p0 + p1 x + p2 y + p3 x^2 + p4 x y + p5 y^2
    y' = q0 + q1 x + q2 y + q3 x^2 + q4 x y + q5 y^2

Everything is interval arithmetic (mpmath.iv, outward rounded, arbitrary
precision).  A step is validated in the standard two-phase way:

  1. Rough enclosure.  Find a box B with z0 + [0,h] f(B) subset of B.  By
     Picard-Lindelof this proves the solution exists on [0,h] and stays in B.
  2. Taylor step with remainder.  Expand to order K from the current
     enclosure and bound the order-(K+1) term by running the same
     coefficient recurrence started from B.

The result encloses the true solution.  There is no floating-point
"approximately" anywhere in the returned quantities.
"""
from mpmath import iv, mp


def set_precision(dps):
    iv.dps = dps
    mp.dps = dps + 10


def ival(lo, hi=None):
    return iv.mpf([lo, lo if hi is None else hi])


def sqr(u):
    """Tight interval square: correct for intervals straddling zero."""
    a, b = u.a, u.b
    if a >= 0:
        return iv.mpf([a * a, b * b])
    if b <= 0:
        return iv.mpf([b * b, a * a])
    m = max(-a, b)
    return iv.mpf([0, m * m])


def contains(outer, inner):
    return outer.a <= inner.a and inner.b <= outer.b


def width(u):
    return u.b - u.a


class QuadraticField:
    """Coefficient vectors p, q of length 6: (1, x, y, x^2, x y, y^2)."""

    def __init__(self, p, q):
        self.p = [ival(c) if not hasattr(c, 'a') else c for c in p]
        self.q = [ival(c) if not hasattr(c, 'a') else c for c in q]

    def eval(self, x, y):
        p, q = self.p, self.q
        x2, xy, y2 = sqr(x), x * y, sqr(y)
        return (p[0] + p[1] * x + p[2] * y + p[3] * x2 + p[4] * xy + p[5] * y2,
                q[0] + q[1] * x + q[2] * y + q[3] * x2 + q[4] * xy + q[5] * y2)

    def taylor(self, x0, y0, order):
        """Taylor coefficients of the solution through (x0, y0), as intervals."""
        p, q = self.p, self.q
        X, Y = [x0], [y0]
        for k in range(order):
            x2 = sqr(X[0]) if k == 0 else None
            xy = None
            y2 = sqr(Y[0]) if k == 0 else None
            if k:
                x2 = ival(0)
                y2 = ival(0)
                half = k // 2
                for i in range(half + 1):
                    if 2 * i == k:
                        x2 = x2 + sqr(X[i])
                        y2 = y2 + sqr(Y[i])
                    else:
                        x2 = x2 + ival(2) * X[i] * X[k - i]
                        y2 = y2 + ival(2) * Y[i] * Y[k - i]
            xy = ival(0)
            for i in range(k + 1):
                xy = xy + X[i] * Y[k - i]
            px = p[1] * X[k] + p[2] * Y[k] + p[3] * x2 + p[4] * xy + p[5] * y2
            qy = q[1] * X[k] + q[2] * Y[k] + q[3] * x2 + q[4] * xy + q[5] * y2
            if k == 0:
                px = px + p[0]
                qy = qy + q[0]
            X.append(px / ival(k + 1))
            Y.append(qy / ival(k + 1))
        return X, Y


def rough_enclosure(fld, x0, y0, h, tries=60):
    """Box B with z0 + [0,h] f(B) subset of B; None if h is too large.

    Inclusion certifies, by Picard-Lindelof, that the solution exists on
    [0,h] and does not leave B.
    """
    hh = iv.mpf([0, h])
    fx, fy = fld.eval(x0, y0)
    mag = max(abs(fx.a), abs(fx.b), abs(fy.a), abs(fy.b), mp.mpf(1))
    r = mp.mpf(h) * mag + mp.mpf(h) * mp.mpf(h)
    for _ in range(tries):
        Bx = iv.mpf([x0.a - r, x0.b + r])
        By = iv.mpf([y0.a - r, y0.b + r])
        gx, gy = fld.eval(Bx, By)
        Nx, Ny = x0 + hh * gx, y0 + hh * gy
        if contains(Bx, Nx) and contains(By, Ny):
            for _ in range(3):
                gx, gy = fld.eval(Bx, By)
                Nx, Ny = x0 + hh * gx, y0 + hh * gy
                if contains(Bx, Nx) and contains(By, Ny):
                    Bx, By = Nx, Ny
                else:
                    break
            return Bx, By
        r = r * 2
    return None


def step(fld, x0, y0, h, order):
    """One validated step of length h.

    Returns (x1, y1, rough enclosure B, (X, Y, XB, YB)) where X, Y are the
    Taylor coefficients from the current point and XB, YB the coefficients
    over B, whose order-K term bounds the remainder.
    """
    enc = rough_enclosure(fld, x0, y0, h)
    if enc is None:
        return None
    Bx, By = enc
    X, Y = fld.taylor(x0, y0, order)
    XB, YB = fld.taylor(Bx, By, order)
    hi = ival(h)
    pw = [ival(1)]
    for _ in range(order):
        pw.append(pw[-1] * hi)
    x1 = X[0]
    y1 = Y[0]
    for k in range(1, order):
        x1 = x1 + X[k] * pw[k]
        y1 = y1 + Y[k] * pw[k]
    x1 = x1 + XB[order] * pw[order]
    y1 = y1 + YB[order] * pw[order]
    return x1, y1, (Bx, By), (X, Y, XB, YB)


def poly_eval(coeffs, t, order, remainder):
    """Sum_{k<order} c_k t^k + remainder t^order, all intervals."""
    out = coeffs[0]
    tp = ival(1)
    for k in range(1, order):
        tp = tp * t
        out = out + coeffs[k] * tp
    return out + remainder * (tp * t)
