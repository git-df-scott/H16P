"""Validated half returns and displacement enclosures for the seed family.

    P = (b-2)/4 + (1-b) y + a x^2 + b y^2 + e1 x + e2 x y
    Q = e0 - 2 x y

Section {x = 0}, coordinate s = log|y|, side = sign(y).  The displacement

    D(s) = log|y_forward| - log|y_backward|

is the mismatch between the forward and backward half returns; D(s) = 0 means
the two half orbits close up, i.e. the orbit through (0, side e^s) is
periodic.  Everything returned here is an interval enclosing the true value.

The crossing is localised rigorously: on the step whose x-enclosure straddles
zero, the Taylor polynomial in t (with its validated remainder) is used, the
derivative is checked to be of one sign over the whole step so the crossing is
unique, and the crossing time is then bisected with interval evaluations.
"""
import math

from mpmath import iv, mp

from taylor import (QuadraticField, contains, ival, poly_eval, set_precision,
                    step, width)


class Unresolved(RuntimeError):
    """A validated return that could not be completed. Never read as a zero."""


def seed_field(a, b, e0, e1, e2):
    """Coefficient vectors in the basis (1, x, y, x^2, x y, y^2)."""
    p = [(ival(b) - 2) / 4, ival(e1), 1 - ival(b), ival(a), ival(e2), ival(b)]
    q = [ival(e0), ival(0), ival(0), ival(0), ival(-2), ival(0)]
    return QuadraticField(p, q)


def _series_at(coeffs, remainder, t, order):
    return poly_eval(coeffs, t, order, remainder)


def _crossing(fld, X, Y, XB, YB, h, order):
    """Rigorous crossing time in [0, h] where x vanishes, and y there."""
    # x'(t) over the whole step must keep one sign, so the crossing is unique.
    dX = [X[k] * ival(k) for k in range(1, order)]
    dXB = XB[order] * ival(order)
    dcheck = _series_at(dX, dXB, iv.mpf([0, h]), order - 1)
    if dcheck.a <= 0 <= dcheck.b:
        raise Unresolved('x-derivative not of one sign across the step')
    lo, hi = mp.mpf(0), mp.mpf(h)
    x_lo = _series_at(X, XB[order], ival(lo), order)
    for _ in range(200):
        mid = (lo + hi) / 2
        xm = _series_at(X, XB[order], ival(mid), order)
        if xm.a <= 0 <= xm.b:
            break
        if (xm.b < 0) == (x_lo.b < 0):
            lo = mid
            x_lo = xm
        else:
            hi = mid
        if hi - lo < mp.mpf(10) ** (-(mp.dps - 6)):
            break
    T = iv.mpf([lo, hi])
    yT = _series_at(Y, YB[order], T, order)
    xT = _series_at(X, XB[order], T, order)
    if not (xT.a <= 0 <= xT.b):
        raise Unresolved('bracketed crossing does not enclose x = 0')
    return T, yT


def half_return(a, b, e0, e1, e2, s0, side, direction, h=0.004, order=16,
                tmax=40.0):
    """Validated half return; returns an interval enclosing log|y| at the
    next crossing of {x = 0}, together with the enclosed crossing time."""
    fld = seed_field(a, b, e0, e1, e2)
    if direction < 0:                      # integrate backwards in time
        fld = QuadraticField([-c for c in fld.p], [-c for c in fld.q])
    s_iv = s0 if hasattr(s0, 'a') else ival(s0)
    y0 = iv.exp(s_iv) if side > 0 else -iv.exp(s_iv)
    x = ival(0)
    y = y0
    fx0, _ = fld.eval(x, y)
    if fx0.a <= 0 <= fx0.b:
        raise Unresolved('launch not transverse')
    launch_positive = fx0.a > 0
    t = mp.mpf(0)
    left = False
    h_cur = mp.mpf(h)
    nsteps = 0
    while t < tmax:
        nsteps += 1
        r = None
        h_try = h_cur
        h_floor = mp.mpf(h) * mp.mpf(2) ** -20
        for _ in range(40):                 # adaptive: shrink until validated
            r = step(fld, x, y, h_try, order)
            if r is not None:
                # accept a validated step whose widths have not blown up
                grow = max(width(r[0]) - width(x), width(r[1]) - width(y))
                if grow < max(width(x), width(y), mp.mpf(1)) * mp.mpf('0.05'):
                    break
                r = None
            r = None
            h_try = h_try / 2
        if r is None:
            raise Unresolved('no validated step at t = %s (h down to %s)'
                             % (t, h_try))
        h_cur = min(mp.mpf(h), h_try * 2)
        h = h  # keep the nominal step for the crossing localisation below
        x1, y1, B, (X, Y, XB, YB) = r
        if not left:
            # wait until the orbit is clearly off the section
            if (x1.a > 0) if launch_positive else (x1.b < 0):
                left = True
            x, y, t = x1, y1, t + h_try
            continue
        crosses = (x1.b < 0) if launch_positive else (x1.a > 0)
        if crosses:
            T, yT = _crossing(fld, X, Y, XB, YB, h_try, order)
            if yT.a <= 0 <= yT.b:
                raise Unresolved('y encloses zero at the crossing')
            ay = yT if yT.a > 0 else iv.mpf([-yT.b, -yT.a])
            return {'s_end': iv.log(ay), 'time': ival(t) + T,
                    'y_end': yT, 'steps': nsteps}
        x, y, t = x1, y1, t + h_try
    raise Unresolved('no crossing within t = %g' % tmax)


def displacement(a, b, e0, e1, e2, s0, side, **kw):
    """Interval enclosing D(s0). Its sign is rigorous when it excludes 0."""
    fwd = half_return(a, b, e0, e1, e2, s0, side, +1, **kw)
    bwd = half_return(a, b, e0, e1, e2, s0, side, -1, **kw)
    D = fwd['s_end'] - bwd['s_end']
    return {'D': D, 'width': width(D), 'forward': fwd, 'backward': bwd,
            'sign': ('negative' if D.b < 0 else
                     'positive' if D.a > 0 else 'indeterminate')}
