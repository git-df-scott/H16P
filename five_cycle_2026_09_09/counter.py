"""Nest-aware limit-cycle counter for planar quadratic fields.

Chart (origin is an antisaddle):

    x' = lam*x - y + l*x^2 + m*x*y + n*y^2
    y' = x + lam*y + a*x^2 + b*x*y + c*y^2

For a focus at `centre` we take the ray {y = cy, x > cx}, integrate one full
revolution, and record D(r) = log(R(r)/r).  Sign changes of D count limit
cycles in that nest.  The crossing direction is taken from the rotation
sense of the linearization, not assumed.

A radius whose return does not resolve ends the domain and is reported as a
missing domain -- never read as a zero.

Evidence class: NUM.  Floating point, no interval arithmetic.
"""
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ESCAPE = 1e3
TMAX = 200.0


def rhs_factory(coeffs):
    lam, l, m, n, a, b, c = coeffs

    def rhs(t, z):
        x, y = z
        if not (abs(x) < ESCAPE and abs(y) < ESCAPE):
            return (0.0, 0.0)
        return (lam * x - y + l * x * x + m * x * y + n * y * y,
                x + lam * y + a * x * x + b * x * y + c * y * y)
    return rhs


def equilibria(coeffs):
    """Real equilibria of the two conics, via the Bezout resultant in y."""
    lam, l, m, n, a, b, c = coeffs
    A = np.array([n]) if n else np.array([0.0])
    B = np.array([m, -1.0])
    C = np.array([l, lam, 0.0])
    Dp = np.array([c]) if c else np.array([0.0])
    E = np.array([b, lam])
    F = np.array([a, 1.0, 0.0])
    mul, sub = np.polymul, np.polysub
    R = sub(mul(sub(mul(A, F), mul(C, Dp)), sub(mul(A, F), mul(C, Dp))),
            mul(sub(mul(A, E), mul(B, Dp)), sub(mul(B, F), mul(C, E))))
    R = np.trim_zeros(np.asarray(R, dtype=float), 'f')
    pts = []
    if len(R) < 2:
        return pts
    for root in np.roots(R):
        if abs(root.imag) > 1e-9 * max(1.0, abs(root.real)):
            continue
        xv = float(root.real)
        for P2, P1, P0 in ((np.polyval(A, xv), np.polyval(B, xv), np.polyval(C, xv)),
                           (np.polyval(Dp, xv), np.polyval(E, xv), np.polyval(F, xv))):
            cand = []
            if abs(P2) > 1e-12:
                disc = P1 * P1 - 4 * P2 * P0
                if disc >= 0:
                    s = math.sqrt(disc)
                    cand = [(-P1 + s) / (2 * P2), (-P1 - s) / (2 * P2)]
            elif abs(P1) > 1e-12:
                cand = [-P0 / P1]
            for yv in cand:
                r1 = lam * xv - yv + l * xv**2 + m * xv * yv + n * yv**2
                r2 = xv + lam * yv + a * xv**2 + b * xv * yv + c * yv**2
                scale = 1.0 + xv * xv + yv * yv
                if max(abs(r1), abs(r2)) < 1e-7 * scale and not any(
                        abs(xv - u) < 1e-7 and abs(yv - v) < 1e-7 for u, v in pts):
                    pts.append((xv, yv))
            if cand:
                break
    return pts


def jacobian(coeffs, p):
    lam, l, m, n, a, b, c = coeffs
    x, y = p
    return np.array([[lam + 2 * l * x + m * y, -1.0 + m * x + 2 * n * y],
                     [1.0 + 2 * a * x + b * y, lam + b * x + 2 * c * y]])


def foci(coeffs):
    """Antisaddles with non-real eigenvalues, each with its rotation sense."""
    out = []
    for p in equilibria(coeffs):
        J = jacobian(coeffs, p)
        if np.linalg.det(J) <= 0:
            continue
        if (J[0, 0] + J[1, 1]) ** 2 - 4 * np.linalg.det(J) >= 0:
            continue                      # node, not a focus
        out.append({'centre': p, 'sense': 1.0 if J[1, 0] > 0 else -1.0,
                    'trace': float(J[0, 0] + J[1, 1])})
    return out


def displacement(coeffs, centre, sense, r, rtol=1e-11):
    """log(R(r)/r) after one revolution; None if the return does not resolve."""
    rhs = rhs_factory(coeffs)
    cx, cy = centre
    lead = solve_ivp(rhs, [0.0, 1e-4], [cx + r, cy], method='DOP853',
                     rtol=rtol, atol=rtol * 1e-3)
    if not lead.success:
        return None
    ev = lambda t, z: z[1] - cy
    ev.terminal = False
    ev.direction = sense
    sol = solve_ivp(rhs, [1e-4, TMAX], lead.y[:, -1], method='DOP853',
                    events=ev, rtol=rtol, atol=rtol * 1e-3)
    if not sol.success:
        return None
    for z in sol.y_events[0]:
        if z[0] - cx > 0:
            R = z[0] - cx
            return math.log(R / r) if R > 0 and np.isfinite(R) else None
    return None


def scan(coeffs, focus, r_lo, r_hi, n=40, rtol=1e-11):
    """D on a logarithmic radial grid; stops at the first unresolved return."""
    centre, sense = focus['centre'], focus['sense']
    grid = np.exp(np.linspace(math.log(r_lo), math.log(r_hi), n))
    vals = []
    for r in grid:
        d = displacement(coeffs, centre, sense, float(r), rtol)
        if d is None:
            break
        vals.append((float(r), d))
    brackets = [(vals[i][0], vals[i + 1][0]) for i in range(len(vals) - 1)
                if vals[i][1] * vals[i + 1][1] < 0]
    return {'evaluated': len(vals), 'resolved_all': len(vals) == n,
            'domain_end': vals[-1][0] if vals else None,
            'sign_changes': len(brackets), 'brackets': brackets, 'values': vals}


def refine(coeffs, focus, bracket, rtol=1e-12):
    f = lambda r: displacement(coeffs, focus['centre'], focus['sense'], r, rtol)
    lo, hi = bracket
    if f(lo) is None or f(hi) is None:
        return None
    return brentq(f, lo, hi, xtol=1e-15, rtol=8.9e-16, maxiter=200)
