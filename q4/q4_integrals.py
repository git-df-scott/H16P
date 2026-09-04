#!/usr/bin/env python3
"""Numerical Q4 Abelian-integral primitives.

The primary evaluator uses arbitrary-precision one-dimensional area quadrature.
The independent control evaluator transports Green-theorem line integrals around
the Hamiltonian orbit with scipy.  Neither routine is interval-rigorous.
"""

from __future__ import annotations

from functools import lru_cache
from math import log, pi, sqrt

import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


def h_of_s(kappa, s):
    return -mp.mpf(2) / 3 * mp.sqrt(mp.mpf(s) / mp.mpf(kappa))


def canonical_mu(mu):
    """Canonical representative of [mu] in RP^3 (infinity norm one)."""
    values = np.asarray(mu, dtype=float)
    scale = float(np.max(np.abs(values)))
    if not scale:
        raise ValueError("mu=0 is the identically-zero degenerate case")
    values = values / scale
    first = int(np.flatnonzero(np.abs(values) > 1.0e-15)[0])
    if values[first] < 0:
        values = -values
    return values


def q4_coefficients(rho: float, sign: int = 1):
    """Real coefficients of the original complex Q4 normal form.

    alpha=b+i*c is rational whenever rho is rational:
      b=2(1-rho^2)/(1+rho^2), c=sign*4rho/(1+rho^2).
    """
    if rho <= 0 or sign not in (-1, 1):
        raise ValueError("require rho>0 and sign in {-1,1}")
    den = 1.0 + rho * rho
    b = 2.0 * (1.0 - rho * rho) / den
    c = sign * 4.0 * rho / den
    return {
        "b": b,
        "c": c,
        "kappa": 1.0 + rho * rho,
        "dx": [0.0, 0.0, 1.0, 6.0 + b, 2.0 * c, -(2.0 + b)],
        "dy": [0.0, -1.0, 0.0, c, 8.0 - 2.0 * b, -c],
        "monomial_order": ["1", "x", "y", "x^2", "xy", "y^2"],
    }


def _bisect_mp(fun, left, right, iterations=220):
    left, right = mp.mpf(left), mp.mpf(right)
    fleft, fright = fun(left), fun(right)
    if fleft == 0:
        return left
    if fright == 0:
        return right
    if fleft * fright > 0:
        raise ValueError("root is not bracketed")
    for _ in range(iterations):
        mid = (left + right) / 2
        fmid = fun(mid)
        if fleft * fmid <= 0:
            right, fright = mid, fmid
        else:
            left, fleft = mid, fmid
    return (left + right) / 2


def _vertical_value_mp(kappa, h, x):
    y = mp.sqrt(((kappa - 1) * x * x + 1) / kappa)
    return (
        mp.mpf(2) / 3 * (kappa - 1) * x**3
        - (kappa - 1) * x * x * y
        + kappa / 3 * y**3
        - y
        - h
    )


def _x_bounds_mp(kappa, h):
    fun = lambda x: _vertical_value_mp(kappa, h, x)
    left = _bisect_mp(fun, 0, 1)
    right_edge = mp.mpf(2)
    while fun(right_edge) < 0:
        right_edge *= 2
    right = _bisect_mp(fun, 1, right_edge)
    return left, right


def _y_pair_mp(kappa, h, x):
    p = -3 * ((kappa - 1) * x * x + 1) / kappa
    q = (2 * (kappa - 1) * x**3 - 3 * h) / kappa
    arg = (3 * q / (2 * p)) * mp.sqrt(-3 / p)
    arg = max(mp.mpf(-1), min(mp.mpf(1), arg))
    theta = mp.acos(arg) / 3
    radius = 2 * mp.sqrt(-p / 3)
    roots = sorted(radius * mp.cos(theta - 2 * mp.pi * j / 3) for j in range(3))
    return roots[1], roots[2]


def basis_mp(kappa, s, dps=60):
    """Return Zhao's four basis functions at (kappa,s)."""
    with mp.workdps(dps):
        kappa, s = mp.mpf(kappa), mp.mpf(s)
        if not (kappa > 1 and 1 < s < kappa):
            raise ValueError("require kappa>1 and 1<s<kappa")
        h = h_of_s(kappa, s)
        xlo, xhi = _x_bounds_mp(kappa, h)
        mid, half = (xlo + xhi) / 2, (xhi - xlo) / 2

        def integrate(fun):
            def transformed(theta):
                x = mid + half * mp.sin(theta)
                lower, upper = _y_pair_mp(kappa, h, x)
                return fun(x, lower, upper) * half * mp.cos(theta)

            return mp.quad(transformed, [-mp.pi / 2, 0, mp.pi / 2])

        i00 = integrate(lambda _x, lo, hi: hi - lo)
        i10 = integrate(lambda x, lo, hi: x * (hi - lo))
        i01 = integrate(lambda _x, lo, hi: (hi * hi - lo * lo) / 2)
        im10 = integrate(lambda x, lo, hi: (hi - lo) / x)
        im11 = integrate(lambda x, lo, hi: (hi * hi - lo * lo) / (2 * x))
        return tuple(
            +v
            for v in (
                h * i00,
                i10,
                i01,
                2 * im10 + 3 * kappa * h * im11,
            )
        )


@lru_cache(maxsize=16)
def _legendre_rule(order):
    return leggauss(order)


def _vertical_value_float(kappa, h, x):
    y = sqrt(((kappa - 1.0) * x * x + 1.0) / kappa)
    return (
        2.0 / 3.0 * (kappa - 1.0) * x**3
        - (kappa - 1.0) * x * x * y
        + kappa / 3.0 * y**3
        - y
        - h
    )


def _x_bounds_float(kappa, h):
    fun = lambda x: _vertical_value_float(kappa, h, x)
    left = brentq(fun, 0.0, 1.0, xtol=1.0e-14)
    right_edge = 2.0
    while fun(right_edge) < 0:
        right_edge *= 2.0
    right = brentq(fun, 1.0, right_edge, xtol=1.0e-14)
    return left, right


def basis_float(kappa, s, order=96):
    """Fast fixed-Gauss version of basis_mp for screening."""
    kappa, s = float(kappa), float(s)
    if not (kappa > 1.0 and 1.0 < s < kappa):
        raise ValueError("require kappa>1 and 1<s<kappa")
    h = -2.0 / 3.0 * sqrt(s / kappa)
    xlo, xhi = _x_bounds_float(kappa, h)
    nodes, weights = _legendre_rule(order)
    theta = (pi / 2.0) * nodes
    x = (xlo + xhi) / 2.0 + (xhi - xlo) / 2.0 * np.sin(theta)
    jac = (xhi - xlo) / 2.0 * np.cos(theta) * (pi / 2.0)
    p = -3.0 * ((kappa - 1.0) * x * x + 1.0) / kappa
    q = (2.0 * (kappa - 1.0) * x**3 - 3.0 * h) / kappa
    arg = np.clip((3.0 * q / (2.0 * p)) * np.sqrt(-3.0 / p), -1.0, 1.0)
    angle = np.arccos(arg) / 3.0
    radius = 2.0 * np.sqrt(-p / 3.0)
    roots = np.sort(
        np.stack([radius * np.cos(angle - 2.0 * pi * j / 3.0) for j in range(3)]),
        axis=0,
    )
    lower, upper = roots[1], roots[2]

    def integral(values):
        return float(np.dot(weights, values * jac))

    i00 = integral(upper - lower)
    i10 = integral(x * (upper - lower))
    i01 = integral((upper**2 - lower**2) / 2.0)
    im10 = integral((upper - lower) / x)
    im11 = integral((upper**2 - lower**2) / (2.0 * x))
    return np.asarray((h * i00, i10, i01, 2.0 * im10 + 3.0 * kappa * h * im11))


def integral_mp(kappa, mu, s, dps=60):
    with mp.workdps(dps):
        coefficients = tuple(mp.mpf(str(v)) for v in mu)
        return +mp.fdot(coefficients, basis_mp(kappa, s, dps=dps))


def basis_orbit_float(kappa, s):
    """Independent ODE/Green-theorem evaluation of the four basis functions."""
    kappa, s = float(kappa), float(s)
    h = -2.0 / 3.0 * sqrt(s / kappa)

    def hamiltonian(x, y):
        return (
            2.0 / 3.0 * (kappa - 1.0) * x**3
            - (kappa - 1.0) * x * x * y
            + kappa / 3.0 * y**3
            - y
        )

    section = lambda x: hamiltonian(x, 1.0) - h
    edge = 2.0
    while section(edge) < 0:
        edge *= 2.0
    x0 = brentq(section, 1.0, edge)

    def rhs(_t, state):
        x, y = state[:2]
        dx = -1.0 - (kappa - 1.0) * x * x + kappa * y * y
        dy = -2.0 * (kappa - 1.0) * x * (x - y)
        qs = (x, x * x / 2.0, x * y, log(x), log(x) * y)
        return np.asarray((dx, dy, *(-q * dy for q in qs)))

    initial = np.asarray((x0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0))
    dt = 1.0e-10
    nudged = initial + dt * rhs(0.0, initial)

    def event(_t, state):
        return state[1] - 1.0

    event.direction = -1
    event.terminal = True
    sol = solve_ivp(
        rhs,
        (dt, 1000.0),
        nudged,
        events=event,
        method="DOP853",
        rtol=2.0e-12,
        atol=2.0e-13,
        max_step=0.01,
    )
    if not len(sol.t_events[0]):
        raise RuntimeError("Hamiltonian orbit did not return")
    values = sol.y_events[0][0]
    i00, i10, i01, im10, im11 = values[2:]
    return np.asarray((h * i00, i10, i01, 2.0 * im10 + 3.0 * kappa * h * im11))


def alpha_beta_from_mu(kappa, mu):
    """Zhao's (alpha_0,alpha_1,alpha_2,beta_0,beta_1).

    ``mu`` is in the Abelian-integral basis returned by :func:`basis_float`.
    Zhao relabels two coefficients between equations (20) and (21); applying
    that change here is essential before using the pruning theorem.
    """
    k, (m1, old_m2, old_m3, m4) = float(kappa), tuple(map(float, mu))
    m2 = -2.0 * old_m2 / 3.0 - 2.0 * (k - 1.0) * old_m3 / (3.0 * k)
    m3 = -2.0 * old_m3 / (3.0 * k)
    a1 = (
        128 * (-117 - 265 * k + 30 * k**2) / (27 * k**2) * m1
        + 16 * (-174 + 5 * k) / (3 * k) * m2
        - 16 * (243 + 121 * k) / (3 * k) * m3
        + 64 * (k - 1) * (-119 + 60 * k) / (9 * k) * m4
    )
    a2 = 1088 * (21 + k) / (27 * k**2) * m1 + 544 / k * m2 + 1088 / k * m3 + 1088 * (k - 1) / (9 * k) * m4
    b0 = (
        -256 * (k - 1) / (3 * k) * m1
        - 32 * (-27 + 25 * k + 15 * k**2) / (3 * k) * m2
        - 16 * (k - 1) * (54 + 31 * k) / (3 * k) * m3
        - 64 * (k - 1) * (-18 + 5 * k) / (3 * k) * m4
    )
    b1 = (
        -64 * (k - 1) * (18 + 77 * k) / (27 * k**2) * m1
        - 16 * (-111 + 137 * k) / (3 * k) * m2
        - 768 * (k - 1) / k * m3
        - 64 * (k - 1) * (-116 + 77 * k) / (9 * k) * m4
    )
    a0 = b0 - k * b1 - k * a1 - k * k * a2
    return np.asarray((a0, a1, a2, b0, b1))
