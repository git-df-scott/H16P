#!/usr/bin/env python3
"""Generate reproducible Q4 positive, negative, and double-root controls."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.optimize import brentq

from q4_integrals import basis_float, basis_mp, basis_orbit_float, canonical_mu, integral_mp


POSITIVE_KAPPA = 4.0
POSITIVE_ROOTS = (2.0, 3.0, 3.7)


def force_roots(kappa, roots, dps=70):
    """Set mu_4=1 and solve the three linear vanishing equations."""
    with mp.workdps(dps):
        matrix = mp.matrix([basis_mp(kappa, s, dps)[:3] for s in roots])
        rhs = mp.matrix([-basis_mp(kappa, s, dps)[3] for s in roots])
        solution = mp.lu_solve(matrix, rhs)
        return tuple(solution) + (mp.mpf(1),)


def double_root_control(kappa=4.0, s0=2.5):
    """Numerical double-root control with mu_3=0 and mu_4=1."""
    step = 2.0e-5
    b0 = basis_float(kappa, s0, order=192)
    derivative = (
        basis_float(kappa, s0 + step, order=192)
        - basis_float(kappa, s0 - step, order=192)
    ) / (2.0 * step)
    matrix = np.asarray([[b0[0], b0[1]], [derivative[0], derivative[1]]])
    rhs = -np.asarray([b0[3], derivative[3]])
    m1, m2 = np.linalg.solve(matrix, rhs)
    return np.asarray((m1, m2, 0.0, 1.0))


def sign_change_roots(kappa, mu, points=1201, margin=2.0e-4):
    left = 1.0 + margin * (kappa - 1.0)
    right = kappa - margin * (kappa - 1.0)
    grid = np.linspace(left, right, points)
    values = np.asarray([np.dot(mu, basis_float(kappa, s)) for s in grid])
    roots = []
    for a, b, fa, fb in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
        if fa * fb < 0:
            roots.append(
                brentq(lambda x: float(np.dot(mu, basis_float(kappa, x))), a, b)
            )
    return roots


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="q4/data/controls.json")
    args = parser.parse_args()

    positive_mu_mp = force_roots(POSITIVE_KAPPA, POSITIVE_ROOTS)
    positive_mu = np.asarray([float(v) for v in positive_mu_mp])
    roots = sign_change_roots(POSITIVE_KAPPA, positive_mu)
    residuals = [integral_mp(POSITIVE_KAPPA, positive_mu_mp, s, 70) for s in POSITIVE_ROOTS]

    compare_s = 2.5
    primary = np.asarray([float(v) for v in basis_mp(4, compare_s, 60)])
    gauss = basis_float(4, compare_s, order=192)
    orbit = basis_orbit_float(4, compare_s)

    double_mu = double_root_control()
    s0, step = 2.5, 2.0e-4
    fminus = float(np.dot(double_mu, basis_float(4, s0 - step, 192)))
    fzero = float(np.dot(double_mu, basis_float(4, s0, 192)))
    fplus = float(np.dot(double_mu, basis_float(4, s0 + step, 192)))
    curvature = (fplus - 2.0 * fzero + fminus) / step**2

    record = {
        "status": "NONRIGOROUS_NUMERICAL_CONTROLS",
        "positive": {
            "provenance": "constructed in this audit; finite instantiation of Zhao Corollary 9",
            "kappa": POSITIVE_KAPPA,
            "forced_s": list(POSITIVE_ROOTS),
            "mu_decimal_60": [mp.nstr(v, 60) for v in positive_mu_mp],
            "forced_root_residuals_decimal_60": [mp.nstr(v, 60) for v in residuals],
            "isolated_sign_change_roots": roots,
            "root_count": len(roots),
        },
        "negative": {
            "kappa": 4.0,
            "mu": [1, 0, 0, 0],
            "reason": "I=h*I00<0 throughout the open annulus",
            "root_count": 0,
        },
        "degenerate": {
            "provenance": "constructed in this audit",
            "kappa": 4.0,
            "s0": s0,
            "mu": double_mu.tolist(),
            "canonical_mu": canonical_mu(double_mu).tolist(),
            "I_s0": fzero,
            "I_s0_minus_delta": fminus,
            "I_s0_plus_delta": fplus,
            "finite_difference_Iss": curvature,
        },
        "independent_evaluator_check": {
            "kappa": 4.0,
            "s": compare_s,
            "arbitrary_precision_area": primary.tolist(),
            "fixed_gauss_area": gauss.tolist(),
            "hamiltonian_orbit_green": orbit.tolist(),
            "max_abs_gauss_difference": float(np.max(np.abs(primary - gauss))),
            "max_abs_orbit_difference": float(np.max(np.abs(primary - orbit))),
        },
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
