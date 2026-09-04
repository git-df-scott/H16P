#!/usr/bin/env python3
"""Bounded, nonrigorous Q4 candidate screen.

This is a ranking tool, never a proof.  It precomputes the four Abelian basis
functions, applies Zhao's necessary five-zero filters, and counts robust sign
changes.  The ``astra`` mode has an explicit CPU-time fuse.
"""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from q4_integrals import alpha_beta_from_mu, basis_float, canonical_mu


KAPPA_MAX = 85.0 / 23.0


def zhao_five_filter(kappa, mu, tolerance=1.0e-12):
    a0, a1, a2, b0, b1 = alpha_beta_from_mu(kappa, mu)
    if abs(b1) <= tolerance:
        return False, "beta1_zero"
    a0, a1, a2, beta0 = a0 / b1, a1 / b1, a2 / b1, b0 / b1
    lower = (23.0 * kappa - 54.0) / 31.0
    if not lower < beta0 < 1.0:
        return False, "beta0_strip"
    p_beta = a2 * beta0**2 + a1 * beta0 + a0
    threshold = 25.0 * (1.0 - beta0) / (432.0 * (kappa - 1.0) ** 2)
    if p_beta < threshold:
        return False, "zhao_comment_bound"
    return True, "survives"


def normalized_grid(kappa, points):
    # Sine clustering samples both singular endpoints without touching either.
    theta = np.linspace(-np.pi / 2.0, np.pi / 2.0, points + 2)[1:-1]
    r = (1.0 + np.sin(theta)) / 2.0
    return 1.0 + (kappa - 1.0) * r


def basis_table(kappa, points, order):
    grid = normalized_grid(kappa, points)
    table = np.asarray([basis_float(kappa, s, order=order) for s in grid])
    return grid, table


def robust_crossings(values, relative_floor=1.0e-10):
    scale = max(float(np.max(np.abs(values))), np.finfo(float).tiny)
    signs = np.sign(np.where(np.abs(values) >= relative_floor * scale, values, 0.0))
    signs = signs[signs != 0]
    return int(np.count_nonzero(signs[:-1] * signs[1:] < 0)) if len(signs) > 1 else 0


def triple_zero_direction(table, rng):
    """Projective mu forced through three separated sampled zeros by SVD."""
    n = len(table)
    indices = np.sort(rng.choice(np.arange(3, n - 3), size=3, replace=False))
    if np.min(np.diff(indices)) < 3:
        return None
    _u, sigma, vh = np.linalg.svd(table[indices], full_matrices=True)
    relative_singular_value = float(sigma[-1] / sigma[0])
    return canonical_mu(vh[-1]), indices, relative_singular_value


def refine_roots(kappa, mu, grid, values):
    roots = []
    slopes = []
    for a, b, fa, fb in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
        if fa * fb >= 0:
            continue
        fun = lambda s: float(np.dot(mu, basis_float(kappa, s, order=96)))
        refined_fa, refined_fb = fun(a), fun(b)
        if refined_fa * refined_fb >= 0:
            continue
        root = brentq(fun, a, b, xtol=2e-13)
        step = min(1e-5 * (kappa - 1), (root - 1) / 4, (kappa - root) / 4)
        slope = (fun(root + step) - fun(root - step)) / (2 * step)
        roots.append(root)
        slopes.append(slope)
    return roots, slopes


def track_parameter_path(path, grid_points=257, quadrature_order=96):
    """Track sign-change roots along an explicit ``[(kappa, mu), ...]`` path.

    This cheap predictor interface matches each new root to the nearest prior
    root in normalized r coordinates.  It is intentionally nonrigorous.
    """
    previous = []
    records = []
    for kappa, mu in path:
        mu = canonical_mu(mu)
        grid, table = basis_table(float(kappa), grid_points, quadrature_order)
        values = table @ mu
        roots, slopes = refine_roots(float(kappa), mu, grid, values)
        normalized = [(root - 1.0) / (float(kappa) - 1.0) for root in roots]
        matches = []
        unused = set(range(len(previous)))
        for root in normalized:
            if not unused:
                matches.append(None)
                continue
            index = min(unused, key=lambda i: abs(previous[i] - root))
            unused.remove(index)
            matches.append(index)
        records.append({"kappa": float(kappa), "mu": mu.tolist(), "roots": roots,
                        "normalized_roots": normalized, "slopes": slopes,
                        "matched_previous_indices": matches})
        previous = normalized
    return records


def run_search(args):
    if args.mode == "smoke" and args.cpu_hours > 0.1:
        raise SystemExit("smoke mode fuse: --cpu-hours must be <= 0.1")
    if not (0.0 < args.cpu_hours <= args.max_cpu_hours):
        raise SystemExit("invalid CPU budget")
    if args.kappa_min <= 1 or args.kappa_max >= KAPPA_MAX:
        raise SystemExit(f"require 1 < kappa_min < kappa_max < 85/23={KAPPA_MAX}")

    rng = np.random.default_rng(args.seed)
    started_wall = time.monotonic()
    started_cpu = time.process_time()
    deadline = started_cpu + 3600.0 * args.cpu_hours
    rejection = {"beta1_zero": 0, "beta0_strip": 0, "zhao_comment_bound": 0}
    distribution = {str(i): 0 for i in range(6)}
    leaders = []
    evaluated = surviving = 0
    kappas_done = []
    five_zero_lead = False

    kappas = np.linspace(args.kappa_min, args.kappa_max, args.kappa_count)
    for kappa in kappas:
        if time.process_time() >= deadline:
            break
        grid, table = basis_table(float(kappa), args.grid_points, args.quad_order)
        kappas_done.append(float(kappa))
        for _ in range(args.samples_per_kappa):
            if time.process_time() >= deadline:
                break
            if args.candidate_mode == "triple":
                direction = triple_zero_direction(table, rng)
                if direction is None:
                    continue
                mu, forced_indices, sv_ratio = direction
            else:
                mu = canonical_mu(rng.normal(size=4))
                forced_indices, sv_ratio = None, None
            evaluated += 1
            keep, reason = zhao_five_filter(float(kappa), mu)
            if not keep:
                rejection[reason] += 1
                continue
            surviving += 1
            values = table @ mu
            count = min(5, robust_crossings(values))
            distribution[str(count)] += 1
            scale = max(np.max(np.abs(values)), 1e-300)
            away = np.abs(values) > 1e-10 * scale
            condition = float(np.min(np.abs(values[away])) / scale) if np.any(away) else 0.0
            leader = {
                    "crossings": count,
                    "minimum_relative_sample_value": condition,
                    "kappa": float(kappa),
                    "mu": mu.tolist(),
                }
            if forced_indices is not None:
                leader["forced_grid_indices"] = forced_indices.tolist()
                leader["constraint_matrix_relative_singular_value"] = sv_ratio
            if count >= 3:
                roots, slopes = refine_roots(float(kappa), mu, grid, values)
                leader["refined_roots"] = roots
                leader["root_slopes"] = slopes
                if len(roots) >= 5:
                    five_zero_lead = True
            leaders.append(leader)
            leaders.sort(key=lambda x: (x["crossings"], -x["minimum_relative_sample_value"]), reverse=True)
            del leaders[20:]
            if five_zero_lead:
                break
        if five_zero_lead:
            break

    cpu = time.process_time() - started_cpu
    wall = time.monotonic() - started_wall
    return {
        "status": "NONRIGOROUS_SCREEN_ONLY",
        "mode": args.mode,
        "seed": args.seed,
        "requested_cpu_hours": args.cpu_hours,
        "actual_cpu_seconds": cpu,
        "wall_seconds": wall,
        "host_cpu_count": os.cpu_count(),
        "parameters": {
            "kappa_interval": [args.kappa_min, args.kappa_max],
            "kappa_count": args.kappa_count,
            "samples_per_kappa": args.samples_per_kappa,
            "grid_points": args.grid_points,
            "quadrature_order": args.quad_order,
            "candidate_mode": args.candidate_mode,
        },
        "kappas_completed": kappas_done,
        "projective_points_evaluated": evaluated,
        "zhao_filter_survivors": surviving,
        "rejections": rejection,
        "root_count_distribution_among_survivors": distribution,
        "leaders": leaders,
        "stopped_on_five_zero_lead": five_zero_lead,
        "warning": "sign changes are leads, not certified Abelian zeros or limit cycles",
    }


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("smoke", "astra"), default="smoke")
    parser.add_argument("--cpu-hours", type=float, default=0.02)
    parser.add_argument("--max-cpu-hours", type=float, default=24.0)
    parser.add_argument("--seed", type=int, default=160926)
    parser.add_argument("--kappa-min", type=float, default=1.05)
    parser.add_argument("--kappa-max", type=float, default=3.68)
    parser.add_argument("--kappa-count", type=int, default=9)
    parser.add_argument("--samples-per-kappa", type=int, default=256)
    parser.add_argument("--grid-points", type=int, default=121)
    parser.add_argument("--quad-order", type=int, default=64)
    parser.add_argument("--candidate-mode", choices=("triple", "random"), default="triple")
    parser.add_argument("--output", default="q4/data/smoke.json")
    return parser.parse_args()


def main():
    args = parse_args()
    record = run_search(args)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
