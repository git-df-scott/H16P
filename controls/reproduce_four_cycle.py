#!/usr/bin/env python3
"""Non-rigorous four-cycle regression for the Kuznetsov--Kuznetsova--Leonov field.

This is a discovery/control calculation, not a computer-assisted proof.  It
finds fixed points of a downward Poincare return map on y=0 with scipy's
DOP853 integrator, evaluates planar Floquet multipliers through the divergence
integral, and writes a JSON record and a phase portrait.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


BRACKETS = ((0.428, 0.785), (1.438, 2.637), (8.858, 16.238), (-3945.0, -2893.0))


def field(_t: float, z: np.ndarray) -> tuple[float, float]:
    x, y = z
    return (
        y + x * x + x * y,
        -10.0 * x * x + 2.2 * x * y + 0.7 * y * y - 72.7778 * x + 0.0015 * y,
    )


def downward_section(_t: float, z: np.ndarray) -> float:
    return float(z[1])


downward_section.direction = -1
downward_section.terminal = True


def return_data(x0: float) -> tuple[float, float]:
    # Move a tiny positive time away from the initial point on the section so
    # solve_ivp does not report t=0 as the return.
    dt = 1.0e-11 if abs(x0) > 100.0 else 1.0e-8
    z0 = np.asarray((x0, 0.0), dtype=float)
    z0 = z0 + dt * np.asarray(field(0.0, z0))
    sol = solve_ivp(
        field,
        (dt, 1000.0),
        z0,
        events=downward_section,
        method="DOP853",
        rtol=2.0e-11,
        atol=1.0e-12,
        max_step=0.02,
    )
    if not len(sol.t_events[0]):
        raise RuntimeError(f"no downward return from x={x0}: {sol.message}")
    return float(sol.y_events[0][0, 0]), float(sol.t_events[0][0])


def displacement(x0: float) -> float:
    return return_data(x0)[0] - x0


def orbit_and_multiplier(x0: float, period: float, samples: int = 2000):
    # For a planar autonomous flow, the nontrivial Floquet multiplier is
    # exp(integral div(F) dt) around the periodic orbit.
    def augmented(_t: float, z: np.ndarray) -> tuple[float, float, float]:
        x, y, _w = z
        dx, dy = field(_t, np.asarray((x, y)))
        divergence = 4.2 * x + 2.4 * y + 0.0015
        return dx, dy, divergence

    sol = solve_ivp(
        augmented,
        (0.0, period),
        (x0, 0.0, 0.0),
        method="DOP853",
        rtol=1.0e-12,
        atol=1.0e-13,
        max_step=0.002,
        dense_output=True,
    )
    ts = np.linspace(0.0, period, samples)
    vals = sol.sol(ts)
    multiplier = float(np.exp(sol.y[2, -1]))
    closure_error = float(np.linalg.norm(sol.y[:2, -1] - np.asarray((x0, 0.0))))
    return vals[0], vals[1], multiplier, closure_error


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    cycles = []
    trajectories = []
    for index, (left, right) in enumerate(BRACKETS, start=1):
        root = brentq(displacement, left, right, xtol=1.0e-10, rtol=1.0e-12, maxiter=100)
        returned, period = return_data(root)
        xs, ys, multiplier, closure_error = orbit_and_multiplier(root, period)
        trajectories.append((xs, ys))
        cycles.append(
            {
                "index": index,
                "section": "y=0, downward crossing",
                "section_x": root,
                "return_residual": returned - root,
                "period": period,
                "floquet_multiplier_nonrigorous": multiplier,
                "stability": "stable" if multiplier < 1.0 else "unstable",
                "closure_error": closure_error,
            }
        )

    record = {
        "status": "NON-RIGOROUS NUMERICAL CONTROL; NOT A CERTIFICATE",
        "system": {
            "dx": "y + x^2 + x*y",
            "dy": "-10*x^2 + 2.2*x*y + 0.7*y^2 - 72.7778*x + 0.0015*y",
        },
        "method": "DOP853 plus Brent roots of a downward y=0 return map",
        "scipy_version": scipy.__version__,
        "cycles": cycles,
    }
    result_path = args.output_dir / "four_cycle_control.json"
    result_path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    colors = ("#2563eb", "#dc2626", "#16a34a")
    for (xs, ys), color, cycle in zip(trajectories[:3], colors, cycles[:3]):
        axes[0].plot(xs, ys, color=color, lw=1.25, label=f"x0={cycle['section_x']:.6g}")
    axes[0].scatter([0.0], [0.0], color="black", s=15)
    axes[0].set_title("Three nested cycles around (0, 0)")
    axes[0].set_aspect("equal", adjustable="datalim")
    axes[0].legend(fontsize=8)

    xs, ys = trajectories[3]
    axes[1].plot(xs, ys, color="#7c3aed", lw=1.1, label=f"x0={cycles[3]['section_x']:.6g}")
    axes[1].scatter([-6.2596], [7.4498], color="black", s=15)
    axes[1].set_title("Cycle around the second focus")
    axes[1].legend(fontsize=8)
    for axis in axes:
        axis.set_xlabel("x")
        axis.set_ylabel("y")
        axis.grid(alpha=0.2)
    fig.suptitle("Four-cycle regression (floating point; not proof)")
    fig.tight_layout()
    fig.savefig(args.output_dir / "four_cycle_control.png", dpi=180)

    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
