"""Step 1 of the outer-fold strike: re-establish the four-cycle seed.

Two independent solvers on the SAME field:
  * Cartesian: state (x, y), as in reversible_reseed/verify_control.py.
  * Logarithmic: state (x, w) with w = log|y|, a different ODE.

Section: {x = 0}. Section coordinate s = log|y|. Displacement
  D(s) = s_forward(s) - s_backward(s),
the forward/backward log-height mismatch of the reversible half returns.

The four saved sign brackets are refined to isolated roots, each root is
given a section derivative, and every half return is gated on a complete
itinerary (exactly one transverse section crossing, no equilibrium
approach, both solvers agreeing).

Evidence class: NUM. No interval certificate, no fifth cycle asserted.
"""
import json, math, time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parent

A = -7.0 / 4.0
B = 1.0 / 3.0
C = 31379.0 / 25000.0
M = 7517.0 / 5000.0
TAU = 1e-4
E0 = TAU / (A - 1.0)          # = -4*tau/11
TMAX = 300.0
EQ_GUARD = 1e-3               # minimum allowed distance to any equilibrium

# Saved section witnesses (|y| at x = 0) from reversible_reseed/data.
SEED = {
    1: [1.650729931857833789905966566951675683954534825223876907,
        3.484771466591614075358386897082783768858247531297605364,
        7.31241833621300096982590430940074288047331783457839216,
        28.17882337367149048058515156640056103908504001403481686],
    -1: [3161.394763622164286466009270150545516911614143570237057,
         20733.33317025880876267330185953580126290330530638437481],
}


def equilibria():
    """Real solutions of P = Q = 0 for the perturbed field."""
    pts = []
    # Q = 0: -2xy + E0 = 0.  x = 0 is impossible (E0 != 0), so y = E0/(2x).
    # Substitute into P and clear denominators -> quintic in x.
    k = E0 / 2.0
    # P = (B-2)/4 + (1-B)*k/x + A x^2 + B k^2/x^2 - TAU*C*x - TAU*M*k
    # multiply by x^2:
    coeffs = [A, -TAU * C, (B - 2.0) / 4.0 - TAU * M * k, (1.0 - B) * k, B * k * k]
    for r in np.roots(coeffs):
        if abs(r.imag) < 1e-12 * max(1.0, abs(r.real)):
            x = float(r.real)
            pts.append((x, k / x))
    return pts


EQ = equilibria()


def rhs_cartesian(t, z, sign):
    x, y = z
    return sign * np.array([
        (B - 2.0) / 4.0 + (1.0 - B) * y + A * x * x + B * y * y - TAU * (C * x + M * x * y),
        -2.0 * x * y + E0,
    ])


def rhs_log(t, z, sign, side):
    """side = +1 -> y = exp(w) > 0;  side = -1 -> y = -exp(w) < 0."""
    x, w = z
    y = side * math.exp(w)
    dx = (B - 2.0) / 4.0 + (1.0 - B) * y + A * x * x + B * y * y - TAU * (C * x + M * x * y)
    dy = -2.0 * x * y + E0
    return sign * np.array([dx, dy / y])


def _event_x(t, z, *args):
    return z[0]


def _launch(y0, sign, rtol):
    """Leave the section along the flow before arming the return event."""
    z = np.array([0.0, y0])
    speed = max(1.0, float(np.linalg.norm(rhs_cartesian(0.0, z, sign))))
    dt = min(1e-6, 1e-5 / speed)
    sol = solve_ivp(rhs_cartesian, [0.0, dt], z, args=(sign,),
                    method='DOP853', rtol=rtol, atol=rtol * 0.01)
    if not sol.success:
        raise RuntimeError(sol.message)
    return dt, sol.y[:, -1]


def half_return(s0, side, sign, rtol, solver):
    """Half return from (0, side*exp(s0)); returns log|y| at the next crossing."""
    y0 = side * math.exp(s0)
    dt, z1 = _launch(y0, sign, rtol)
    ev = lambda t, z, *a: z[0]
    ev.terminal = True
    ev.direction = -sign
    if solver == 'cartesian':
        f, z_start, args = rhs_cartesian, z1, (sign,)
    else:
        f, z_start, args = rhs_log, np.array([z1[0], math.log(abs(z1[1]))]), (sign, side)
    sol = solve_ivp(f, [dt, TMAX], z_start, args=args, method='DOP853',
                    events=ev, rtol=rtol, atol=rtol * 0.01, max_step=0.2,
                    dense_output=False)
    if not sol.success:
        raise RuntimeError('integration failed: ' + sol.message)
    if len(sol.t_events[0]) != 1:
        raise RuntimeError('incomplete itinerary: %d section crossings' % len(sol.t_events[0]))
    xs = sol.y[0]
    ys = sol.y[1] if solver == 'cartesian' else side * np.exp(sol.y[1])
    dmin = min(math.hypot(xs[i] - ex, ys[i] - ey)
               for i in range(len(xs)) for ex, ey in EQ)
    core = xs[1:-1]  # drop the launch point and the terminal event point
    interior_sign_changes = int(np.sum(core[:-1] * core[1:] < 0)) if len(core) > 1 else 0
    if interior_sign_changes:
        raise RuntimeError('interior section crossing before the return')
    zend = sol.y_events[0][0]
    s_end = math.log(abs(zend[1])) if solver == 'cartesian' else zend[1]
    return {'s_end': s_end, 'time': float(sol.t_events[0][0]),
            'eq_distance': dmin, 'max_abs_y': float(np.max(np.abs(ys)))}


def displacement(s0, side, rtol, solver):
    fwd = half_return(s0, side, +1, rtol, solver)
    bwd = half_return(s0, side, -1, rtol, solver)
    if min(fwd['eq_distance'], bwd['eq_distance']) < EQ_GUARD:
        raise RuntimeError('orbit passes within %g of an equilibrium' % EQ_GUARD)
    return fwd['s_end'] - bwd['s_end'], fwd, bwd


def main():
    begin = time.perf_counter()
    out = {
        'evidence': 'NUM; two independent non-interval solvers, no fifth cycle asserted',
        'field': {'a': '-7/4', 'b': '1/3', 'tau': TAU,
                  'epsilon0': '-4*tau/11', 'epsilon1': '-31379*tau/25000',
                  'epsilon2': '-7517*tau/5000'},
        'section': 'x = 0, coordinate s = log|y|',
        'equilibria': [{'x': p[0], 'y': p[1]} for p in EQ],
        'samples': [], 'roots': [], 'gates': {},
    }

    rtol = 2e-13
    # --- sampled displacements at the saved witnesses, both solvers ---
    for side, ys in SEED.items():
        for y in ys:
            s = math.log(y)
            row = {'side': side, 's': s, 'abs_y': y}
            for solver in ('cartesian', 'log'):
                D, fwd, bwd = displacement(s, side, rtol, solver)
                row[solver] = {'D': D, 'forward_time': fwd['time'],
                               'backward_time': bwd['time'],
                               'eq_distance': min(fwd['eq_distance'], bwd['eq_distance']),
                               'max_abs_y': max(fwd['max_abs_y'], bwd['max_abs_y'])}
            row['solver_gap'] = row['cartesian']['D'] - row['log']['D']
            out['samples'].append(row)
            print(json.dumps(row), flush=True)

    # --- brackets: consecutive sign changes per side ---
    brackets = []
    for side in (1, -1):
        rows = [r for r in out['samples'] if r['side'] == side]
        for lo, hi in zip(rows, rows[1:]):
            if lo['log']['D'] * hi['log']['D'] < 0:
                brackets.append((side, lo['s'], hi['s']))
    out['gates']['brackets_found'] = len(brackets)
    out['gates']['upper_brackets'] = sum(1 for b in brackets if b[0] == 1)
    out['gates']['lower_brackets'] = sum(1 for b in brackets if b[0] == -1)

    # --- refine each root with both solvers, derivative by central differences ---
    for side, slo, shi in brackets:
        entry = {'side': side, 'bracket': [slo, shi]}
        for solver in ('cartesian', 'log'):
            f = lambda s: displacement(s, side, rtol, solver)[0]
            s_root = brentq(f, slo, shi, xtol=1e-12, rtol=8.9e-16, maxiter=200)
            derivs = {}
            for h in (1e-4, 1e-5):
                derivs['h=%g' % h] = (f(s_root + h) - f(s_root - h)) / (2 * h)
            _, fwd, bwd = displacement(s_root, side, rtol, solver)
            entry[solver] = {
                's_root': s_root, 'abs_y_root': math.exp(s_root),
                'residual_D': f(s_root), 'derivative': derivs,
                'forward_time': fwd['time'], 'backward_time': bwd['time'],
                'period_estimate': fwd['time'] + bwd['time'],
                'eq_distance': min(fwd['eq_distance'], bwd['eq_distance']),
            }
        entry['root_gap'] = entry['cartesian']['s_root'] - entry['log']['s_root']
        entry['derivative_gap'] = (entry['cartesian']['derivative']['h=0.0001']
                                   - entry['log']['derivative']['h=0.0001'])
        out['roots'].append(entry)
        print(json.dumps(entry), flush=True)

    # --- disjointness and nondegeneracy gates ---
    upper = sorted(r['log']['s_root'] for r in out['roots'] if r['side'] == 1)
    lower = sorted(r['log']['s_root'] for r in out['roots'] if r['side'] == -1)
    out['gates']['upper_roots'] = len(upper)
    out['gates']['lower_roots'] = len(lower)
    out['gates']['min_upper_separation'] = (min(b - a for a, b in zip(upper, upper[1:]))
                                            if len(upper) > 1 else None)
    out['gates']['all_derivatives_nonzero'] = all(
        abs(r[s]['derivative']['h=0.0001']) > 1e-8 for r in out['roots']
        for s in ('cartesian', 'log'))
    out['gates']['max_solver_root_gap'] = max(abs(r['root_gap']) for r in out['roots'])
    out['gates']['four_cycles_reproduced'] = (
        out['gates']['upper_roots'] == 3 and out['gates']['lower_roots'] == 1
        and out['gates']['all_derivatives_nonzero'])
    out['wall_seconds'] = time.perf_counter() - begin

    (ROOT / 'data' / 'step1_seed.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out['gates'], indent=2), flush=True)
    print('seconds', out['wall_seconds'], flush=True)
    if not out['gates']['four_cycles_reproduced']:
        raise SystemExit('CONTROL FAILED: do not start a fold search')


if __name__ == '__main__':
    main()
