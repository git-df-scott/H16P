"""Step 2 reconnaissance for the outer-fold strike.

Part A: at the seed field, scan the section outward beyond the outermost
tracked upper root and map the displacement D(s) until the return domain
ends.  Records where the itinerary gate fails -- a missing domain, not a
zero.

Part B: two locally distinct continuation branches, capped at 50 accepted
steps each.  At every step the four seed roots are re-corrected and the
outer region is re-scanned for a fold signature (D and d_s D near zero
together, beyond the outermost tracked root).

Evidence class: NUM.  Nothing here is an interval certificate, and a small
outer |D| alone is not a cycle.
"""
import json, math, time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parent
TAU = 1e-4
TMAX = 300.0
EQ_GUARD = 1e-3
RTOL = 2e-13

# theta = (a, b, e0, e1, e2)
THETA0 = np.array([-7.0 / 4.0, 1.0 / 3.0, TAU / (-7.0 / 4.0 - 1.0),
                   -31379.0 * TAU / 25000.0, -7517.0 * TAU / 5000.0])

SEED_ROOTS = {1: [0.9802427079261901, 1.3899647291580621, 2.0452066338801704],
              -1: [9.143711643029466]}
OUTER = 2.0452066338801704   # outermost tracked upper root at the seed


def equilibria(th):
    a, b, e0, e1, e2 = th
    k = e0 / 2.0
    coeffs = [a, e1, (b - 2.0) / 4.0 + e2 * k, (1.0 - b) * k, b * k * k]
    pts = []
    for r in np.roots(coeffs):
        if abs(r.imag) < 1e-12 * max(1.0, abs(r.real)) and abs(r.real) > 1e-300:
            x = float(r.real)
            pts.append((x, k / x))
    return pts


def make_field(th):
    a, b, e0, e1, e2 = th
    eq = equilibria(th)

    def rhs_log(t, z, sign, side):
        x, w = z
        y = side * math.exp(w)
        dx = (b - 2.0) / 4.0 + (1.0 - b) * y + a * x * x + b * y * y + e1 * x + e2 * x * y
        dy = e0 - 2.0 * x * y
        return sign * np.array([dx, dy / y])

    def rhs_cart(t, z, sign):
        x, y = z
        return sign * np.array([(b - 2.0) / 4.0 + (1.0 - b) * y + a * x * x + b * y * y
                                + e1 * x + e2 * x * y, e0 - 2.0 * x * y])
    return rhs_log, rhs_cart, eq


def half_return(th, s0, side, sign):
    rhs_log, rhs_cart, eq = make_field(th)
    y0 = side * math.exp(s0)
    z = np.array([0.0, y0])
    speed = max(1.0, float(np.linalg.norm(rhs_cart(0.0, z, sign))))
    dt = min(1e-6, 1e-5 / speed)
    pre = solve_ivp(rhs_cart, [0.0, dt], z, args=(sign,), method='DOP853',
                    rtol=RTOL, atol=RTOL * 0.01)
    if not pre.success:
        raise RuntimeError('launch failed')
    x1, y1 = pre.y[:, -1]
    ev = lambda t, zz, *a_: zz[0]
    ev.terminal = True
    ev.direction = -sign
    sol = solve_ivp(rhs_log, [dt, TMAX], np.array([x1, math.log(abs(y1))]),
                    args=(sign, side), method='DOP853', events=ev,
                    rtol=RTOL, atol=RTOL * 0.01, max_step=0.2)
    if not sol.success:
        raise RuntimeError('integration failed')
    if len(sol.t_events[0]) != 1:
        raise RuntimeError('no return: %d crossings' % len(sol.t_events[0]))
    xs = sol.y[0]
    ys = side * np.exp(sol.y[1])
    core = xs[1:-1]
    if len(core) > 1 and np.any(core[:-1] * core[1:] < 0):
        raise RuntimeError('interior section crossing')
    if eq:
        dmin = min(math.hypot(xs[i] - ex, ys[i] - ey)
                   for i in range(len(xs)) for ex, ey in eq)
        if dmin < EQ_GUARD:
            raise RuntimeError('equilibrium approach %.3g' % dmin)
    return sol.y_events[0][0][1], float(sol.t_events[0][0])


def D(th, s, side):
    sf, _ = half_return(th, s, side, +1)
    sb, _ = half_return(th, s, side, -1)
    return sf - sb


def safe_D(th, s, side):
    try:
        return D(th, s, side), None
    except RuntimeError as exc:
        return None, str(exc)


def track_root(th, s_prev, side, window=0.05):
    """Re-isolate a root near s_prev; returns (s_root, derivative)."""
    lo, hi = s_prev - window, s_prev + window
    for _ in range(6):
        dlo, elo = safe_D(th, lo, side)
        dhi, ehi = safe_D(th, hi, side)
        if elo or ehi:
            raise RuntimeError('domain failure while bracketing: %s' % (elo or ehi))
        if dlo * dhi < 0:
            break
        lo -= window
        hi += window
    else:
        raise RuntimeError('root lost: no sign change within +-%.3g' % (hi - s_prev))
    r = brentq(lambda s: D(th, s, side), lo, hi, xtol=1e-12, rtol=8.9e-16, maxiter=200)
    h = 1e-4
    d = (D(th, r + h, side) - D(th, r - h, side)) / (2 * h)
    return r, d


def outer_profile(th, s_start, s_stop=6.0, n=48):
    """D on the upper section beyond s_start, plus the fold signature."""
    grid = np.linspace(s_start, s_stop, n)
    vals, fails = [], []
    for s in grid:
        d, err = safe_D(th, float(s), 1)
        if err is None:
            vals.append((float(s), d))
        else:
            fails.append({'s': float(s), 'reason': err})
    prof = {'evaluated': len(vals), 'failed': len(fails),
            'domain_end': vals[-1][0] if vals else None,
            'first_failure': fails[0] if fails else None}
    if len(vals) >= 3:
        ss = np.array([v[0] for v in vals])
        dd = np.array([v[1] for v in vals])
        prof['sign_changes'] = int(np.sum(dd[:-1] * dd[1:] < 0))
        prof['min_abs_D'] = float(np.min(np.abs(dd)))
        prof['s_at_min_abs_D'] = float(ss[int(np.argmin(np.abs(dd)))])
        slope = np.diff(dd) / np.diff(ss)
        prof['min_abs_slope'] = float(np.min(np.abs(slope)))
        # fold signature: |D| and |dD/ds| both small at the same s
        mid = 0.5 * (ss[:-1] + ss[1:])
        dmid = 0.5 * (dd[:-1] + dd[1:])
        scale = float(np.max(np.abs(dd)))
        score = np.abs(dmid) / scale + np.abs(slope) / (np.max(np.abs(slope)) + 1e-300)
        prof['fold_score'] = float(np.min(score))
        prof['s_at_fold_score'] = float(mid[int(np.argmin(score))])
        prof['profile'] = [{'s': float(s), 'D': float(d)} for s, d in vals]
    return prof


def main():
    begin = time.perf_counter()
    out = {'evidence': 'NUM; log-solver scan and bounded continuation',
           'theta0': list(THETA0), 'step_cap': 50, 'branches': []}

    print('--- Part A: seed outer scan ---', flush=True)
    out['seed_outer_profile'] = outer_profile(THETA0, OUTER + 1e-3)
    print(json.dumps({k: v for k, v in out['seed_outer_profile'].items()
                      if k != 'profile'}, indent=2), flush=True)

    print('--- Part B: continuation branches ---', flush=True)
    directions = {
        # grow the reversibility-breaking perturbation, coefficients only
        'perturbation_growth': np.array([0.0, 0.0, THETA0[2], THETA0[3], THETA0[4]]) / TAU,
        # move the quadratic coefficient a, holding the perturbation fixed
        'a_sweep': np.array([1.0, 0.0, 0.0, 0.0, 0.0]),
    }
    for name, d in directions.items():
        branch = {'name': name, 'direction': list(d), 'steps': [], 'stop_reason': None}
        th = THETA0.copy()
        roots = {1: list(SEED_ROOTS[1]), -1: list(SEED_ROOTS[-1])}
        step = 0.05 if name == 'a_sweep' else 0.25 * TAU
        for k in range(50):
            th_next = th + step * d
            try:
                new = {1: [], -1: []}
                for side in (1, -1):
                    for r in roots[side]:
                        new[side].append(track_root(th_next, r, side))
                up = sorted(r for r, _ in new[1])
                if len(up) > 1 and min(b - a for a, b in zip(up, up[1:])) < 1e-3:
                    raise RuntimeError('upper isolating intervals merged')
                if abs(th_next[2]) < 1e-12:
                    raise RuntimeError('e0 reached the invariant-line center limit')
            except RuntimeError as exc:
                branch['stop_reason'] = 'step %d rejected: %s' % (k + 1, exc)
                break
            th = th_next
            roots = {s: [r for r, _ in new[s]] for s in (1, -1)}
            prof = outer_profile(th, max(roots[1]) + 1e-3, n=24)
            row = {'step': k + 1, 'theta': list(th),
                   'upper_roots': roots[1], 'lower_roots': roots[-1],
                   'upper_derivatives': [dv for _, dv in new[1]],
                   'lower_derivatives': [dv for _, dv in new[-1]],
                   'outer': {kk: vv for kk, vv in prof.items() if kk != 'profile'}}
            branch['steps'].append(row)
            if prof.get('sign_changes'):
                branch['stop_reason'] = 'EXTRA OUTER SIGN CHANGE at step %d' % (k + 1)
                print(json.dumps(row), flush=True)
                break
        else:
            branch['stop_reason'] = 'reached the 50-step cap'
        branch['accepted_steps'] = len(branch['steps'])
        if branch['steps']:
            last = branch['steps'][-1]
            branch['final_theta'] = last['theta']
            branch['final_upper_roots'] = last['upper_roots']
            branch['final_outer'] = last['outer']
            branch['best_fold_score'] = min(
                s['outer'].get('fold_score', float('inf')) for s in branch['steps'])
        out['branches'].append(branch)
        print(json.dumps({k: v for k, v in branch.items() if k != 'steps'},
                         indent=2), flush=True)

    out['wall_seconds'] = time.perf_counter() - begin
    out['conclusion'] = {
        'extra_outer_cycle_found': any('EXTRA OUTER' in (b['stop_reason'] or '')
                                       for b in out['branches']),
        'four_cycles_held_throughout': all(
            all(len(s['upper_roots']) == 3 and len(s['lower_roots']) == 1
                for s in b['steps']) for b in out['branches']),
    }
    (ROOT / 'data' / 'step2_outer_scan.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out['conclusion'], indent=2), flush=True)
    print('seconds', out['wall_seconds'], flush=True)


if __name__ == '__main__':
    main()
