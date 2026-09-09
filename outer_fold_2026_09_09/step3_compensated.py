"""Repaired seed audit, repaired root tracking, and compensated continuation.

Stage 1  Re-audit the four seed roots with event-corrected variational
         derivatives cross-checked against finite differences at several
         step sizes and tolerances, and report the spread as an empirical
         uncertainty rather than a nonvanishing certificate.

Stage 2  Track the four roots with tangent prediction from the analytic
         sensitivities and adaptive step reduction, so that a rejected step
         is attributed to the tracker before it is attributed to the field.
         Trackers are kept apart by an explicit separation test.

Stage 3  Look beyond the outermost upper root for a DOUBLE zero of the
         return mismatch: D = 0 and dD/ds = 0 together.  A double zero need
         not change sign, so nothing here waits for a sampled sign change.
         The continuation drives `a` and compensates with (b, e0, e1, e2)
         along the scaled gradient of the outer slope functional; the fold
         Jacobian's rank and conditioning are reported at every solve.

All residuals are reported against the fixed scales declared in
machinery.py, never against a scale recomputed from the profile at hand.

Evidence class: NUM throughout.  Nothing here certifies existence or
isolation of any periodic orbit.
"""
import json, math, os, sys, time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from machinery import (D_SCALE, DS_SCALE, THETA_SCALE, ReturnFailure,
                       derivative_audit, displacement, log_cross_check,
                       seed_theta)

ROOT = Path(__file__).resolve().parent

SEED_ROOTS = {1: [0.9802427079261901, 1.3899647291580621, 2.0452066338801704],
              -1: [9.143711643029466]}

MAX_ACCEPTED = 50
MAX_REJECTED = 60
MAX_HALVINGS = 10
MAX_EVALS = 40000
MIN_SEPARATION = 1e-3
OUTER_MARGIN = 5e-3       # start of the outer window past the outermost root
OUTER_WIDTH = 4.0         # width of the outer window in s
OUTER_POINTS = 22

EVALS = {'n': 0}


def ev_displacement(th, s, side, rtol=1e-12, variational=True):
    EVALS['n'] += 1
    if EVALS['n'] > MAX_EVALS:
        raise RuntimeError('return-evaluation cap reached')
    return displacement(th, s, side, rtol, variational=variational)


# Newton cannot drive |D| below the integrator's own noise, so convergence is
# declared at the noise floor rather than at an arbitrary tight tolerance.
D_NOISE = 2e-14
S_TOL = 1e-10


def guarded_newton(th, s_pred, side, guard, tol=S_TOL, max_iter=25):
    """Newton on D confined to [s_pred - guard, s_pred + guard].

    The guard is what keeps one tracker from falling into another root's
    basin; leaving it is reported as a tracker failure, not as a statement
    about the field.
    """
    lo, hi = s_pred - guard, s_pred + guard
    s = s_pred
    for _ in range(max_iter):
        r = ev_displacement(th, s, side)
        if abs(r['dD_ds']) < 1e-14:
            raise ReturnFailure('section derivative underflow during correction')
        step = -r['D'] / r['dD_ds']
        s_new = s + step
        if s_new < lo or s_new > hi:
            raise ReturnFailure('correction left the guard interval '
                                '(|step| = %.3e, guard = %.3e)' % (abs(step), guard))
        converged = abs(s_new - s) < tol or abs(r['D']) < D_NOISE
        s = s_new
        if converged:
            break
    else:
        raise ReturnFailure('correction did not converge inside the guard '
                            '(residual stalled at the integrator noise floor)')
    final = ev_displacement(th, s, side)
    return {'s_root': s, 'D': final['D'], 'D_scaled': final['D'] / D_SCALE,
               'dD_ds': final['dD_ds'], 'dD_dtheta': final['dD_dtheta'].tolist(),
               'period_estimate': final['period_estimate'],
               'eq_distance': final['eq_distance'],
               'terminal_transversality': final['min_terminal_transversality'],
               'predicted_s': s_pred, 'prediction_error': s - s_pred}


def seed_records(th, roots):
    """Locate the roots at the seed with a wide guard, and record sensitivities."""
    recs = {1: [], -1: []}
    for side in (1, -1):
        for s in roots[side]:
            recs[side].append(guarded_newton(th, s, side, guard=0.05))
    return recs


def predict(recs, delta):
    """Tangent prediction s -> s - (dD/dtheta . delta) / (dD/ds)."""
    pred = {1: [], -1: []}
    for side in (1, -1):
        for rec in recs[side]:
            g = np.array(rec['dD_dtheta'])
            pred[side].append(rec['s_root'] - float(g @ delta) / rec['dD_ds'])
    return pred


def track_all(th, recs, delta):
    """Predict from sensitivities, then correct inside identity-preserving guards."""
    pred = predict(recs, delta)
    new = {1: [], -1: []}
    for side in (1, -1):
        order = sorted(range(len(pred[side])), key=lambda i: pred[side][i])
        for i, s_pred in enumerate(pred[side]):
            neighbours = [pred[side][j] for j in range(len(pred[side])) if j != i]
            gap = min((abs(s_pred - o) for o in neighbours), default=1.0)
            guard = min(0.35 * gap, 0.15)
            if guard < MIN_SEPARATION:
                raise ReturnFailure('predicted roots too close to guard apart')
            new[side].append(guarded_newton(th, s_pred, side, guard))
        got = [r['s_root'] for r in new[side]]
        if sorted(range(len(got)), key=lambda i: got[i]) != order:
            raise ReturnFailure('tracker order changed: identities not preserved')
        if len(got) > 1:
            gaps = [b - a for a, b in zip(sorted(got), sorted(got)[1:])]
            if min(gaps) < MIN_SEPARATION:
                raise ReturnFailure('two trackers converged: min gap %.3e' % min(gaps))
    return new


def outer_profile(th, s_start, width=OUTER_WIDTH, n=OUTER_POINTS):
    """Full outer profile: D and the analytic dD/ds, with failures kept."""
    grid = np.linspace(s_start, s_start + width, n)
    values, failures = [], []
    for s in grid:
        try:
            r = ev_displacement(th, float(s), 1)
            values.append({'s': float(s), 'D': r['D'], 'dD_ds': r['dD_ds']})
        except (ReturnFailure, RuntimeError) as exc:
            failures.append({'s': float(s), 'reason': str(exc)})
    return {'values': values, 'failures': failures,
            'evaluated': len(values), 'failed': len(failures)}


def slope_functional(profile):
    """sup of dD/ds over the outer window: the quantity a fold must lift to 0."""
    if not profile['values']:
        return None
    best = max(profile['values'], key=lambda r: r['dD_ds'])
    return {'sup_dD_ds': best['dD_ds'], 's_at_sup': best['s'], 'D_at_sup': best['D'],
            'sup_scaled': best['dD_ds'] / DS_SCALE,
            'D_at_sup_scaled': best['D'] / D_SCALE}


def slope_gradient(th, s_probe, rel=1e-3):
    """d(dD/ds)/dtheta by differencing the ANALYTIC dD/ds. Scaled by THETA_SCALE."""
    g = np.zeros(5)
    for i in range(5):
        h = rel * THETA_SCALE[i]
        tp, tm = th.copy(), th.copy()
        tp[i] += h
        tm[i] -= h
        gp = ev_displacement(tp, s_probe, 1)['dD_ds']
        gm = ev_displacement(tm, s_probe, 1)['dD_ds']
        g[i] = (gp - gm) / (2 * h)
    return g * THETA_SCALE          # gradient with respect to scaled coefficients


def fold_jacobian(th, s, rel=1e-3):
    """Jacobian of F = (D, dD/ds) in (s, scaled theta); returns rank and cond."""
    base = ev_displacement(th, s, 1)
    hs = 1e-5
    dp = ev_displacement(th, s + hs, 1)
    dm = ev_displacement(th, s - hs, 1)
    d2D_ds2 = (dp['dD_ds'] - dm['dD_ds']) / (2 * hs)
    cols = [np.array([base['dD_ds'] / D_SCALE, d2D_ds2 / DS_SCALE])]
    for i in range(5):
        h = rel * THETA_SCALE[i]
        tp, tm = th.copy(), th.copy()
        tp[i] += h
        tm[i] -= h
        rp, rm = ev_displacement(tp, s, 1), ev_displacement(tm, s, 1)
        cols.append(np.array([(rp['D'] - rm['D']) / (2 * h) * THETA_SCALE[i] / D_SCALE,
                              (rp['dD_ds'] - rm['dD_ds']) / (2 * h)
                              * THETA_SCALE[i] / DS_SCALE]))
    J = np.column_stack(cols)
    sv = np.linalg.svd(J, compute_uv=False)
    return {'jacobian': J.tolist(), 'singular_values': sv.tolist(),
            'rank': int(np.sum(sv > sv[0] * 1e-10)) if sv[0] > 0 else 0,
            'condition': float(sv[0] / sv[-1]) if sv[-1] > 0 else float('inf'),
            'residual': [base['D'] / D_SCALE, base['dD_ds'] / DS_SCALE],
            'd2D_ds2': d2D_ds2,
            'note': 'columns are (s, a, b, e0, e1, e2) in scaled coefficients; '
                    'residuals use the fixed scales from machinery.py'}


def run_branch(name, th0, recs0, driver_step, compensate, out):
    """One continuation branch. `a` is the driver, (b, e0, e1, e2) compensate.

    The compensating direction is the scaled gradient of the outer slope
    functional, so the free parameter is the driver step alone; nothing is
    imposed as an overdetermined system on the four tracked roots. They are
    tracked, not constrained.
    """
    branch = {'name': name, 'driver_step': driver_step, 'compensate': compensate,
              'constraints': 'a driven; (b,e0,e1,e2) move along the unit scaled '
                             'gradient of sup dD/ds over the outer window; the '
                             'four roots are tracked, never constrained',
              'steps': [], 'rejected': [], 'stop_reason': None,
              'caps': {'accepted': MAX_ACCEPTED, 'rejected': MAX_REJECTED,
                       'halvings': MAX_HALVINGS, 'return_evaluations': MAX_EVALS}}
    th = th0.copy()
    recs = {s: list(v) for s, v in recs0.items()}
    accepted = rejected = halvings = 0
    step = driver_step

    while accepted < MAX_ACCEPTED and rejected < MAX_REJECTED:
        roots_now = [r['s_root'] for r in recs[1]]
        prof = outer_profile(th, max(roots_now) + OUTER_MARGIN)
        sf = slope_functional(prof)
        if sf is None:
            branch['stop_reason'] = ('outer window unresolved after %d accepted steps'
                                     % accepted)
            break
        try:
            g = slope_gradient(th, sf['s_at_sup'])
        except (ReturnFailure, RuntimeError) as exc:
            branch['stop_reason'] = 'gradient unavailable: %s' % exc
            break

        direction = np.zeros(5)
        direction[0] = 1.0
        if compensate:
            comp = g[1:]
            nrm = float(np.linalg.norm(comp))
            if nrm > 0:
                direction[1:] = comp / nrm
        delta = step * direction * THETA_SCALE

        th_try = th + delta
        try:
            new = track_all(th_try, recs, delta)
        except (ReturnFailure, RuntimeError) as exc:
            rejected += 1
            halvings += 1
            branch['rejected'].append({'attempt': rejected, 'step': step,
                                       'halvings': halvings,
                                       'theta': list(th_try), 'reason': str(exc)})
            step *= 0.5
            if halvings > MAX_HALVINGS:
                branch['stop_reason'] = (
                    'step halved %d times without an accepted step at a = %.6f; '
                    'this is a tracker limit, not evidence that a cycle ceased to '
                    'exist' % (MAX_HALVINGS, th[0]))
                break
            continue

        th = th_try
        recs = new
        accepted += 1
        halvings = 0
        roots_new = [r['s_root'] for r in recs[1]]
        prof_new = outer_profile(th, max(roots_new) + OUTER_MARGIN)
        sf_new = slope_functional(prof_new)
        branch['steps'].append({
            'step': accepted, 'theta': list(th), 'step_size': step,
            'gradient_scaled': g.tolist(),
            'upper_records': recs[1], 'lower_records': recs[-1],
            'upper_roots': roots_new,
            'lower_roots': [r['s_root'] for r in recs[-1]],
            'outer': {'evaluated': prof_new['evaluated'], 'failed': prof_new['failed'],
                      'failures': prof_new['failures'], 'profile': prof_new['values']},
            'slope_functional': sf_new})
        print(json.dumps({'branch': name, 'step': accepted, 'a': th[0],
                          'sup_dD_ds_scaled': None if sf_new is None
                          else sf_new['sup_scaled'],
                          'upper_roots': [round(r, 6) for r in roots_new],
                          'evals': EVALS['n']}), flush=True)
        step = math.copysign(min(abs(step) * 1.4, abs(driver_step)), driver_step)

        if sf_new and sf_new['sup_dD_ds'] > 0:
            branch['stop_reason'] = ('outer slope turned positive at step %d; '
                                     'fold solve warranted' % accepted)
            break
    else:
        branch['stop_reason'] = ('accepted-step cap reached' if accepted >= MAX_ACCEPTED
                                 else 'rejected-step cap reached')

    branch['accepted_steps'] = accepted
    branch['rejected_steps'] = rejected
    scored = [s for s in branch['steps'] if s['slope_functional']]
    if scored:
        best = max(scored, key=lambda s: s['slope_functional']['sup_dD_ds'])
        branch['best_step'] = best['step']
        branch['best_slope_functional'] = best['slope_functional']
        try:
            branch['fold_jacobian_at_best'] = fold_jacobian(
                np.array(best['theta']), best['slope_functional']['s_at_sup'])
        except (ReturnFailure, RuntimeError) as exc:
            branch['fold_jacobian_at_best'] = {'error': str(exc)}
    out['branches'].append(branch)
    return branch


def main():
    t0 = time.perf_counter()
    th0 = seed_theta()
    out = {'evidence': 'NUM; double precision, no interval arithmetic, '
                       'no validated integration, no certified isolation',
           'independence': 'The Cartesian and logarithmic formulations share the '
                           'DOP853 stepper and the same event machinery. Their '
                           'agreement is a coordinate cross-check, not an '
                           'independent numerical verification.',
           'fixed_scales': {'D_SCALE': D_SCALE, 'DS_SCALE': DS_SCALE,
                            'THETA_SCALE': THETA_SCALE.tolist()},
           'theta0': th0.tolist(), 'seed_audit': [], 'branches': []}

    print('--- Stage 1: repaired seed audit ---', flush=True)
    recs0 = seed_records(th0, SEED_ROOTS)
    for side in (1, -1):
        for rec in recs0[side]:
            s_root = rec['s_root']
            audit = derivative_audit(th0, s_root, side)
            try:
                gap = ev_displacement(th0, s_root, side, variational=False)['D'] \
                    - log_cross_check(th0, s_root, side)
            except ReturnFailure as exc:
                gap = 'unresolved: %s' % exc
            row = {'side': side, 's_root': s_root, 'abs_y': math.exp(s_root),
                   'D': rec['D'], 'D_scaled': rec['D'] / D_SCALE,
                   'dD_ds_variational': rec['dD_ds'],
                   'derivative_estimate': audit['estimate'],
                   'derivative_spread': audit['spread'],
                   'derivative_relative_spread': audit['relative_spread'],
                   'derivative_rows': audit['rows'],
                   'period_estimate': rec['period_estimate'],
                   'eq_distance': rec['eq_distance'],
                   'terminal_transversality': rec['terminal_transversality'],
                   'log_coordinate_cross_check_gap': gap}
            out['seed_audit'].append(row)
            print(json.dumps({k: row[k] for k in
                              ('side', 's_root', 'D_scaled', 'dD_ds_variational',
                               'derivative_relative_spread',
                               'log_coordinate_cross_check_gap')}), flush=True)

    print('--- Stage 2/3: compensated continuation ---', flush=True)
    run_branch('a_up_compensated', th0, recs0, +2e-3, True, out)
    run_branch('a_down_compensated', th0, recs0, -2e-3, True, out)

    out['return_evaluations'] = EVALS['n']
    out['wall_seconds'] = time.perf_counter() - t0
    out['summary'] = {
        'branches': [{'name': b['name'], 'accepted': b['accepted_steps'],
                      'rejected': b['rejected_steps'], 'stop_reason': b['stop_reason'],
                      'best_sup_dD_ds_scaled': (b.get('best_slope_functional') or {})
                      .get('sup_scaled')} for b in out['branches']],
        'double_zero_found': any(
            (b.get('best_slope_functional') or {}).get('sup_dD_ds', -1) > 0
            for b in out['branches']),
    }
    (ROOT / 'data' / 'step3_compensated.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out['summary'], indent=2), flush=True)
    print('return evaluations', EVALS['n'], 'seconds %.1f' % out['wall_seconds'],
          flush=True)


if __name__ == '__main__':
    main()
