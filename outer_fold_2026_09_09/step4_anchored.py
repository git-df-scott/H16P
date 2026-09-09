"""Main experiment: preserve four roots, solve for an additional outer fold.

Replaces the slope-maximisation run (step3), whose success criterion was
wrong. Stages:

  A  Reproduce four disjoint sign brackets with the repaired gates, audit the
     section derivative AND all five parameter sensitivities, cross-check the
     log coordinate, and establish DISTINCTNESS -- the upper section
     coordinate is two-to-one on cycles, so a sign change is not a cycle.
  B  Build the anchored system F_i(theta) = D(s_i; theta, side_i) and inspect
     its scaled 4x5 Jacobian: rank, conditioning, kernel.
  C  Follow the resulting curve by pseudo-arclength predictor-corrector in
     both kernel directions, with explicit caps.
  D  At each accepted point solve the genuine six-equation fold system
     (four anchors, D(s*) = 0, dD/ds(s*) = 0) in six unknowns (theta, s*)
     by damped Newton, multi-started across the outer window. No sampled sign
     change is required; a positive sampled slope is not a criterion.
  E  Any converged fold is then tested for a nonzero second section
     derivative and a transverse unfolding derivative, and perturbed to both
     sides while the four original cycles are re-verified.

Evidence class: NUM. Nothing here certifies a periodic orbit.
"""
import json, math, os, sys, time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from machinery import (D_SCALE, DS_SCALE, SECTION_EQ_GUARD, THETA_SCALE,
                       ReturnFailure, candidate_gate, derivative_audit,
                       distinct_cycles, log_cross_check,
                       parameter_sensitivity_audit, seed_theta)
from anchored import (Budget, anchor_system, correct, evaluate, find_bracket,
                      fold_solve, fold_system, refine_root)

ROOT = Path(__file__).resolve().parent

SEED_GUESSES = [(0.9802427, 1), (1.3899488, 1), (2.0452069, 1), (9.1437116, -1)]
MAX_ACCEPTED = 50
MAX_REJECTED = 40
MAX_EVALUATIONS = int(os.environ.get('MAX_EVALUATIONS', '250000'))
MAX_SECONDS = float(os.environ.get('MAX_SECONDS', '5400'))
E0_FLOOR = 1e-9           # invariant-line class guard, see note below
OUTER_MARGINS = (0.05, 0.15, 0.30)
FOLD_STARTS = (0.1, 0.4, 0.9, 1.6, 2.6)

# Why e0 must stay away from zero: at e0 = 0 the second equation becomes
# y' = -2xy, so {y = 0} is an invariant straight line. Quadratic systems with
# an invariant straight line have at most one limit cycle, so that class
# cannot host the target and is excluded rather than searched. (It is NOT
# excluded because D vanishes there -- on the seed direction D = 4.6e-3 at
# e0 = 0, so "centre limit" was the wrong reason to give.)


def stage_a(bud, th):
    roots, audits = [], []
    for s, side in SEED_GUESSES:
        rec = refine_root(bud, th, s, side)
        rec['section_eq_distance'] = candidate_gate(th, rec['s_root'], side)[
            'section_eq_distance']
        try:
            rec['log_cross_check_gap'] = rec['D'] - log_cross_check(
                th, rec['s_root'], side)
        except ReturnFailure as exc:
            rec['log_cross_check_gap'] = exc.as_record()
        roots.append(rec)
        audits.append({'s': rec['s_root'], 'side': side,
                       'section_derivative': derivative_audit(th, rec['s_root'], side),
                       'parameter_sensitivities':
                           parameter_sensitivity_audit(th, rec['s_root'], side)})
    ups = sorted(r['s_root'] for r in roots if r['side'] == 1)
    sep = min(b - a for a, b in zip(ups, ups[1:])) if len(ups) > 1 else None
    dist = distinct_cycles(th, roots)
    return {'roots': roots, 'audits': audits,
            'min_upper_separation': sep,
            'worst_position_uncertainty': max(r['position_uncertainty'] for r in roots),
            'uncertainty_over_separation':
                max(r['position_uncertainty'] for r in roots) / sep if sep else None,
            'distinctness': dist,
            'independence_note':
                'The Cartesian and logarithmic formulations share the DOP853 '
                'stepper and the same event machinery: a coordinate '
                'cross-check, not independent verification.'}


def verify_point(bud, th, anchors):
    """Everything an accepted continuation point must satisfy."""
    out = {'theta': list(th), 'anchor_checks': [], 'ok': True, 'reasons': []}
    if abs(th[2]) < E0_FLOOR:
        out['ok'] = False
        out['reasons'].append('e0 below the invariant-line floor')
    for s, side in anchors:
        try:
            r = evaluate(bud, th, s, side)
            br = find_bracket(bud, th, s, side)
            chk = {'s': s, 'side': side, 'D_scaled': r['D'] / D_SCALE,
                   'dD_ds': r['dD_ds'], 'eq_distance': r['eq_distance'],
                   'period_estimate': r['period_estimate'],
                   'bracket': br}
            if br is None:
                out['ok'] = False
                out['reasons'].append('anchor %.6f lost its sign bracket' % s)
            if r['eq_distance'] < 1e-2:
                out['ok'] = False
                out['reasons'].append('anchor %.6f orbit near an equilibrium' % s)
        except ReturnFailure as exc:
            chk = {'s': s, 'side': side, 'failure': exc.as_record()}
            out['ok'] = False
            out['reasons'].append('anchor %.6f unresolved' % s)
        out['anchor_checks'].append(chk)
    if out['ok']:
        d = distinct_cycles(th, [{'s_root': s, 'side': side} for s, side in anchors])
        out['distinct_count'] = d['distinct_count']
        if d['distinct_count'] != len(anchors):
            out['ok'] = False
            out['reasons'].append('anchors ceased to be distinct cycles')
    return out


def outer_profile(bud, th, s_start, width=3.0, n=16):
    vals, fails = [], []
    for s in np.linspace(s_start, s_start + width, n):
        try:
            r = evaluate(bud, th, float(s), 1)
            vals.append({'s': float(s), 'D': r['D'], 'dD_ds': r['dD_ds']})
        except (ReturnFailure, RuntimeError) as exc:
            rec = exc.as_record() if isinstance(exc, ReturnFailure) else {'reason': str(exc)}
            rec['s'] = float(s)
            fails.append(rec)
    return {'values': vals, 'failures': fails}


def try_folds(bud, th, anchors, s_out, margins=OUTER_MARGINS):
    """Actual double-zero solves, multi-started, with margin sensitivity."""
    attempts = []
    for margin in margins:
        exclusion = s_out + margin
        for off in FOLD_STARTS:
            try:
                res = fold_solve(bud, th, anchors, exclusion + off, exclusion)
            except RuntimeError as exc:
                return attempts, str(exc)
            res['margin'] = margin
            res['s_star_start'] = exclusion + off
            attempts.append(res)
    return attempts, None


def run_branch(name, bud, th0, anchors, tangent0, ds, out):
    branch = {'name': name, 'arclength_step': ds, 'steps': [], 'rejected': [],
              'stop_reason': None,
              'caps': {'accepted': MAX_ACCEPTED, 'rejected': MAX_REJECTED,
                       'evaluations': MAX_EVALUATIONS, 'seconds': MAX_SECONDS}}
    th = np.array(th0, dtype=float)
    tangent = np.array(tangent0, dtype=float)
    accepted = rejected = 0
    s_out = max(s for s, side in anchors if side == 1)
    while accepted < MAX_ACCEPTED and rejected < MAX_REJECTED:
        th_pred = th + ds * tangent * THETA_SCALE
        try:
            info = correct(bud, th_pred, anchors, tangent, ds)
        except (ReturnFailure, RuntimeError) as exc:
            rejected += 1
            rec = exc.as_record() if isinstance(exc, ReturnFailure) else {'reason': str(exc)}
            branch['rejected'].append({'attempt': rejected, 'theta': list(th_pred),
                                       'failure': rec})
            if isinstance(exc, RuntimeError):
                branch['stop_reason'] = str(exc)
                break
            ds *= 0.5
            continue
        th_new = np.array(info['theta'], dtype=float)
        if not info['converged']:
            rejected += 1
            branch['rejected'].append({'attempt': rejected, 'theta': list(th_new),
                                       'failure': {'reason': 'corrector did not converge',
                                                   'residual': info['final_residual_norm']}})
            ds *= 0.5
            continue
        try:
            check = verify_point(bud, th_new, anchors)
        except RuntimeError as exc:
            branch['stop_reason'] = str(exc)
            break
        if not check['ok']:
            rejected += 1
            branch['rejected'].append({'attempt': rejected, 'theta': list(th_new),
                                       'failure': {'reason': '; '.join(check['reasons']),
                                                   'checks': check}})
            ds *= 0.5
            continue
        th = th_new
        accepted += 1
        if info['rank'] == 4:
            new_tan = np.array(info['null_vector_scaled'], dtype=float)
            if float(new_tan @ tangent) < 0:
                new_tan = -new_tan
            tangent = new_tan
        try:
            prof = outer_profile(bud, th, s_out + OUTER_MARGINS[0])
            folds, halt = try_folds(bud, th, anchors, s_out)
        except RuntimeError as exc:
            branch['stop_reason'] = str(exc)
            break
        converged = [f for f in folds if f.get('status') == 'converged']
        row = {'step': accepted, 'theta': list(th), 'arclength_step': ds,
               'anchor_rank': info['rank'],
               'anchor_singular_values': info['singular_values'],
               'anchor_condition_on_range': info['condition_on_range'],
               'anchor_residual_norm_scaled': info['residual_norm_scaled'],
               'verification': check, 'outer_profile': prof,
               'fold_attempts': folds,
               'converged_folds': converged}
        branch['steps'].append(row)
        print(json.dumps({'branch': name, 'step': accepted,
                          'theta': [float('%.6g' % v) for v in th],
                          'rank': info['rank'],
                          'fold_attempts': len(folds),
                          'converged_folds': len(converged),
                          'best_fold_residual': min(
                              [f.get('final_residual_norm_scaled')
                               or (f['history'][-1]['residual_norm_scaled']
                                   if f.get('history') else None)
                               for f in folds if f.get('history')] or [None],
                              key=lambda v: (v is None, v)),
                          'evals': bud.evaluations}), flush=True)
        if converged:
            branch['stop_reason'] = 'fold solve converged at step %d' % accepted
            break
        ds = math.copysign(min(abs(ds) * 1.3, abs(out['arclength_base'])), ds)
    else:
        branch['stop_reason'] = ('accepted cap reached' if accepted >= MAX_ACCEPTED
                                 else 'rejected cap reached')
    branch['accepted_steps'] = accepted
    branch['rejected_steps'] = rejected
    out['branches'].append(branch)
    return branch


def main():
    t0 = time.perf_counter()
    bud = Budget(MAX_EVALUATIONS, MAX_SECONDS)
    th0 = seed_theta()
    out = {'evidence': 'NUM; double precision, no interval arithmetic, no '
                       'validated integration, no certified periodic orbit',
           'fixed_scales': {'D_SCALE': D_SCALE, 'DS_SCALE': DS_SCALE,
                            'THETA_SCALE': THETA_SCALE.tolist()},
           'caps': {'accepted_per_branch': MAX_ACCEPTED,
                    'rejected_per_branch': MAX_REJECTED,
                    'return_evaluations': MAX_EVALUATIONS,
                    'wall_seconds': MAX_SECONDS},
           'theta0': th0.tolist(), 'branches': []}

    print('--- Stage A: seed reproduction, audits, distinctness ---', flush=True)
    A = stage_a(bud, th0)
    out['stage_a'] = A
    print(json.dumps({'distinct_cycles': A['distinctness']['distinct_count'],
                      'min_upper_separation': A['min_upper_separation'],
                      'worst_position_uncertainty': A['worst_position_uncertainty'],
                      'uncertainty_over_separation': A['uncertainty_over_separation'],
                      'worst_parameter_relative_spread':
                          max(a['parameter_sensitivities']['worst_relative_spread']
                              for a in A['audits'])}, indent=2), flush=True)

    anchors = [(r['s_root'], r['side']) for r in A['roots']]
    print('--- Stage B: anchored 4x5 system ---', flush=True)
    info = anchor_system(bud, th0, anchors)
    out['stage_b'] = {'rank': info['rank'],
                      'singular_values': info['singular_values'],
                      'condition_on_range': info['condition_on_range'],
                      'null_vector_scaled': info['null_vector_scaled'],
                      'residual_norm_scaled': info['residual_norm_scaled'],
                      'interpretation':
                          'rank 4 gives a one-dimensional kernel, so the '
                          'anchor-preserving set is locally a curve. Rank alone '
                          'says nothing about the six-equation fold system.'}
    print(json.dumps(out['stage_b'], indent=2), flush=True)
    if info['rank'] != 4:
        out['halt'] = 'anchored Jacobian is not rank 4; continuation not attempted'
    else:
        tangent = np.array(info['null_vector_scaled'], dtype=float)
        out['arclength_base'] = 0.05
        print('--- Stage C/D: pseudo-arclength continuation + fold solves ---',
              flush=True)
        run_branch('kernel_forward', bud, th0, anchors, tangent, 0.05, out)
        run_branch('kernel_backward', bud, th0, anchors, -tangent, 0.05, out)

    out['evaluations'] = bud.evaluations
    out['wall_seconds'] = time.perf_counter() - t0
    all_folds = [f for b in out['branches'] for s in b['steps']
                 for f in s['fold_attempts']]
    out['summary'] = {
        'branches': [{'name': b['name'], 'accepted': b['accepted_steps'],
                      'rejected': b['rejected_steps'],
                      'stop_reason': b['stop_reason']} for b in out['branches']],
        'fold_solves_attempted': len(all_folds),
        'fold_solves_converged': sum(1 for f in all_folds
                                     if f.get('status') == 'converged'),
        'fold_solve_statuses': sorted({f.get('status') for f in all_folds}),
        'criterion': 'a fold candidate requires a CONVERGED six-equation solve '
                     'with small residual in the fixed scales, then a nonzero '
                     'second section derivative and transverse unfolding '
                     'derivative, then additional DISTINCT brackets in one '
                     'frozen field. A positive sampled slope is not a criterion.',
        'five_distinct_cycles_found': False,
    }
    (ROOT / 'data' / 'step4_anchored.json').write_text(
        json.dumps(out, indent=2, default=str) + '\n')
    print(json.dumps(out['summary'], indent=2), flush=True)
    print('evaluations', bud.evaluations, 'wall %.1f s' % out['wall_seconds'])


if __name__ == '__main__':
    main()
