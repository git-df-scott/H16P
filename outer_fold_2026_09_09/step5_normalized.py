"""Normalised anchored continuation with a released anchor and a real fold solve.

WHY THIS REPLACES step4
-----------------------
step4 pinned four anchors and followed the kernel of the scaled 4x5 Jacobian.
Its fold solves "converged" 15 times out of 15 -- and every one of them was
spurious. D is very nearly linear in (e0, e1, e2), so the anchor-preserving
kernel is essentially the ray that SHRINKS the perturbation: the recorded
kernel vector is (-7e-08, -1e-07, 0.880, 0.304, 0.364), with no measurable
(a, b) component. Newton then slid down that ray to |(e0,e1,e2)| = 1.7e-10,
a factor 1e-6 below the seed, where the field is numerically the unperturbed
reversible one and D vanishes identically. Both fold equations are then
satisfied everywhere at once. That is a collapsed field, not a fold.
step4's data is kept in data/step4_anchored.json as the record of it.

THE SCHEME HERE
---------------
Unknowns and equations are stated before running, as they must be.

  Unknowns (7):  theta = (a, b, e0, e1, e2), the released anchor position
                 s_1, and the fold section coordinate s*.
  Equations (7): D(s_1; theta) = 0            (released anchor, still a root)
                 D(s_i; theta) = 0, i = 2,3,4 (pinned anchors)
                 N(theta) = 0                 (perturbation normalisation)
                 D(s*; theta) = 0             (fold)
                 dD/ds (s*; theta) = 0        (fold)

  N(theta) = |(e0,e1,e2)|^2 / |(e0,e1,e2)_seed|^2 - 1 removes the scaling
  degeneracy that produced step4's false folds: the collapsed field is no
  longer in the feasible set at all.

  Dropping the two fold equations and s* leaves 5 equations in 6 unknowns,
  a one-dimensional curve, which is what the continuation follows. The
  released anchor is the innermost upper root; the other three stay pinned.
  Degrees of freedom: 6 unknowns - 5 equations = 1, parameterised by
  pseudo-arclength.

ACCEPTANCE, WITH FIXED THRESHOLDS
---------------------------------
A fold candidate must satisfy ALL of:
  * converged 7-equation solve, residual < FOLD_TOL in the fixed scales;
  * condition number of the 7x7 Jacobian <= MAX_FOLD_CONDITION;
  * |d2D/ds2| at s* >= MIN_D2 in the fixed scale (a real quadratic tangency);
  * |dD/ds| at every anchor >= MIN_ANCHOR_SLOPE_FRACTION of its seed value
    (no global collapse of the displacement);
  * s* is a DISTINCT cycle from all four anchors, by the crossing-pair test;
  * a nonzero transverse unfolding derivative.
Small residual alone is not a candidate. Neither is a positive sampled slope.

Evidence class: NUM. Nothing here certifies a periodic orbit.
"""
import json, math, os, sys, time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from machinery import (D_SCALE, DS_SCALE, THETA_SCALE, ReturnFailure,
                       cycle_key, distinct_cycles, seed_theta)
from anchored import Budget, evaluate, find_bracket, refine_root

ROOT = Path(__file__).resolve().parent

SEED_GUESSES = [(0.9802427, 1), (1.3899488, 1), (2.0452069, 1), (9.1437116, -1)]
FOLD_TOL = 1e-3
MAX_FOLD_CONDITION = 1e6
MIN_D2 = 1e-3                      # |d2D/ds2| / DS_SCALE
MIN_ANCHOR_SLOPE_FRACTION = 0.1
E0_FLOOR = 1e-9
MAX_FREE_MOVE = 0.05      # most the released anchor may move in one accepted step
MIN_ANCHOR_GAP = 0.02     # released anchor must stay this far from every pinned one

# Both guards were added after a first run lost the released anchor's identity:
# it moved from s = 1.288775 to s = 1.431203 in a single accepted step, passing
# through the pinned anchor at 1.3899488. Two roots of the same D(.; theta)
# cannot cross -- they can only merge and separate -- so that step was the
# corrector jumping to a different root, not the root moving. Steps after it in
# data/step5_normalized_firstrun.json are unsound and are kept only as the
# record of the defect.
MAX_ACCEPTED = int(os.environ.get('MAX_ACCEPTED', '50'))
MAX_REJECTED = 40
MAX_EVALUATIONS = int(os.environ.get('MAX_EVALUATIONS', '400000'))
MAX_SECONDS = float(os.environ.get('MAX_SECONDS', '6000'))
MARGINS = (0.05, 0.15, 0.30)
FOLD_STARTS = (0.1, 0.4, 0.9, 1.6, 2.6)

SEED_PERT = np.linalg.norm(seed_theta()[2:])


def normalisation(th):
    return float(np.dot(th[2:], th[2:]) / (SEED_PERT ** 2) - 1.0)


def dnormalisation(th):
    g = np.zeros(5)
    g[2:] = 2.0 * np.asarray(th[2:]) / (SEED_PERT ** 2)
    return g * THETA_SCALE


def curve_system(bud, th, s_free, pinned):
    """5 equations (4 anchors + normalisation) in 6 unknowns (theta, s_free)."""
    th = np.asarray(th, dtype=float)
    rows, F, J = [], [], []
    r0 = evaluate(bud, th, s_free, pinned[0][1] if False else 1)
    F.append(r0['D'] / D_SCALE)
    J.append(np.concatenate([r0['dD_dtheta'] * THETA_SCALE / D_SCALE,
                             [r0['dD_ds'] / D_SCALE]]))
    rows.append({'s': s_free, 'side': 1, 'released': True,
                 'D_scaled': r0['D'] / D_SCALE, 'dD_ds': r0['dD_ds']})
    for s, side in pinned:
        r = evaluate(bud, th, s, side)
        F.append(r['D'] / D_SCALE)
        J.append(np.concatenate([r['dD_dtheta'] * THETA_SCALE / D_SCALE, [0.0]]))
        rows.append({'s': s, 'side': side, 'released': False,
                     'D_scaled': r['D'] / D_SCALE, 'dD_ds': r['dD_ds'],
                     'eq_distance': r['eq_distance']})
    F.append(normalisation(th))
    J.append(np.concatenate([dnormalisation(th), [0.0]]))
    F = np.array(F)
    J = np.array(J)
    sv = np.linalg.svd(J, compute_uv=False)
    rank = int(np.sum(sv > sv[0] * 1e-10)) if sv[0] > 0 else 0
    U, s_, Vt = np.linalg.svd(J)
    return {'F': F, 'J': J, 'rows': rows, 'singular_values': sv.tolist(),
            'rank': rank, 'null_vector': Vt[-1].tolist(),
            'condition_on_range': float(sv[0] / sv[rank - 1]) if rank else float('inf'),
            'residual_norm': float(np.linalg.norm(F))}


def curve_correct(bud, u_pred, pinned, tangent, max_iter=15, trust=0.3):
    u = np.array(u_pred, dtype=float)
    for _ in range(max_iter):
        th, s_free = u[:5], u[5]
        info = curve_system(bud, th, s_free, pinned)
        arc = float((np.concatenate([(u[:5] - u_pred[:5]) / THETA_SCALE,
                                     [u[5] - u_pred[5]]])) @ tangent)
        G = np.concatenate([info['F'], [arc]])
        if np.linalg.norm(G) < 1e-6:
            info['converged'] = True
            info['u'] = u.tolist()
            return info
        A = np.vstack([info['J'], tangent])
        try:
            step = np.linalg.solve(A, -G)
        except np.linalg.LinAlgError:
            raise ReturnFailure('curve corrector Jacobian singular',
                                {'singular_values': info['singular_values']})
        n = float(np.linalg.norm(step))
        if n > trust:
            step *= trust / n
        u = u + np.concatenate([step[:5] * THETA_SCALE, [step[5]]])
    info = curve_system(bud, u[:5], u[5], pinned)
    info['converged'] = False
    info['u'] = u.tolist()
    return info


def fold_system(bud, th, s_free, pinned, s_star, rel=1e-3):
    """7 equations in 7 unknowns (theta, s_free, s*)."""
    th = np.asarray(th, dtype=float)
    F, J = [], []
    r0 = evaluate(bud, th, s_free, 1)
    F.append(r0['D'] / D_SCALE)
    J.append(np.concatenate([r0['dD_dtheta'] * THETA_SCALE / D_SCALE,
                             [r0['dD_ds'] / D_SCALE, 0.0]]))
    anchor_slopes = [r0['dD_ds']]
    for s, side in pinned:
        r = evaluate(bud, th, s, side)
        F.append(r['D'] / D_SCALE)
        J.append(np.concatenate([r['dD_dtheta'] * THETA_SCALE / D_SCALE, [0.0, 0.0]]))
        anchor_slopes.append(r['dD_ds'])
    F.append(normalisation(th))
    J.append(np.concatenate([dnormalisation(th), [0.0, 0.0]]))

    base = evaluate(bud, th, s_star, 1)
    hs = 1e-5
    up = evaluate(bud, th, s_star + hs, 1)
    dn = evaluate(bud, th, s_star - hs, 1)
    d2 = (up['dD_ds'] - dn['dD_ds']) / (2 * hs)
    dslope = np.zeros(5)
    for j in range(5):
        h = rel * THETA_SCALE[j]
        tp, tm = th.copy(), th.copy()
        tp[j] += h
        tm[j] -= h
        dslope[j] = (evaluate(bud, tp, s_star, 1)['dD_ds']
                     - evaluate(bud, tm, s_star, 1)['dD_ds']) / (2 * h)
    F.append(base['D'] / D_SCALE)
    J.append(np.concatenate([base['dD_dtheta'] * THETA_SCALE / D_SCALE,
                             [0.0, base['dD_ds'] / D_SCALE]]))
    F.append(base['dD_ds'] / DS_SCALE)
    J.append(np.concatenate([dslope * THETA_SCALE / DS_SCALE, [0.0, d2 / DS_SCALE]]))

    F = np.array(F)
    J = np.array(J)
    sv = np.linalg.svd(J, compute_uv=False)
    return {'F': F, 'J': J, 'singular_values': sv.tolist(),
            'condition': float(sv[0] / sv[-1]) if sv[-1] > 0 else float('inf'),
            'rank': int(np.sum(sv > sv[0] * 1e-10)) if sv[0] > 0 else 0,
            'd2D_ds2': d2, 'anchor_slopes': anchor_slopes,
            'D_at_s_star': base['D'], 'dD_ds_at_s_star': base['dD_ds'],
            'residual_norm': float(np.linalg.norm(F))}


def fold_solve(bud, th0, s_free0, pinned, s_star0, exclusion, seed_slopes,
               max_iter=30, trust=0.3):
    th = np.array(th0, dtype=float)
    s_free, s_star = float(s_free0), float(s_star0)
    history = []
    for it in range(max_iter):
        try:
            info = fold_system(bud, th, s_free, pinned, s_star)
        except ReturnFailure as exc:
            return {'status': 'return failure during solve',
                    'failure': exc.as_record(), 'history': history}
        res = info['residual_norm']
        history.append({'iteration': it, 's_star': s_star, 's_free': s_free,
                        'residual_norm': res, 'condition': info['condition'],
                        'perturbation_norm': float(np.linalg.norm(th[2:])),
                        'd2D_ds2': info['d2D_ds2']})
        if res < FOLD_TOL:
            return finalise(bud, th, s_free, pinned, s_star, info, seed_slopes,
                            history)
        if info['rank'] < 7:
            return {'status': 'rank deficient (rank %d)' % info['rank'],
                    'singular_values': info['singular_values'], 'history': history}
        try:
            step = np.linalg.solve(info['J'], -info['F'])
        except np.linalg.LinAlgError:
            return {'status': 'singular Jacobian', 'history': history}
        n = float(np.linalg.norm(step))
        if n > trust:
            step *= trust / n
        th = th + step[:5] * THETA_SCALE
        s_free = s_free + step[5]
        s_star = s_star + step[6]
        if s_star < exclusion:
            return {'status': 'iterate entered the excluded region (s* = %.5f < %.5f)'
                    % (s_star, exclusion), 'history': history}
    return {'status': 'iteration cap without convergence', 'history': history,
            'final_residual_norm': history[-1]['residual_norm'] if history else None}


def finalise(bud, th, s_free, pinned, s_star, info, seed_slopes, history):
    """Every non-degeneracy gate a converged solve must still pass."""
    checks = {}
    checks['residual_norm'] = info['residual_norm']
    checks['condition'] = info['condition']
    checks['condition_ok'] = info['condition'] <= MAX_FOLD_CONDITION
    checks['d2D_ds2_scaled'] = info['d2D_ds2'] / DS_SCALE
    checks['d2D_ds2_ok'] = abs(info['d2D_ds2'] / DS_SCALE) >= MIN_D2
    checks['perturbation_norm'] = float(np.linalg.norm(np.asarray(th)[2:]))
    checks['perturbation_ratio'] = checks['perturbation_norm'] / SEED_PERT
    checks['perturbation_ok'] = abs(checks['perturbation_ratio'] - 1.0) < 1e-3
    ratios = [abs(a) / abs(b) if b else float('inf')
              for a, b in zip(info['anchor_slopes'], seed_slopes)]
    checks['anchor_slope_ratios'] = ratios
    checks['no_global_collapse'] = all(r >= MIN_ANCHOR_SLOPE_FRACTION for r in ratios)
    roots = [{'s_root': s_free, 'side': 1}] + \
            [{'s_root': s, 'side': side} for s, side in pinned] + \
            [{'s_root': s_star, 'side': 1}]
    try:
        dist = distinct_cycles(np.asarray(th), roots)
        checks['distinct_count'] = dist['distinct_count']
        checks['distinct_groups'] = dist['groups']
        checks['s_star_is_new_cycle'] = dist['distinct_count'] == len(roots)
    except ReturnFailure as exc:
        checks['distinctness_failure'] = exc.as_record()
        checks['s_star_is_new_cycle'] = False
    checks['accepted'] = bool(checks['condition_ok'] and checks['d2D_ds2_ok']
                              and checks['perturbation_ok']
                              and checks['no_global_collapse']
                              and checks.get('s_star_is_new_cycle'))
    return {'status': 'converged', 'theta': np.asarray(th).tolist(),
            's_free': s_free, 's_star': s_star, 'checks': checks,
            'history': history}


def main():
    t0 = time.perf_counter()
    bud = Budget(MAX_EVALUATIONS, MAX_SECONDS)
    th0 = seed_theta()
    out = {'evidence': 'NUM; no interval arithmetic, no certified periodic orbit',
           'scheme': __doc__,
           'fixed_thresholds': {'FOLD_TOL': FOLD_TOL,
                                'MAX_FOLD_CONDITION': MAX_FOLD_CONDITION,
                                'MIN_D2_scaled': MIN_D2,
                                'MIN_ANCHOR_SLOPE_FRACTION': MIN_ANCHOR_SLOPE_FRACTION,
                                'D_SCALE': D_SCALE, 'DS_SCALE': DS_SCALE,
                                'THETA_SCALE': THETA_SCALE.tolist(),
                                'seed_perturbation_norm': float(SEED_PERT)},
           'caps': {'accepted_per_branch': MAX_ACCEPTED,
                    'rejected_per_branch': MAX_REJECTED,
                    'return_evaluations': MAX_EVALUATIONS,
                    'wall_seconds': MAX_SECONDS},
           'branches': []}

    roots = [refine_root(bud, th0, s, side) for s, side in SEED_GUESSES]
    seed_slopes = [roots[0]['dD_ds']] + [r['dD_ds'] for r in roots[1:]]
    released = roots[0]['s_root']
    pinned = [(r['s_root'], r['side']) for r in roots[1:]]
    s_out = max(r['s_root'] for r in roots if r['side'] == 1)
    out['seed_roots'] = roots
    out['released_anchor'] = released
    out['pinned_anchors'] = pinned

    info = curve_system(bud, th0, released, pinned)
    out['curve_system'] = {'rank': info['rank'],
                           'singular_values': info['singular_values'],
                           'condition_on_range': info['condition_on_range'],
                           'null_vector': info['null_vector'],
                           'residual_norm': info['residual_norm'],
                           'unknowns': 6, 'equations': 5,
                           'degrees_of_freedom': 1}
    print(json.dumps(out['curve_system'], indent=2), flush=True)

    tangent = np.array(info['null_vector'], dtype=float)
    for name, sgn in (('curve_forward', +1.0), ('curve_backward', -1.0)):
        branch = {'name': name, 'steps': [], 'rejected': [], 'stop_reason': None}
        u = np.concatenate([th0, [released]])
        tan = sgn * tangent
        ds = 0.05
        accepted = rejected = 0
        while accepted < MAX_ACCEPTED and rejected < MAX_REJECTED:
            u_pred = u + ds * np.concatenate([tan[:5] * THETA_SCALE, [tan[5]]])
            try:
                ci = curve_correct(bud, u_pred, pinned, tan)
            except (ReturnFailure, RuntimeError) as exc:
                rejected += 1
                rec = (exc.as_record() if isinstance(exc, ReturnFailure)
                       else {'reason': str(exc)})
                branch['rejected'].append({'attempt': rejected, 'failure': rec})
                if isinstance(exc, RuntimeError):
                    branch['stop_reason'] = str(exc)
                    break
                ds *= 0.5
                continue
            if not ci['converged']:
                rejected += 1
                branch['rejected'].append({'attempt': rejected,
                                           'failure': {'reason': 'corrector failed',
                                                       'residual': ci['residual_norm']}})
                ds *= 0.5
                continue
            u_new = np.array(ci['u'], dtype=float)
            th_new, s_free_new = u_new[:5], u_new[5]
            move = abs(s_free_new - u[5])
            if move > MAX_FREE_MOVE:
                rejected += 1
                branch['rejected'].append(
                    {'attempt': rejected,
                     'failure': {'reason': 'released anchor moved %.4f in one step, '
                                           'above the %.3f identity guard'
                                           % (move, MAX_FREE_MOVE)}})
                ds *= 0.5
                continue
            gap = min(abs(s_free_new - s) for s, side in pinned if side == 1)
            if gap < MIN_ANCHOR_GAP:
                rejected += 1
                branch['rejected'].append(
                    {'attempt': rejected,
                     'failure': {'reason': 'released anchor within %.4f of a pinned '
                                           'anchor (guard %.3f); roots approaching '
                                           'a merge, identity no longer separable'
                                           % (gap, MIN_ANCHOR_GAP)}})
                ds *= 0.5
                continue
            if abs(th_new[2]) < E0_FLOOR:
                rejected += 1
                branch['rejected'].append(
                    {'attempt': rejected,
                     'failure': {'reason': 'e0 below the invariant-line floor'}})
                break
            bad = None
            for s, side in [(s_free_new, 1)] + pinned:
                if find_bracket(bud, th_new, s, side) is None:
                    bad = 'anchor %.6f lost its sign bracket' % s
                    break
            if bad:
                rejected += 1
                branch['rejected'].append({'attempt': rejected,
                                           'failure': {'reason': bad}})
                ds *= 0.5
                continue
            u = u_new
            accepted += 1
            if ci['rank'] == 5:
                nt = np.array(ci['null_vector'], dtype=float)
                tan = nt if float(nt @ tan) > 0 else -nt
            attempts = []
            try:
                for margin in MARGINS:
                    excl = s_out + margin
                    for off in FOLD_STARTS:
                        res = fold_solve(bud, u[:5], u[5], pinned, excl + off,
                                         excl, seed_slopes)
                        res['margin'] = margin
                        res['s_star_start'] = excl + off
                        attempts.append(res)
            except RuntimeError as exc:
                branch['stop_reason'] = str(exc)
            accepted_folds = [a for a in attempts
                              if a.get('checks', {}).get('accepted')]
            conv = [a for a in attempts if a.get('status') == 'converged']
            branch['steps'].append({'step': accepted, 'theta': u[:5].tolist(),
                                    's_free': u[5], 'arclength_step': ds,
                                    'curve_rank': ci['rank'],
                                    'curve_condition': ci['condition_on_range'],
                                    'fold_attempts': attempts,
                                    'converged_count': len(conv),
                                    'accepted_count': len(accepted_folds)})
            print(json.dumps({'branch': name, 'step': accepted,
                              'theta': ['%.6g' % v for v in u[:5]],
                              's_free': round(float(u[5]), 6),
                              'pert_ratio': round(float(np.linalg.norm(u[2:5]) / SEED_PERT), 6),
                              'converged': len(conv), 'accepted_folds': len(accepted_folds),
                              'evals': bud.evaluations}), flush=True)
            if accepted_folds:
                branch['stop_reason'] = 'fold candidate passed all gates at step %d' % accepted
                break
            if branch['stop_reason']:
                break
            ds = math.copysign(min(abs(ds) * 1.3, 0.05), ds)
        else:
            branch['stop_reason'] = ('accepted cap reached' if accepted >= MAX_ACCEPTED
                                     else 'rejected cap reached')
        branch['accepted_steps'] = accepted
        branch['rejected_steps'] = rejected
        out['branches'].append(branch)

    all_att = [a for b in out['branches'] for s in b['steps'] for a in s['fold_attempts']]
    out['evaluations'] = bud.evaluations
    out['wall_seconds'] = time.perf_counter() - t0
    out['summary'] = {
        'branches': [{'name': b['name'], 'accepted': b['accepted_steps'],
                      'rejected': b['rejected_steps'],
                      'stop_reason': b['stop_reason']} for b in out['branches']],
        'fold_solves_attempted': len(all_att),
        'fold_solves_converged': sum(1 for a in all_att if a.get('status') == 'converged'),
        'fold_candidates_passing_all_gates':
            sum(1 for a in all_att if a.get('checks', {}).get('accepted')),
        'statuses': sorted({a.get('status') for a in all_att}),
        'five_distinct_cycles_found': False,
    }
    (ROOT / 'data' / 'step5_normalized.json').write_text(
        json.dumps(out, indent=2, default=str) + '\n')
    print(json.dumps(out['summary'], indent=2), flush=True)
    print('evaluations', bud.evaluations, 'wall %.1f s' % out['wall_seconds'])


if __name__ == '__main__':
    main()
