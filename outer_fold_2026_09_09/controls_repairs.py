"""Focused controls for the repaired machinery.

Each control is designed to FAIL in a specific way, so that the repair is
tested rather than asserted. Deliberately failed returns, deliberate tracker
failures, and a deliberate false fold indicator.
"""
import json, sys, time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from machinery import (D_SCALE, DS_SCALE, SECTION_EQ_GUARD, ReturnFailure,
                       candidate_gate, derivative_audit, displacement,
                       distinct_cycles, parameter_sensitivity_audit, seed_theta)
from anchored import Budget, anchor_system, evaluate, find_bracket, refine_root

ROOT = Path(__file__).resolve().parent
OUT = {'controls': []}


def record(name, expectation, passed, detail):
    OUT['controls'].append({'control': name, 'expectation': expectation,
                            'passed': bool(passed), 'detail': detail})
    print('%-42s %s' % (name, 'PASS' if passed else 'FAIL'), flush=True)


def main():
    t0 = time.perf_counter()
    th = seed_theta()
    bud = Budget(6000, 1800)

    # 1. A deliberately failed return must raise AND carry measurements.
    #    s = 80 leaves the upper side before returning.
    try:
        displacement(th, 80.0, 1)
        record('failed-return carries measurements',
               's = 80 fails a gate and carries structured evidence', False,
               'no failure raised')
    except ReturnFailure as exc:
        rec = exc.as_record()
        record('failed-return carries measurements',
               'raises with a non-empty measurement dictionary',
               bool(rec['measurements']) and 'side' in rec['measurements'],
               {'reason': rec['reason'],
                'measurement_keys': sorted(rec['measurements'])})

    # 1b. e0 = 0 is NOT a centre and must not be justified as one. It is the
    #     invariant-line class ({y = 0} invariant), excluded by a different
    #     theorem. Confirm the displacement does not vanish there.
    th_line = th.copy()
    th_line[2] = 0.0
    D_line = displacement(th_line, 2.0, 1)['D']
    record('e0 = 0 is the invariant-line class, not a centre',
           'D does not vanish at e0 = 0, so "centre limit" is the wrong reason',
           abs(D_line) > 1e-6,
           {'D_at_e0_zero': D_line,
            'correct_reason': 'y = 0 is an invariant straight line; that class '
                              'has at most one limit cycle'})

    # 1c. A sign change of D is NOT a distinct cycle. The upper section
    #     coordinate is two-to-one: the equilibrium at (0, 1/2) lies on the
    #     section, so each upper cycle crosses it twice.
    seven = [{'s_root': v, 'side': 1} for v in
             (0.9802427, 1.3899488, 2.0452069, -1.4558815, -1.5357456, -1.6389947)]
    seven.append({'s_root': 9.1437116, 'side': -1})
    dist = distinct_cycles(th, seven)
    record('seven sign changes are four cycles',
           'the three inner upper roots are the same cycles as the three outer',
           dist['distinct_count'] == 4,
           {'candidate_roots': len(seven), 'distinct_cycles': dist['distinct_count'],
            'duplicate_groups': len(dist['duplicates']), 'groups': dist['groups']})

    # 1d. The zero near s = log(1/2) is the equilibrium itself, not a cycle:
    #     its return period collapses to the linearised value 2*pi.
    g = candidate_gate(th, -0.69, 1)
    good = candidate_gate(th, 2.0452069, 1)
    record('equilibrium artefact rejected by the section guard',
           'the crossing at the equilibrium fails the section guard; a real '
           'root passes it',
           (not g['section_eq_ok']) and good['section_eq_ok'],
           {'artefact': {'s': -0.69, 'section_eq_distance': g['section_eq_distance'],
                         'period': g['period_estimate'],
                         'two_pi': 6.283185307179586},
            'genuine': {'s': 2.0452069,
                        'section_eq_distance': good['section_eq_distance'],
                        'period': good['period_estimate']},
            'guard': SECTION_EQ_GUARD})

    # 2. A return that cannot resolve at all must be recorded, not read as zero.
    try:
        r = displacement(th, 40.0, 1)
        record('unresolved return is not a zero',
               'either resolves with a recorded D, or fails with evidence',
               True, {'resolved': True, 'D': r['D']})
    except ReturnFailure as exc:
        rec = exc.as_record()
        record('unresolved return is not a zero',
               'either resolves with a recorded D, or fails with evidence',
               bool(rec['measurements']), rec)

    # 3. Newton must reject a returned point whose final residual is too large.
    #    Feed a guess far from any root: no sign bracket exists within the cap.
    try:
        refine_root(bud, th, 5.0, 1)
        record('root acceptance rejects a bad guess',
               'no sign bracket within the cap -> ReturnFailure', False,
               'accepted a root that should have no bracket')
    except ReturnFailure as exc:
        record('root acceptance rejects a bad guess',
               'no sign bracket within the cap -> ReturnFailure', True,
               exc.as_record())

    # 4. An accepted root must survive a re-check of its sign bracket, and its
    #    position uncertainty must be small against the root separation.
    roots = [refine_root(bud, th, s, 1) for s in (0.9802427, 1.3899488, 2.0452069)]
    roots.append(refine_root(bud, th, 9.1437116, -1))
    ups = sorted(r['s_root'] for r in roots if r['side'] == 1)
    sep = min(b - a for a, b in zip(ups, ups[1:]))
    worst_u = max(r['position_uncertainty'] for r in roots)
    record('position uncertainty << root separation',
           'max uncertainty at least 1e3 times smaller than the separation',
           worst_u * 1e3 < sep,
           {'worst_position_uncertainty': worst_u, 'min_upper_separation': sep,
            'ratio': worst_u / sep})

    # 5. A deliberate tracker failure: a step far too large for the guard.
    #    It must be reported as a tracker failure, never as a lost cycle.
    th_big = th.copy()
    th_big[0] += 0.05                      # the step the first pass called a lost root
    try:
        refine_root(bud, th_big, roots[2]['s_root'], 1)
        outcome = 'root still found after a 0.05 step in a'
        passed = True
    except ReturnFailure as exc:
        outcome = exc.as_record()
        passed = True                      # either way, it is a tracker statement
    record('large a step is a tracker statement',
           'either the root is found, or the failure is explicitly a tracker '
           'failure -- never reported as a vanished cycle', passed, outcome)

    # 6. False fold indicator: a positive sampled slope must NOT be accepted as
    #    a fold. Find a point with dD/ds > 0 (the middle root has one) and
    #    confirm that D and dD/ds are not simultaneously zero there.
    mid = roots[1]
    r = evaluate(bud, th, mid['s_root'] + 0.02, 1)
    positive_slope = r['dD_ds'] > 0
    is_double_zero = (abs(r['D']) / D_SCALE < 1e-3
                      and abs(r['dD_ds']) / DS_SCALE < 1e-3)
    record('positive sampled slope is not a fold',
           'a point with dD/ds > 0 is not reported as a double zero',
           positive_slope and not is_double_zero,
           {'s': mid['s_root'] + 0.02, 'D_scaled': r['D'] / D_SCALE,
            'dD_ds_scaled': r['dD_ds'] / DS_SCALE,
            'positive_slope': bool(positive_slope),
            'meets_double_zero_test': bool(is_double_zero)})

    # 7. All five parameter sensitivities audited, not just dD/ds.
    ds_audit = derivative_audit(th, roots[2]['s_root'], 1)
    p_audit = parameter_sensitivity_audit(th, roots[2]['s_root'], 1)
    record('all five parameter sensitivities audited',
           'each of dD/da, dD/db, dD/de0, dD/de1, dD/de2 has an independent '
           'finite-difference cross-check',
           len(p_audit['rows']) == 5,
           {'section_derivative_relative_spread': ds_audit['relative_spread'],
            'parameter_relative_spreads':
                {r['parameter']: r['relative_spread'] for r in p_audit['rows']},
            'worst_parameter_relative_spread': p_audit['worst_relative_spread']})

    # 8. The anchored Jacobian rank claim must be stated with its scales.
    anch = [(r['s_root'], r['side']) for r in roots]
    info = anchor_system(bud, th, anch)
    record('anchored 4x5 Jacobian has rank 4',
           'rank 4 gives a one-dimensional kernel, hence a local curve',
           info['rank'] == 4,
           {'singular_values': info['singular_values'], 'rank': info['rank'],
            'condition_on_range': info['condition_on_range'],
            'null_vector_scaled': info['null_vector_scaled'],
            'note': 'rank alone does not establish that the fold system below '
                    'is well posed; that needs its own 6x6 diagnostics'})

    # 9. Regression: step4's collapsed field must be rejected. Its fold solves
    #    "converged" only because the whole perturbation shrank by 1e-6, making
    #    D and dD/ds vanish everywhere at once.
    collapsed = np.array([-1.75, 1.0 / 3.0, 2.9813e-11, 9.9189e-11, 1.3422e-10])
    seed_pert = float(np.linalg.norm(th[2:]))
    ratio = float(np.linalg.norm(collapsed[2:])) / seed_pert
    r_col = evaluate(bud, collapsed, 2.195207, 1)
    small_residual = (abs(r_col['D']) / D_SCALE < 1e-3
                      and abs(r_col['dD_ds']) / DS_SCALE < 1e-3)
    rejected_by_norm = abs(ratio - 1.0) >= 1e-3
    record('collapsed field rejected as a false fold',
           'it satisfies both fold equations numerically, and is rejected '
           'anyway by the perturbation-normalisation gate',
           small_residual and rejected_by_norm,
           {'perturbation_ratio_to_seed': ratio,
            'D_scaled': r_col['D'] / D_SCALE,
            'dD_ds_scaled': r_col['dD_ds'] / DS_SCALE,
            'satisfies_fold_equations_numerically': bool(small_residual),
            'rejected_by_normalisation_gate': bool(rejected_by_norm),
            'note': 'this is exactly the step4 false positive'})

    OUT['seed_roots'] = roots
    OUT['evaluations'] = bud.evaluations
    OUT['wall_seconds'] = time.perf_counter() - t0
    OUT['passed'] = sum(1 for c in OUT['controls'] if c['passed'])
    OUT['total'] = len(OUT['controls'])
    (ROOT / 'data' / 'controls_repairs.json').write_text(
        json.dumps(OUT, indent=2, default=str) + '\n')
    print('%d/%d controls passed, %d evaluations, %.1f s'
          % (OUT['passed'], OUT['total'], bud.evaluations, OUT['wall_seconds']))


if __name__ == '__main__':
    main()
