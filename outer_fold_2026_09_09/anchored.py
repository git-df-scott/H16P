"""Anchor-preserving continuation and an actual double-zero solve.

This replaces the slope-maximisation experiment, which was wrong in ways
recorded in OUTER_FOLD_RESULT_2026_09_09.md.

Scheme
------
Fix four section anchors s_1..s_4 at refined seed roots, with the correct
side for each, and treat

    F_i(theta) = D(s_i; theta, side_i),    i = 1..4

as four equations in the five coefficients theta = (a, b, e0, e1, e2).  The
seed is a solution.  If the scaled 4x5 Jacobian has rank 4 its kernel is one
dimensional, so the solution set is locally a curve; we follow it by
pseudo-arclength predictor-corrector.  Every point of that curve has all
four cycles pinned at exactly the same four section coordinates, so no
tracker identity can drift -- the anchors cannot swap with each other.

Beyond the outermost upper anchor we then solve the genuine fold system

    F_1..F_4 = 0,   D(s*; theta) = 0,   dD/ds (s*; theta) = 0

-- six equations in six unknowns (theta, s*) -- with a damped Newton
corrector.  A fold candidate requires a converged solve with small residual
in the FIXED scales, a nonzero second section derivative, and a nonzero
transverse unfolding derivative.  A positive sampled slope is not a
criterion and is not used as one; a double zero need not change sign.

Evidence class: NUM throughout.  Nothing here certifies a periodic orbit.
"""
import math

import numpy as np

from machinery import (D_SCALE, DS_SCALE, THETA_SCALE, ReturnFailure,
                       displacement, D_only)

ACCEPT_D = 5e-14          # accepted residual at an anchor, absolute
BRACKET_MAX = 0.2         # widest half-width searched for a sign bracket


class Budget:
    """Explicit caps. Exceeding one raises, and the run records why."""

    def __init__(self, evaluations, seconds):
        import time
        self.max_evaluations = evaluations
        self.max_seconds = seconds
        self.evaluations = 0
        self._t0 = time.perf_counter()
        self._time = time.perf_counter

    def spend(self, n=1):
        self.evaluations += n
        if self.evaluations > self.max_evaluations:
            raise RuntimeError('return-evaluation cap %d reached'
                               % self.max_evaluations)
        if self._time() - self._t0 > self.max_seconds:
            raise RuntimeError('wall-time cap %.0f s reached' % self.max_seconds)

    def elapsed(self):
        return self._time() - self._t0


def evaluate(bud, th, s, side, rtol=1e-12, variational=True):
    bud.spend()
    return displacement(np.asarray(th, dtype=float), s, side, rtol,
                        variational=variational)


def find_bracket(bud, th, s, side, w0=1e-3):
    """Smallest sampled half-width whose endpoints straddle zero."""
    w = w0
    while w <= BRACKET_MAX:
        lo = evaluate(bud, th, s - w, side, variational=False)['D']
        hi = evaluate(bud, th, s + w, side, variational=False)['D']
        if lo * hi < 0:
            return {'half_width': w, 'D_low': lo, 'D_high': hi}
        w *= 2
    return None


def refine_root(bud, th, s_guess, side):
    """Bracket-preserving safeguarded Newton with a FINAL acceptance check.

    The previous implementation accepted on a pre-update residual and never
    checked the value it returned. This one keeps a sign bracket throughout,
    bisects whenever Newton would leave it, and accepts only after
    re-evaluating at the returned point.
    """
    br = find_bracket(bud, th, s_guess, side)
    if br is None:
        raise ReturnFailure('no sign bracket within +-%g of the guess' % BRACKET_MAX,
                            {'s_guess': s_guess, 'side': side})
    w = br['half_width']
    lo, hi = s_guess - w, s_guess + w
    D_lo = br['D_low']
    s = s_guess
    for _ in range(60):
        r = evaluate(bud, th, s, side)
        D, dD = r['D'], r['dD_ds']
        if D * D_lo > 0:
            lo = s
            D_lo = D
        else:
            hi = s
        if abs(dD) > 1e-16:
            s_new = s - D / dD
        else:
            s_new = 0.5 * (lo + hi)
        if not (lo < s_new < hi):
            s_new = 0.5 * (lo + hi)
        if abs(s_new - s) < 1e-14 or abs(D) < ACCEPT_D * 0.2:
            s = s_new
            break
        s = s_new
    final = evaluate(bud, th, s, side)
    if not abs(final['D']) <= ACCEPT_D:
        raise ReturnFailure('final residual %.3e exceeds the acceptance level %.1e'
                            % (final['D'], ACCEPT_D),
                            {'s': s, 'D': final['D'], 'dD_ds': final['dD_ds']})
    check = find_bracket(bud, th, s, side)
    if check is None:
        raise ReturnFailure('accepted root has no surviving sign bracket',
                            {'s': s, 'D': final['D']})
    uncertainty = max(abs(final['D']), ACCEPT_D) / max(abs(final['dD_ds']), 1e-300)
    return {'s_root': s, 'side': side, 'D': final['D'],
            'D_scaled': final['D'] / D_SCALE, 'dD_ds': final['dD_ds'],
            'dD_dtheta': final['dD_dtheta'].tolist(),
            'bracket_half_width': check['half_width'],
            'bracket_D_low': check['D_low'], 'bracket_D_high': check['D_high'],
            'position_uncertainty': uncertainty,
            'period_estimate': final['period_estimate'],
            'eq_distance': final['eq_distance'],
            'terminal_transversality': final['min_terminal_transversality']}


def anchor_system(bud, th, anchors):
    """F (scaled) and its scaled 4x5 Jacobian, with SVD diagnostics."""
    F, J = [], []
    rows = []
    for s, side in anchors:
        r = evaluate(bud, th, s, side)
        F.append(r['D'] / D_SCALE)
        J.append(r['dD_dtheta'] * THETA_SCALE / D_SCALE)
        rows.append({'s': s, 'side': side, 'D': r['D'], 'D_scaled': r['D'] / D_SCALE,
                     'dD_ds': r['dD_ds'], 'eq_distance': r['eq_distance']})
    F = np.array(F)
    J = np.array(J)
    U, sv, Vt = np.linalg.svd(J)
    rank = int(np.sum(sv > sv[0] * 1e-10)) if sv[0] > 0 else 0
    return {'F': F, 'J': J, 'rows': rows, 'singular_values': sv.tolist(),
            'rank': rank,
            'condition_on_range': float(sv[0] / sv[rank - 1]) if rank else float('inf'),
            'null_vector_scaled': Vt[-1].tolist(),
            'residual_norm_scaled': float(np.linalg.norm(F))}


def correct(bud, th, anchors, tangent, ds, max_iter=12, trust=0.5):
    """Pseudo-arclength corrector: 4 anchor equations + 1 arclength condition."""
    th = np.array(th, dtype=float)
    th_pred = th.copy()
    for _ in range(max_iter):
        info = anchor_system(bud, th, anchors)
        arc = float(((th - th_pred) / THETA_SCALE) @ tangent)
        G = np.concatenate([info['F'], [arc]])
        if np.linalg.norm(G) < 1e-6:
            info['converged'] = True
            info['theta'] = th.tolist()
            info['final_residual_norm'] = float(np.linalg.norm(G))
            return info
        A = np.vstack([info['J'], tangent])
        try:
            delta_scaled = np.linalg.solve(A, -G)
        except np.linalg.LinAlgError:
            raise ReturnFailure('corrector Jacobian singular',
                                {'singular_values': info['singular_values']})
        n = np.linalg.norm(delta_scaled)
        if n > trust:
            delta_scaled = delta_scaled * (trust / n)
        th = th + delta_scaled * THETA_SCALE
    info = anchor_system(bud, th, anchors)
    info['converged'] = False
    info['theta'] = th.tolist()
    info['final_residual_norm'] = float(np.linalg.norm(info['F']))
    return info


def fold_system(bud, th, anchors, s_star, rel=1e-3):
    """G and its scaled 6x6 Jacobian for the genuine double-zero system."""
    th = np.array(th, dtype=float)
    F, JF = [], []
    for s, side in anchors:
        r = evaluate(bud, th, s, side)
        F.append(r['D'] / D_SCALE)
        JF.append(np.concatenate([r['dD_dtheta'] * THETA_SCALE / D_SCALE, [0.0]]))
    base = evaluate(bud, th, s_star, 1)
    hs = 1e-5
    up = evaluate(bud, th, s_star + hs, 1)
    dn = evaluate(bud, th, s_star - hs, 1)
    d2D = (up['dD_ds'] - dn['dD_ds']) / (2 * hs)

    dslope_dtheta = np.zeros(5)
    for j in range(5):
        h = rel * THETA_SCALE[j]
        tp, tm = th.copy(), th.copy()
        tp[j] += h
        tm[j] -= h
        dslope_dtheta[j] = (evaluate(bud, tp, s_star, 1)['dD_ds']
                            - evaluate(bud, tm, s_star, 1)['dD_ds']) / (2 * h)

    row5 = np.concatenate([base['dD_dtheta'] * THETA_SCALE / D_SCALE,
                           [base['dD_ds'] / D_SCALE]])
    row6 = np.concatenate([dslope_dtheta * THETA_SCALE / DS_SCALE,
                           [d2D / DS_SCALE]])
    G = np.array(F + [base['D'] / D_SCALE, base['dD_ds'] / DS_SCALE])
    J = np.vstack(JF + [row5, row6])
    sv = np.linalg.svd(J, compute_uv=False)
    return {'G': G, 'J': J, 'singular_values': sv.tolist(),
            'condition': float(sv[0] / sv[-1]) if sv[-1] > 0 else float('inf'),
            'rank': int(np.sum(sv > sv[0] * 1e-10)) if sv[0] > 0 else 0,
            'd2D_ds2': d2D, 'D_at_s_star': base['D'],
            'dD_ds_at_s_star': base['dD_ds'],
            'residual_norm_scaled': float(np.linalg.norm(G))}


def fold_solve(bud, th0, anchors, s_star0, exclusion, max_iter=25, trust=0.4,
               tol=1e-3):
    """Damped Newton on the six-equation fold system. Bounded and reported."""
    th = np.array(th0, dtype=float)
    s_star = float(s_star0)
    history = []
    for it in range(max_iter):
        try:
            info = fold_system(bud, th, anchors, s_star)
        except ReturnFailure as exc:
            return {'status': 'return failure during solve',
                    'failure': exc.as_record(), 'iterations': it,
                    'history': history}
        res = float(np.linalg.norm(info['G']))
        history.append({'iteration': it, 's_star': s_star,
                        'residual_norm_scaled': res,
                        'condition': info['condition'], 'rank': info['rank'],
                        'D_scaled': info['D_at_s_star'] / D_SCALE,
                        'dD_ds_scaled': info['dD_ds_at_s_star'] / DS_SCALE})
        if res < tol:
            return {'status': 'converged', 'theta': th.tolist(), 's_star': s_star,
                    'residual_norm_scaled': res, 'condition': info['condition'],
                    'rank': info['rank'], 'd2D_ds2': info['d2D_ds2'],
                    'iterations': it, 'history': history}
        if info['rank'] < 6:
            return {'status': 'rank deficient at iteration %d (rank %d)'
                    % (it, info['rank']), 'singular_values': info['singular_values'],
                    'iterations': it, 'history': history}
        try:
            step = np.linalg.solve(info['J'], -info['G'])
        except np.linalg.LinAlgError:
            return {'status': 'singular Jacobian', 'iterations': it,
                    'history': history}
        n = float(np.linalg.norm(step))
        if n > trust:
            step = step * (trust / n)
        th = th + step[:5] * THETA_SCALE
        s_star = s_star + step[5]
        if s_star < exclusion:
            return {'status': 'iterate entered the excluded region near the '
                            'outermost anchor (s* = %.6f < %.6f)' % (s_star, exclusion),
                    'iterations': it, 'history': history}
    return {'status': 'iteration cap reached without convergence',
            'iterations': max_iter, 'history': history,
            'final_residual_norm_scaled': history[-1]['residual_norm_scaled']
            if history else None}
