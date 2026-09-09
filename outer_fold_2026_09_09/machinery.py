"""Repaired return-map machinery for the reversible-reseed family.

Family (five coefficients theta = (a, b, e0, e1, e2)):

    x' = (b-2)/4 + (1-b) y + a x^2 + b y^2 + e1 x + e2 x y
    y' = e0 - 2 x y

Section {x = 0}, section coordinate s = log|y|, side = sign(y).  A half
return runs from (0, side*e^s) to the next crossing of {x = 0}; the
displacement is the forward/backward log-height mismatch

    D(s; theta) = s_forward - s_backward.

WHAT THIS MODULE IS AND IS NOT
------------------------------
Everything here is double-precision numerics.  The itinerary checks below
are *numerical safeguards*: they can reject a bad return, but passing them
is not an enclosure and proves nothing about existence or isolation of a
periodic orbit.  Interval arithmetic and validated integration would be
required for that, and are not implemented.

The Cartesian and logarithmic formulations share the same DOP853 stepper
and the same event machinery.  They differ only in the coordinate the flow
is written in, so agreement between them is a COORDINATE CROSS-CHECK.  It
is not an independent numerical verification and must not be reported as
one.

FIXED RESIDUAL SCALES
---------------------
Declared once here, never rescaled from the profile being measured, so that
residuals from different runs stay comparable:

    D_SCALE      = 1e-6   typical |D| on the seed's upper section at tau=1e-4
    DS_SCALE     = 1e-6   typical |dD/ds| there, per unit s
    THETA_SCALE  = (1, 1, 1e-4, 1e-3, 1e-3)   coefficient scales at the seed

A residual divided by these is dimensionless and comparable across runs.
"""
import math

import numpy as np
from scipy.integrate import solve_ivp

TMAX = 300.0
D_SCALE = 1e-6
DS_SCALE = 1e-6
THETA_SCALE = np.array([1.0, 1.0, 1e-4, 1e-3, 1e-3])

# Numerical safeguard thresholds. These are heuristics, not enclosures.
EQ_GUARD = 1e-3          # minimum accepted distance to an equilibrium
TRANSVERSALITY = 1e-6    # minimum accepted |x'| at the return crossing
DENSE_SAMPLES = 400      # dense-output samples used by the gates


class ReturnFailure(RuntimeError):
    """A half return that did not resolve. Recorded, never read as a zero."""


def seed_theta(tau=1e-4):
    return np.array([-7.0 / 4.0, 1.0 / 3.0, -4.0 * tau / 11.0,
                     -31379.0 * tau / 25000.0, -7517.0 * tau / 5000.0])


def vector_field(th, x, y):
    a, b, e0, e1, e2 = th
    return ((b - 2.0) / 4.0 + (1.0 - b) * y + a * x * x + b * y * y + e1 * x + e2 * x * y,
            e0 - 2.0 * x * y)


def jac(th, x, y):
    a, b, e0, e1, e2 = th
    return np.array([[2.0 * a * x + e1 + e2 * y, (1.0 - b) + 2.0 * b * y + e2 * x],
                     [-2.0 * y, -2.0 * x]])


def dfdtheta(th, x, y):
    """Columns: d f / d(a, b, e0, e1, e2)."""
    return np.array([[x * x, 0.25 - y + y * y, 0.0, x, x * y],
                     [0.0, 0.0, 1.0, 0.0, 0.0]])


def equilibria(th):
    a, b, e0, e1, e2 = th
    if abs(e0) < 1e-300:
        return []
    k = e0 / 2.0            # y = k / x on {y' = 0}
    coeffs = [a, e1, (b - 2.0) / 4.0 + e2 * k, (1.0 - b) * k, b * k * k]
    pts = []
    for r in np.roots(coeffs):
        if abs(r.imag) < 1e-12 * max(1.0, abs(r.real)) and abs(r.real) > 1e-300:
            x = float(r.real)
            pts.append((x, k / x))
    return pts


def _augmented_rhs(th, sign, want_variational):
    """z = (x, y[, v(2), W(2x5)]); v = dz/ds0, W = dz/dtheta."""
    def rhs(t, u):
        x, y = u[0], u[1]
        fx, fy = vector_field(th, x, y)
        out = np.empty_like(u)
        out[0], out[1] = sign * fx, sign * fy
        if want_variational:
            J = jac(th, x, y)
            v = u[2:4]
            out[2:4] = sign * (J @ v)
            W = u[4:14].reshape(2, 5)
            out[4:14] = (sign * (J @ W + dfdtheta(th, x, y))).ravel()
        return out
    return rhs


def half_return(th, s0, side, sign, rtol=1e-12, variational=True, gates=True):
    """One half return with itinerary gates and event-corrected derivatives.

    Returns a dict with the terminal log-height, its derivatives with respect
    to s0 and to theta, the flight time, and every gate measurement.
    """
    y0 = side * math.exp(s0)
    fx0, _ = vector_field(th, 0.0, y0)
    launch = sign * fx0                      # signed departure speed off {x=0}
    if abs(launch) < TRANSVERSALITY:
        raise ReturnFailure('launch not transverse: x_dot = %.3e' % launch)
    launch_side = math.copysign(1.0, launch)

    n = 14 if variational else 2
    u0 = np.zeros(n)
    u0[1] = y0
    if variational:
        u0[3] = y0                           # d y0 / d s0 = y0

    rhs = _augmented_rhs(th, sign, variational)

    # The departure at t=0 has x' = launch; the return crossing has the
    # opposite sign, so filtering on direction excludes t=0 itself.
    ev = lambda t, u: u[0]
    ev.terminal = True
    ev.direction = -launch_side

    sol = solve_ivp(rhs, [0.0, TMAX], u0, method='DOP853', events=ev,
                    rtol=rtol, atol=rtol * 1e-3, dense_output=True)
    if not sol.success:
        raise ReturnFailure('integration failed: %s' % sol.message)
    if len(sol.t_events[0]) != 1:
        raise ReturnFailure('itinerary incomplete: %d qualifying crossings'
                            % len(sol.t_events[0]))

    T = float(sol.t_events[0][0])
    uT = sol.y_events[0][0]
    xT, yT = uT[0], uT[1]
    fxT, fyT = vector_field(th, xT, yT)
    xdotT = sign * fxT

    report = {'time': T, 'launch_xdot': launch, 'launch_side': launch_side,
              'terminal_xdot': xdotT, 'y_end': yT}

    if gates:
        ts = np.linspace(0.0, T, DENSE_SAMPLES)
        traj = sol.sol(ts)
        xs, ys = traj[0], traj[1]
        interior = xs[1:-1]
        report['interior_min_abs_x'] = float(np.min(np.abs(interior)))
        report['interior_sign_changes'] = int(np.sum(interior[:-1] * interior[1:] < 0))
        report['side_preserved'] = bool(np.all(np.sign(ys) == side))
        report['min_abs_y'] = float(np.min(np.abs(ys)))
        report['max_abs_y'] = float(np.max(np.abs(ys)))
        eq = equilibria(th)
        report['equilibria'] = eq
        report['eq_distance'] = (min(math.hypot(xs[i] - ex, ys[i] - ey)
                                     for i in range(len(xs)) for ex, ey in eq)
                                 if eq else float('inf'))
        # gate order matters: report every measurement before raising
        if report['interior_sign_changes']:
            raise ReturnFailure('unintended interior crossing of the section')
        if not report['side_preserved']:
            raise ReturnFailure('side not preserved along the return segment')
        if report['eq_distance'] < EQ_GUARD:
            raise ReturnFailure('equilibrium approach %.3e' % report['eq_distance'])
        if not np.all(np.sign(interior) == launch_side):
            raise ReturnFailure('return segment left the launch side of the section')
    if abs(xdotT) < TRANSVERSALITY:
        raise ReturnFailure('terminal crossing not transverse: %.3e' % xdotT)

    report['s_end'] = math.log(abs(yT))
    if variational:
        v = uT[2:4]
        W = uT[4:14].reshape(2, 5)
        # Event correction: T depends on s0 and theta through x(T) = 0.
        dT_ds0 = -v[0] / xdotT
        dT_dth = -W[0] / xdotT
        dy_ds0 = v[1] + sign * fyT * dT_ds0
        dy_dth = W[1] + sign * fyT * dT_dth
        report['ds_end_ds0'] = dy_ds0 / yT
        report['ds_end_dtheta'] = dy_dth / yT
        report['dT_ds0'] = dT_ds0
    return report


def displacement(th, s, side, rtol=1e-12, variational=True, gates=True):
    """D(s; theta) with event-corrected derivatives, or a recorded failure."""
    fwd = half_return(th, s, side, +1.0, rtol, variational, gates)
    bwd = half_return(th, s, side, -1.0, rtol, variational, gates)
    out = {'D': fwd['s_end'] - bwd['s_end'],
           'forward': fwd, 'backward': bwd,
           'period_estimate': fwd['time'] + bwd['time'],
           'eq_distance': min(fwd.get('eq_distance', float('inf')),
                              bwd.get('eq_distance', float('inf'))),
           'min_terminal_transversality': min(abs(fwd['terminal_xdot']),
                                              abs(bwd['terminal_xdot']))}
    if variational:
        out['dD_ds'] = fwd['ds_end_ds0'] - bwd['ds_end_ds0']
        out['dD_dtheta'] = fwd['ds_end_dtheta'] - bwd['ds_end_dtheta']
    return out


def D_only(th, s, side, rtol=1e-12):
    return displacement(th, s, side, rtol, variational=False)['D']


def derivative_audit(th, s, side, rtol_list=(1e-12, 1e-13),
                     steps=(1e-3, 1e-4, 1e-5)):
    """Variational dD/ds against central differences; reports the spread.

    The spread is an empirical uncertainty estimate, not an error bound.
    """
    rows = []
    for rtol in rtol_list:
        analytic = displacement(th, s, side, rtol)['dD_ds']
        fds = {}
        for h in steps:
            try:
                fds['h=%g' % h] = (D_only(th, s + h, side, rtol)
                                   - D_only(th, s - h, side, rtol)) / (2 * h)
            except ReturnFailure as exc:
                fds['h=%g' % h] = 'failed: %s' % exc
        rows.append({'rtol': rtol, 'variational': analytic, 'finite_differences': fds})
    numeric = [r['variational'] for r in rows]
    fd_numeric = [v for r in rows for v in r['finite_differences'].values()
                  if isinstance(v, float)]
    allv = numeric + fd_numeric
    spread = (max(allv) - min(allv)) if allv else float('nan')
    centre = float(np.mean(allv)) if allv else float('nan')
    return {'rows': rows, 'estimate': centre, 'spread': spread,
            'relative_spread': abs(spread / centre) if centre else float('inf'),
            'note': 'spread across tolerances and step sizes; '
                    'an empirical uncertainty, not an error bound'}


def log_cross_check(th, s, side, rtol=1e-12):
    """Same flow written in (x, log|y|). Shares DOP853: a coordinate check."""
    a, b, e0, e1, e2 = th

    def rhs(t, z, sign):
        x, w = z
        y = side * math.exp(w)
        fx, fy = vector_field(th, x, y)
        return (sign * fx, sign * fy / y)

    ends = []
    for sign in (+1.0, -1.0):
        fx0, _ = vector_field(th, 0.0, side * math.exp(s))
        launch_side = math.copysign(1.0, sign * fx0)
        ev = lambda t, z, sg: z[0]
        ev.terminal = True
        ev.direction = -launch_side
        sol = solve_ivp(rhs, [0.0, TMAX], [0.0, s], args=(sign,), method='DOP853',
                        events=ev, rtol=rtol, atol=rtol * 1e-3)
        if not sol.success or len(sol.t_events[0]) != 1:
            raise ReturnFailure('log-coordinate cross-check did not resolve')
        ends.append(sol.y_events[0][0][1])
    return ends[0] - ends[1]
