"""Second, independent integrator (PROTOCOL rule 2).

scipy DOP853 in GLOBAL (un-translated) coordinates with root-finding events on
the section line.  This differs from the lane engine in method (scipy's
Runge-Kutta 8(5,3) with its own step controller vs a hand-rolled DP5(4)), in
coordinates (global vs focus-local), and in the crossing predicate (two
terminal half-turn events vs an accumulated-angle target).
"""
import numpy as np
from scipy.integrate import solve_ivp

def _rhs_factory(c):
    c = np.asarray(c, float)
    def rhs(t, z):
        x, y = z
        return [c[0] + c[1]*x + c[2]*y + c[3]*x*x + c[4]*x*y + c[5]*y*y,
                c[6] + c[7]*x + c[8]*y + c[9]*x*x + c[10]*x*y + c[11]*y*y]
    return rhs

def rotate(c, b):
    c = np.asarray(c, float)
    cb, sb = np.cos(b), np.sin(b)
    out = np.empty(12)
    out[:6] = c[:6]*cb - c[6:]*sb
    out[6:] = c[:6]*sb + c[6:]*cb
    return out

def _leg(rhs, z, t0, Tmax, dirn, fx, fy, nx, ny, rtol, atol, Rmax):
    def ev(t, zz):
        return (zz[0]-fx)*nx + (zz[1]-fy)*ny
    ev.terminal = True
    ev.direction = dirn
    def esc(t, zz):
        return Rmax - np.hypot(zz[0]-fx, zz[1]-fy)
    esc.terminal = True
    esc.direction = -1.0
    sol = solve_ivp(rhs, (t0, t0+Tmax), z, method="DOP853",
                    rtol=rtol, atol=atol, events=[ev, esc])
    if len(sol.t_events[0]) == 0:
        return None, None
    return float(sol.t_events[0][0]), np.asarray(sol.y_events[0][0], float)

def _ret_radius(cglob, focus, direction, s, b, rtol, atol, Tmax, Rmax):
    c = rotate(cglob, b)
    fx, fy = float(focus[0]), float(focus[1])
    ux, uy = direction
    nx, ny = -uy, ux
    rhs = _rhs_factory(c)
    z = np.array([fx + s*ux, fy + s*uy], float)
    v = rhs(0.0, z)
    sp = v[0]*nx + v[1]*ny
    if sp == 0.0:
        return None
    sense = 1.0 if sp > 0 else -1.0
    t0 = 0.0
    at = atol*max(1.0, s)
    for leg, dirn in ((0, -sense), (1, sense)):
        t0, z = _leg(rhs, z, t0, Tmax, dirn, fx, fy, nx, ny, rtol, at, Rmax)
        if z is None:
            return None
        dot = (z[0]-fx)*ux + (z[1]-fy)*uy
        want = -1.0 if leg == 0 else 1.0
        if np.sign(dot) != want:
            return None
    return float(np.hypot(z[0]-fx, z[1]-fy))

def ret_radius(cglob, focus, direction, s, b=0.0, rtol=1e-12, atol=1e-16,
               Tmax=None, Rmax=1e6):
    """Radius at the first return of focus + s*direction to the ray, or None.

    scipy's cost grows with the nominal t_span even when a terminal event fires
    early, so the time bound is escalated rather than set large up front.
    """
    if Tmax is not None:
        return _ret_radius(cglob, focus, direction, s, b, rtol, atol, Tmax, Rmax)
    for tm in (60.0, 1200.0, 30000.0):
        R = _ret_radius(cglob, focus, direction, s, b, rtol, atol, tm, Rmax)
        if R is not None:
            return R
    return None

def displacement(cglob, focus, direction, svals, b=0.0, **kw):
    out = []
    for s in svals:
        R = ret_radius(cglob, focus, direction, float(s), b=b, **kw)
        out.append(np.nan if R is None else R - float(s))
    return np.array(out)
