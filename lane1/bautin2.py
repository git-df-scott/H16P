"""The cusp manifold, calibrated numerically -- no printed focal values.

`bautin.py` imposes the triple-root ratios on Cherkas--Artes--Llibre's PRINTED
focal values.  That is wrong: the printed `V1, V3, V5, V7` carry per-field
normalisation constants (measured `c7/V7` of `-6.7e-05` and `-9.0e-06`,
differing in magnitude and in sign), so the ratios do not transfer.

Here the true Taylor coefficients of the displacement are measured instead:

        D(s) = c1 s + c3 s^3 + c5 s^5 + c7 s^7 + O(s^9)

fitted in binary128 under the two-tolerance noise gate, and the cusp conditions

        c1 = -r0^6 c7,   c3 = 3 r0^4 c7,   c5 = -3 r0^2 c7

are solved on those.  No normalisation constant can survive that.

The continuation has an exact starting point.  At `r0 = 0` the three conditions
collapse to `c1 = c3 = c5 = 0`, which is precisely a third-order weak focus --
and Cherkas et al. hand us a closed-form two-parameter family of those
(`a11 = 4 - 2a`, `a01 = 2a + 1 - a11`,
`a10 = (6(a^2 - a - 2) + a20(6a - 7))/(1 - 3a)`), verified against this engine
to slope 7.10.  So the `r0 = 0` slice of the cusp manifold is known in closed
form, and growing `r0` from there is exactly "continue the cusp manifold to
normal amplitude", which PROTOCOL section (c) says nobody has done.
"""
import math
import numpy as np
import engine as E
import bautin as B


# ---------------------------------------------------------------- the fit
def taylor_D(L, phi, s_lo, s_hi, npts=16, rtol=1e-18, rtol2=1e-16, gate=30.0,
             deg=5):
    """(c1, c3, c5, c7) of D(s) near the focus, or None.

    D/s is a cubic in w = s^2, so the fit is a cubic least squares in w with w
    rescaled to [-1,1] for conditioning.  Points whose |D| is not `gate` times
    its own two-tolerance noise estimate are dropped: at these amplitudes |D|
    runs down to 1e-28 and the binary128 floor is around 1e-24, so without the
    gate the fit reads noise and reports the wrong order of the weak focus.
    """
    s = np.ascontiguousarray(np.geomspace(s_lo, s_hi, npts))
    D1, st1, _ = E.d_curve_quad(L, phi, s, rtol=rtol)
    D2, st2, _ = E.d_curve_quad(L, phi, s, rtol=rtol2)
    noise = 10.0 * np.abs(D1 - D2) + 1e-300
    ok = (st1 == 0) & (st2 == 0) & (np.abs(D1) > gate * noise)
    if ok.sum() < deg + 4:
        return None
    ss, DD = s[ok], D1[ok]
    w = ss ** 2
    lo, hi = w.min(), w.max()
    z = 2.0 * (w - lo) / (hi - lo) - 1.0
    # RELATIVE least squares.  D/s spans seven decades across the window, so an
    # unweighted fit is decided entirely by the largest points at the top of it,
    # where the O(s^9) tail lives -- c5 then absorbs that tail and comes out
    # wrong by orders of magnitude while the residual still looks small.  The
    # small-s points are exactly the ones that pin c5, so they must carry equal
    # weight: minimise (model/y - 1) instead of (model - y).
    y = DD / ss
    Wt = 1.0 / np.abs(y)
    # Fit DEGREE 5 in w, not 3: the O(s^9) and O(s^11) tail is what corrupted c5
    # before.  Giving the tail its own coefficients costs nothing and stops c5
    # from absorbing it; only c1, c3, c5, c7 are read off afterwards.
    Vd = np.vander(z, deg + 1)                     # highest power first
    cz, *_ = np.linalg.lstsq(Vd * Wt[:, None], y * Wt, rcond=None)
    # z = A w + Bq  ->  express as a polynomial in w
    A = 2.0 / (hi - lo)
    Bq = -2.0 * lo / (hi - lo) - 1.0
    p = np.polynomial.polynomial.Polynomial([Bq, A])
    q = sum(cz[::-1][k] * p ** k for k in range(deg + 1))
    c = q.convert().coef
    if len(c) < 4:
        c = np.concatenate([c, np.zeros(4 - len(c))])
    resid = float(np.max(np.abs(np.polyval(cz, z) / y - 1.0)))
    return dict(c1=float(c[0]), c3=float(c[1]), c5=float(c[2]), c7=float(c[3]),
                n_used=int(ok.sum()), rel_resid=resid,
                s_lo=float(ss.min()), s_hi=float(ss.max()))


def fit_window(L, phi, s_ref, gate=30.0):
    """Pick a fit window that is above the noise and still inside the O(s^9)
    truncation: try a few and keep the one with the smallest cubic residual.

    binary128 cost is steeply nonlinear in the tolerance -- 0.12 s per return at
    rtol 1e-18 against 4.4 s at 1e-20 -- so the pair is 1e-18 / 1e-16 and the
    window list is short.  It is the noise gate, not the tolerance, that keeps
    the fit honest.

    Calibrated on the verified third-order family (a = -2, a20 = -1), where the
    direct measurement gives c7 = 4.0758 and c5 = 0 exactly:

        window                 deg  rel_resid   c5         c7
        [5.1e-3, 2.0e-2]        5    2.8e-05   -2.8e-06   4.1286
        [4.7e-3, 1.3e-2]        5    3.3e-06   -1.2e-06   4.1063
        [4.5e-3, 2.5e-2]        3    2.8e-02   -5.3e-05   4.4447

    A narrow window with the tail fitted gets c7 to 0.7% and drives the spurious
    c5 down to 1.4% of the c7 term at the bottom of the window.  A wide window
    with a bare cubic is off by 9% in c7 and by a factor of 40 in c5, which is
    what made the first pass at this useless.
    """
    best = None
    for (lo_mult, hi_mult) in ((1e-2, 2.5e-2), (1e-2, 4e-2), (2e-2, 6e-2),
                               (5e-3, 2e-2), (3e-2, 1e-1)):
        f = taylor_D(L, phi, lo_mult * s_ref, hi_mult * s_ref, gate=gate)
        if f is None:
            continue
        if best is None or f["rel_resid"] < best["rel_resid"]:
            best = f
    return best


# ------------------------------------------------------------ the solve
def params_to_L(a, a20, a11, a01, a10):
    v = B.vec12(dict(a=a, a20=a20, a11=a11, a01=a01, a10=a10))
    J = E.jac(v, (1.0, -1.0))
    tr, dt = J[0, 0] + J[1, 1], float(np.linalg.det(J))
    if not (dt > 0 and tr * tr < 4 * dt):
        return None, v, tr, dt
    return E.local10(v, (1.0, -1.0)), v, tr, dt


def residual(a, a20, x, r0, phi=0.0, s_ref=None, window=None):
    """F(x) for x = (a11, a01, a10), scaled by c7 so it is dimensionless.

    `window` MUST be held fixed along a Newton solve.  `fit_window` chooses
    among a discrete list by residual, so letting it re-choose makes the fitted
    c5 a discontinuous function of x, and no finite-difference Jacobian can
    survive that -- the symptom is a Newton that reports no descent direction
    however far it backtracks.
    """
    a11, a01, a10 = x
    L, v, tr, dt = params_to_L(a, a20, a11, a01, a10)
    if L is None:
        return None, None
    if s_ref is None:
        import ahsweep as A
        s_ref = A.length_scale(L)
    if window is None:
        f = fit_window(L, phi, s_ref)
    else:
        f = taylor_D(L, phi, window[0], window[1])
    if f is None or f["c7"] == 0.0 or not fit_is_sane(f, L, phi):
        return None, None
    c7 = f["c7"]
    # Scale each condition by the size of its OWN two terms, never by |c7|.
    # Dividing by |c7| makes every component small whenever the fit blows c7 up,
    # so a solver can "descend" by driving the fit degenerate instead of by
    # finding a cusp: a run that reported F = (4e-15, -8e-10, +4.8e-05) ->
    # (-8e-11, 6e-09, -4.4e-05), an apparent sign change bracketing a cusp, was
    # at a point where the fit had returned c7 = 1.8e+08 and c5 = -1.6e+04.
    # Scaled this way each component lies in [-1, 1] and is small only when the
    # two terms genuinely cancel.
    def rel(u, v):
        d = abs(u) + abs(v)
        return 0.0 if d == 0.0 else (u + v) / d
    F = np.array([rel(f["c1"], r0 ** 6 * c7),
                  rel(f["c3"], -3.0 * r0 ** 4 * c7),
                  rel(f["c5"], 3.0 * r0 ** 2 * c7)])
    return F, f


def fit_is_sane(f, L, phi, tol_resid=1e-3, dom=0.3):
    """Reject a fit whose s^7 term is not actually visible in the window.

    Away from a weak focus, D is dominated by c1 s across the whole usable
    window and c3, c5, c7 are unconstrained: a degree-5 least squares will
    happily return c7 = 1.8e+08.  Require the s^7 term to carry at least `dom`
    of |D| at the top of the window, and the fit residual to be small.
    """
    if f is None or f["rel_resid"] > tol_resid:
        return False
    s = f["s_hi"]
    terms = [abs(f["c1"]) * s, abs(f["c3"]) * s ** 3,
             abs(f["c5"]) * s ** 5, abs(f["c7"]) * s ** 7]
    tot = sum(terms)
    return tot > 0 and terms[3] / tot >= dom


def third_order_seed(a, a20):
    """The exact r0 = 0 point of the cusp manifold."""
    if abs(1 - 3 * a) < 1e-12 or abs(a - 2) < 1e-12:
        return None
    a11 = 4 - 2 * a
    a01 = 2 * a + 1 - a11
    a10 = (6 * (a * a - a - 2) + a20 * (6 * a - 7)) / (1 - 3 * a)
    return np.array([a11, a01, a10], float)


def newton(a, a20, x0, r0, phi=0.0, iters=12, tol=1e-6, hstep=2e-3,
           cap_rel=0.02, backtrack=6, window=None):
    """Damped Newton on the three calibrated cusp conditions in (a11, a01, a10).

    Step control matters more than usual here because the Jacobian is finite
    differenced from a FITTED quantity: the coefficients carry roughly 1e-6 of
    relative noise, so a differencing step of 1e-5 gives derivative estimates
    with order-one noise and the first "Newton" step is then arbitrary.  A
    1e-5 step moved a11 from 8.0 to 5.4 in one go and left the domain.  Hence
    `hstep = 2e-3` (well above the fit noise, still inside the linear regime),
    a step capped at 2% of ||x||, and backtracking that only accepts a decrease
    in ||F||.
    """
    x = np.array(x0, float)
    F, f = residual(a, a20, x, r0, phi, window=window)
    if F is None:
        return dict(ok=False, why="unresolved_start", x=x.tolist(), it=0)
    for it in range(iters):
        if np.max(np.abs(F)) < tol:
            return dict(ok=True, x=x.tolist(), F=F.tolist(), fit=f, it=it)
        J = np.zeros((3, 3))
        for j in range(3):
            xp = x.copy()
            dx = hstep * max(1.0, abs(x[j]))
            xp[j] += dx
            Fp, _ = residual(a, a20, xp, r0, phi, window=window)
            if Fp is None:
                return dict(ok=False, why="jacobian", x=x.tolist(), it=it,
                            F=F.tolist())
            J[:, j] = (Fp - F) / dx
        try:
            step = np.linalg.solve(J, -F)
        except np.linalg.LinAlgError:
            return dict(ok=False, why="singular", x=x.tolist(), it=it,
                        F=F.tolist())
        nrm = np.linalg.norm(step)
        cap = cap_rel * max(1.0, np.linalg.norm(x))
        if nrm > cap:
            step *= cap / nrm
        n0 = np.linalg.norm(F)
        moved = False
        for _ in range(backtrack):
            xt = x + step
            Ft, ft = residual(a, a20, xt, r0, phi, window=window)
            if Ft is not None and np.linalg.norm(Ft) < n0:
                x, F, f, moved = xt, Ft, ft, True
                break
            step *= 0.4
        if not moved:
            return dict(ok=False, why="no_descent", x=x.tolist(), it=it,
                        F=F.tolist())
    return dict(ok=False, why="no_convergence", x=x.tolist(),
                F=F.tolist(), it=iters)


def bracket_solve(a, a20, xa, xb, r0, window, phi=0.0, iters=40, tol=1e-7):
    """Land the cusp by bisecting the governing residual along a segment.

    Newton brackets the cusp and then stalls, because its Jacobian is finite
    differenced from a fitted quantity and stops being reliable once |F3| is
    down at the fit's own noise on c5.  But F1 and F2 are already at 1e-11 and
    6e-9 by then and the whole residual is F3, which is smooth and monotone
    along the path -- so a one-dimensional bisection lands it with no Jacobian
    at all.  `xa` and `xb` must give F3 of opposite signs.
    """
    xa, xb = np.array(xa, float), np.array(xb, float)
    Fa, _ = residual(a, a20, xa, r0, phi, window=window)
    Fb, _ = residual(a, a20, xb, r0, phi, window=window)
    if Fa is None or Fb is None:
        return dict(ok=False, why="endpoint_unresolved")
    if Fa[2] * Fb[2] >= 0:
        return dict(ok=False, why="not_bracketed",
                    F3=[float(Fa[2]), float(Fb[2])])
    for it in range(iters):
        xm = 0.5 * (xa + xb)
        Fm, fm = residual(a, a20, xm, r0, phi, window=window)
        if Fm is None:
            return dict(ok=False, why="unresolved", x=xm.tolist(), it=it)
        if abs(Fm[2]) < tol:
            return dict(ok=True, x=xm.tolist(), F=Fm.tolist(), fit=fm, it=it)
        if Fa[2] * Fm[2] < 0:
            xb, Fb = xm, Fm
        else:
            xa, Fa = xm, Fm
        if np.linalg.norm(xb - xa) < 1e-13 * (1 + np.linalg.norm(xa)):
            break
    xm = 0.5 * (xa + xb)
    Fm, fm = residual(a, a20, xm, r0, phi, window=window)
    return dict(ok=(Fm is not None and abs(Fm[2]) < 10 * tol),
                why="tolerance_floor", x=xm.tolist(),
                F=(Fm.tolist() if Fm is not None else None), fit=fm, it=iters)
