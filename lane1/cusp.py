"""Direct Newton continuation onto the CUSP of the Andronov-Hopf curve.

A hill climb on "how close is beta* to growing a third extremum" stalls: the
event is a codimension-one variety in the 8-dimensional live coefficient space,
and a (1+1)-ES walks past it.  It is far better to solve for it.

Write u = log s and A(u, lambda) = beta*(s; L(lambda)) for the base field moved
along one direction v,
        L(lambda) = ||L0|| * (L0 + lambda ||L0|| v) / ||L0 + lambda ||L0|| v||.
Interior extrema of beta* are the zeros of A_u; they are folds of the family
(multiplicity-two limit cycles).  A CUSP -- a limit cycle of multiplicity three,
PROTOCOL section (c) -- is

        A_u(u, lambda) = 0,   A_uu(u, lambda) = 0,

two equations in the two unknowns (u, lambda): an isolated point, reachable by
Newton.  Crossing lambda through it creates or destroys a PAIR of extrema of
beta*.  Started in a region OUTSIDE the interval already spanned by the seed's
two extrema, a solution is a cusp that adds two extrema -- four in all, i.e. a
horizontal line meeting beta* five times.  Started between them it is the cusp
where the seed's own two extrema annihilate, which is the useless direction;
both are found and the two cases are separated by counting afterwards.

Derivatives of A come from a degree-4 least-squares fit over a 9-point stencil
in u, which is what makes A_uu usable: beta* carries ~1e-12 of integration
noise, and a bare second difference would amplify it by 1/h^2.
"""
import json, math, os, time
import numpy as np
import engine as E
import ahsweep as A


# ------------------------------------------------------------------ the curve
def field(L0, v, lam):
    """Move along v and renormalise back to the seed's scale."""
    n0 = float(np.linalg.norm(L0))
    L = np.asarray(L0, float) + lam * n0 * np.asarray(v, float)
    nn = float(np.linalg.norm(L))
    if nn == 0:
        return None
    return np.ascontiguousarray(L * (n0 / nn))


def bstar_at(L, phi, us, dirhint, bmax=1.5, btol=1e-11, **kw):
    """beta* at a vector of log-radii; NaN where unresolved."""
    s = np.ascontiguousarray(np.exp(np.asarray(us, float)))
    b, st, d0, nf = E.betastar(L, phi, s, dirhint=dirhint, bmax=bmax,
                               btol=btol, **kw)
    b = np.where(st == 0, b, np.nan)
    return b


STENCIL = 9
FITDEG = 4


def derivs(L, phi, u, h, dirhint, **kw):
    """(A, A_u, A_uu, A_uuu) at u from a degree-4 fit on a 9-point stencil."""
    off = np.linspace(-3 * h, 3 * h, STENCIL)
    b = bstar_at(L, phi, u + off, dirhint, **kw)
    m = np.isfinite(b)
    if m.sum() < FITDEG + 2:
        return None
    c = np.polyfit(off[m], b[m], FITDEG)          # highest power first
    p = np.poly1d(c)
    return (float(p(0.0)), float(p.deriv(1)(0.0)),
            float(p.deriv(2)(0.0)), float(p.deriv(3)(0.0)))


# ------------------------------------------------------------------ the solve
def solve_cusp(L0, phi, v, u0, dirhint, lam0=0.0, h=0.03, dlam=2e-3,
               lam_max=0.8, iters=40, **kw):
    """Damped Newton on (A_u, A_uu) = 0 in (u, lambda).  Returns a dict."""
    u, lam = float(u0), float(lam0)
    hist = []
    for it in range(iters):
        L = field(L0, v, lam)
        if L is None:
            return dict(ok=False, why="degenerate")
        okf, tr, dt = A.is_focus(L)
        if not okf:
            return dict(ok=False, why="not_a_focus", lam=lam, u=u)
        d = derivs(L, phi, u, h, dirhint, **kw)
        if d is None:
            return dict(ok=False, why="unresolved", lam=lam, u=u)
        _a, au, auu, auuu = d
        hist.append((u, lam, au, auu))
        # scale-free residual: beta* moves by `rng` over the curve
        if abs(au) < 1e-13 and abs(auu) < 1e-11:
            return dict(ok=True, u=u, lam=lam, au=au, auu=auu, auuu=auuu,
                        iters=it, hist=hist)
        # Jacobian: d/du analytic from the same fit, d/dlambda by difference
        L2 = field(L0, v, lam + dlam)
        d2 = derivs(L2, phi, u, h, dirhint, **kw) if L2 is not None else None
        if d2 is None:
            return dict(ok=False, why="unresolved_lam", lam=lam, u=u)
        dau_dl = (d2[1] - au) / dlam
        dauu_dl = (d2[2] - auu) / dlam
        J = np.array([[auu, dau_dl], [auuu, dauu_dl]])
        if abs(np.linalg.det(J)) < 1e-300:
            return dict(ok=False, why="singular", lam=lam, u=u)
        step = np.linalg.solve(J, -np.array([au, auu]))
        # trust region in both coordinates
        step[0] = float(np.clip(step[0], -4 * h, 4 * h))
        step[1] = float(np.clip(step[1], -0.05, 0.05))
        u += step[0]
        lam += step[1]
        if abs(lam) > lam_max:
            return dict(ok=False, why="lam_out_of_range", lam=lam, u=u)
    return dict(ok=False, why="no_convergence", lam=lam, u=u, hist=hist)


# --------------------------------------------------------------- exploitation
def cross_cusp(L0, phi, v, lam_star, n=260, deltas=(0.0, 5e-4, 2e-3, 8e-3, 3e-2)):
    """Count interior extrema of beta* on both sides of a converged cusp."""
    out = []
    for d in deltas:
        for sgn in ((0,) if d == 0 else (+1, -1)):
            lam = lam_star + sgn * d
            L = field(L0, v, lam)
            if L is None:
                continue
            r = A.evaluate(L, phi, n=n)
            f = r[0] if isinstance(r, tuple) else r
            out.append(dict(lam=float(lam), delta=float(sgn * d),
                            status=f.get("status"),
                            n_extrema=f.get("n_extrema"),
                            n_extrema_robust=f.get("n_extrema_robust"),
                            n_run=f.get("n_run"),
                            fold_margin=f.get("fold_margin"),
                            wiggle_range=f.get("wiggle_range"),
                            extrema=f.get("extrema"),
                            local10=A.coeff_strings(L)))
    return out


def start_points(feat, s_lo, s_hi):
    """Where to launch Newton: the near-shoulder of each outside region, the
    midpoints of the two outside regions, and the interval between the two
    known extrema (which finds the annihilating cusp -- kept as a control)."""
    us = []
    ex = [e["s"] for e in (feat.get("extrema") or [])]
    lo, hi = math.log(s_lo), math.log(s_hi)
    if feat.get("fold_margin_s"):
        us.append(math.log(feat["fold_margin_s"]))
    if ex:
        e0, e1 = math.log(min(ex)), math.log(max(ex))
        us += [0.5 * (lo + e0), 0.5 * (e1 + hi),
               e0 - 0.25 * (e0 - lo), e1 + 0.25 * (hi - e1),
               0.5 * (e0 + e1)]
    else:
        us += [lo + 0.25 * (hi - lo), 0.5 * (lo + hi), lo + 0.75 * (hi - lo)]
    return [u for u in us if lo + 0.02 * (hi - lo) < u < hi - 0.02 * (hi - lo)]
