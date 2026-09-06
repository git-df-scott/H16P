"""Cycle counting under PROTOCOL rule 1, cusp/swallowtail unfolding, and the
Andronov-Hopf function AH(x) = a11 (Cherkas's rotating parameter).

PROTOCOL rule 1: a cycle is claimed ONLY from a sign change of D on a bracket
[s1,s2] with min(|D(s1)|,|D(s2)|) above a two-tolerance noise estimate
(recompute at looser rtol; noise = 10*|difference| + 5e-12*s).
"""
import json, os, time
import mpmath as mp
from engine import Engine, third_order, V7_of, L_of, V1_of
from cusp import Cusp, solve3, solve4, perko_data, wres

mp.mp.dps = 50


class Noise:
    """Second, looser-tolerance instance of the same integrator, used only for the
    two-tolerance noise estimate of PROTOCOL rule 1."""

    def __init__(self):
        self.loose = Engine(quad=True, env={"CUSP_TOL": "1e-22", "CUSP_ORDER": "18"})

    def est(self, tight, a, a20, mu, x0, side=1):
        r2 = self.loose.D(a, a20, mu[0], mu[1], mu[2], x0, side=side)
        if r2["status"] != "OK":
            return None
        return 10 * abs(tight - r2["D"]) + mp.mpf("5e-12") * abs(x0)

    def close(self):
        self.loose.close()


def sign_changes(c, mu, lo, hi, n, noise=None):
    """Sample D on [lo,hi]; return (samples, certified brackets, uncertified, fails)."""
    xs, ds = [], []
    for i in range(n + 1):
        x = mp.mpf(lo) + (mp.mpf(hi) - mp.mpf(lo)) * i / n
        r = c.val(mu, x)
        xs.append(x)
        ds.append(r["D"] if r["status"] == "OK" else None)
    good, weak, fails = [], [], sum(1 for d in ds if d is None)
    for i in range(n):
        if ds[i] is None or ds[i + 1] is None:
            continue
        if (ds[i] > 0) != (ds[i + 1] > 0):
            m = min(abs(ds[i]), abs(ds[i + 1]))
            if noise is not None:
                nz = max(noise.est(ds[i], c.a, c.a20, mu, xs[i], c.side) or mp.mpf(0),
                         noise.est(ds[i + 1], c.a, c.a20, mu, xs[i + 1], c.side) or mp.mpf(0))
            else:
                nz = mp.mpf(0)
            (good if m > nz else weak).append(
                {"lo": xs[i], "hi": xs[i + 1], "Dlo": ds[i], "Dhi": ds[i + 1],
                 "min_abs": m, "noise": nz})
    return xs, ds, good, weak, fails


def refine_root(c, mu, lo, hi):
    flo = c.val(mu, lo)["D"]
    for _ in range(120):
        m = (lo + hi) / 2
        r = c.val(mu, m)
        if r["status"] != "OK":
            return None
        if (r["D"] > 0) == (flo > 0):
            lo, flo = m, r["D"]
        else:
            hi = m
        if hi - lo < abs(hi) * mp.mpf("1e-30"):
            break
    return (lo + hi) / 2


# ---------------------------------------------------------------- unfolding

def unfold_cusp(c, mu, x0, target):
    """Return mu' such that (D, D_x, D_xx)(x0) = target, to first order.
    target = (q, p, s) as mpf.  Uses the 3x3 parameter Jacobian."""
    M = c.jac_mu(mu, x0)
    if M is None:
        return None
    F, r = c.F(mu, x0)
    if F is None:
        return None
    rhs = [F[k] - mp.mpf(target[k]) for k in range(3)]
    d = solve3(M, rhs)
    if d is None:
        return None
    return [mu[k] - d[k] for k in range(3)]


def triple_confirm(c, mu, x0, noise, frac="0.10", nsample=400):
    """TASK 2 confirmation: perturb INTO the cuspidal region and check that the
    triple root resolves into THREE certified sign changes of D in one nest.

    Near the cusp, D(x0+u) ~ c1 u + (1/6) D_xxx u^3.  Three simple roots
    (u = 0, +-sqrt(-6 c1 / D_xxx)) require c1 * D_xxx < 0.  We choose the
    separation to be `frac` of the amplitude r0 = x0 - 1.
    """
    F, r = c.F(mu, x0)
    if F is None:
        return {"status": "eval-fail"}
    Dxxx = r["Dxxx"]
    u = mp.mpf(frac) * (mp.mpf(x0) - 1)
    c1 = -Dxxx * u * u / 6          # so that c1*Dxxx < 0
    mu2 = unfold_cusp(c, mu, x0, (mp.mpf(0), c1, mp.mpf(0)))
    if mu2 is None:
        return {"status": "unfold-fail"}
    lo, hi = x0 - 3 * u, x0 + 3 * u
    if lo <= 1:
        lo = 1 + (mp.mpf(x0) - 1) / 100
    xs, ds, good, weak, fails = sign_changes(c, mu2, lo, hi, nsample, noise)
    roots = [refine_root(c, mu2, b["lo"], b["hi"]) for b in good]
    return {"status": "OK",
            "mu_perturbed": [mp.nstr(v, 34) for v in mu2],
            "delta_mu": [mp.nstr(mu2[k] - mu[k], 10) for k in range(3)],
            "window": [mp.nstr(lo, 20), mp.nstr(hi, 20)],
            "n_certified_sign_changes": len(good),
            "n_uncertified": len(weak), "fails": fails,
            "Dxxx_at_cusp": mp.nstr(Dxxx, 12),
            "predicted_separation": mp.nstr(u, 10),
            "roots": [mp.nstr(x, 24) for x in roots if x is not None],
            "brackets": [{"lo": mp.nstr(b["lo"], 20), "hi": mp.nstr(b["hi"], 20),
                          "Dlo": mp.nstr(b["Dlo"], 10), "Dhi": mp.nstr(b["Dhi"], 10),
                          "min_abs": mp.nstr(b["min_abs"], 8),
                          "noise": mp.nstr(b["noise"], 8)} for b in good]}


# ------------------------------------------------- Andronov-Hopf function ---

def AH(c, mu, x, lo=None, hi=None, tol="1e-28"):
    """Solve D(x; a11) = 0 for a11 with (a01, a10) fixed.  a11 is Cherkas's
    rotating parameter, so D is strictly monotone in it (Duff/Perko) and
    bisection is safe."""
    a11 = mu[0]
    if lo is None:
        lo, hi = a11 - 4, a11 + 4
    lo, hi = mp.mpf(lo), mp.mpf(hi)
    def f(v):
        r = c.val([v, mu[1], mu[2]], x)
        return r["D"] if r["status"] == "OK" else None
    flo, fhi = f(lo), f(hi)
    if flo is None or fhi is None or (flo > 0) == (fhi > 0):
        return None
    for _ in range(200):
        m = (lo + hi) / 2
        fm = f(m)
        if fm is None:
            return None
        if (fm > 0) == (flo > 0):
            lo, flo = m, fm
        else:
            hi = m
        if hi - lo < mp.mpf(tol) * max(mp.mpf(1), abs(hi)):
            break
    return (lo + hi) / 2


def ah_sweep(c, mu, lo, hi, n=60):
    """Sample AH(x) and count interior extrema.  TWO extrema <=> three cycles in
    the nest; THREE extrema <=> four cycles in the nest (Cherkas's criterion)."""
    xs, vs = [], []
    for i in range(n + 1):
        x = mp.mpf(lo) + (mp.mpf(hi) - mp.mpf(lo)) * i / n
        v = AH(c, mu, x)
        xs.append(x); vs.append(v)
    ext = []
    for i in range(1, len(vs) - 1):
        if vs[i - 1] is None or vs[i] is None or vs[i + 1] is None:
            continue
        if (vs[i] - vs[i - 1] > 0) != (vs[i + 1] - vs[i] > 0):
            ext.append({"x": mp.nstr(xs[i], 12), "a11": mp.nstr(vs[i], 16),
                        "type": "max" if vs[i] > vs[i - 1] else "min"})
    return {"x": [mp.nstr(v, 12) for v in xs],
            "AH": [None if v is None else mp.nstr(v, 20) for v in vs],
            "n_defined": sum(1 for v in vs if v is not None),
            "n_interior_extrema": len(ext), "extrema": ext}
