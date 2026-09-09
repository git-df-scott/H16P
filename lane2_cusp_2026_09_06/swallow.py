"""TASK 4 machinery — the SWALLOW-TAIL (multiplicity-four limit cycle) Newton.

The cusp manifold {D = D_x = D_xx = 0} is 3-dimensional in the 5-dimensional
Cherkas moduli space (a, a20, a11, a01, a10) x section.  Adding D_xxx = 0 is one
more equation, so the swallow-tail set C4 is 2-dimensional -- codimension 1
inside the cusp manifold.

Rather than only WATCHING for a sign change of D_xxx along a fixed-(a,a20) cusp
curve, we can solve for it directly: fix (a, x0) and take

    unknowns  (a11, a01, a10, a20)          4
    equations (D, D_x, D_xx, D_xxx) = 0     4

which is square.  a20 is the second shape parameter, so this lets the Newton
leave the fixed-(a,a20) curve in exactly the one direction that can kill D_xxx.

Perko 1995 Thm 4.3 then requires, with mu1 = a11 (the rotating parameter):
    d_{mu1} != 0, d_{r mu1} != 0, d_{rr mu1} != 0, and for j = 2..n
    d(d,d_r)/d(mu1,mu_j) != 0, d(d,d_rr)/d(mu1,mu_j) != 0, d(d_r,d_rr)/d(mu1,mu_j) != 0.
"""
import json, time
import mpmath as mp
from engine import Engine, V7_of, L_of, V1_of
from cusp import Cusp, solve4, wres

mp.mp.dps = 50

FD = mp.mpf("1e-13")


def wres4(F, x0):
    """Amplitude-weighted residual of (D, D_x, D_xx, D_xxx)."""
    r = abs(mp.mpf(x0) - 1)
    if r == 0:
        r = mp.mpf("1e-300")
    return max(abs(F[0]), abs(F[1]) * r, abs(F[2]) * r * r / 2, abs(F[3]) * r ** 3 / 6)


class Swallow:
    """Unknowns u = (a11, a01, a10, a20) at fixed (a, x0, side)."""

    def __init__(self, eng, a, x0, side=1):
        self.eng = eng
        self.a = mp.mpf(a)
        self.x0 = mp.mpf(x0)
        self.side = side

    def val(self, u, x0=None):
        return self.eng.D(self.a, u[3], u[0], u[1], u[2],
                          self.x0 if x0 is None else x0, side=self.side)

    def F(self, u, x0=None):
        r = self.val(u, x0)
        if r["status"] != "OK":
            return None, r
        return [r["D"], r["Dx"], r["Dxx"], r["Dxxx"]], r

    def jac(self, u, x0=None):
        J = [[mp.mpf(0)] * 4 for _ in range(4)]
        for j in range(4):
            h = FD * max(mp.mpf(1), abs(u[j]))
            up = list(u); up[j] += h
            um = list(u); um[j] -= h
            rp = self.val(up, x0); rm = self.val(um, x0)
            if rp["status"] != "OK" or rm["status"] != "OK":
                return None
            for i, k in enumerate(("D", "Dx", "Dxx", "Dxxx")):
                J[i][j] = (rp[k] - rm[k]) / (2 * h)
        return J

    def newton(self, u, itmax=40, verbose=False, want="1e-24", a20max="5e3",
               parmax="1e5"):
        """Damped Newton with backtracking on wres4 and a divergence guard.

        The guard matters: at SMALL amplitude D_xxx ~ 48 d7 r0^4, and d7 vanishes
        only on the centre variety (V1=V3=V5=V7=0), where D is identically zero.
        So a Newton seeded at small amplitude runs a20 off to infinity instead of
        finding a nondegenerate swallow-tail.  That is a real feature of the
        problem, not a solver failure -- but it must not burn the budget.
        """
        u = list(u)
        best, bestres = list(u), mp.inf
        stall = 0
        for it in range(itmax):
            F, r = self.F(u)
            if F is None:
                return None, r["status"], best
            res = wres4(F, self.x0)
            if res < bestres:
                bestres, best = res, list(u)
                stall = 0
            else:
                stall += 1
            if res < mp.mpf(want):
                break
            J = self.jac(u)
            if J is None:
                return None, "jac-fail", best
            s = solve4(J, F)
            if s is None:
                return None, "singular", best
            sc = max(abs(s[k]) / max(mp.mpf(1), abs(u[k])) for k in range(4))
            lam = mp.mpf(1)
            if sc > mp.mpf("0.3"):
                lam = mp.mpf("0.3") / sc
            # backtracking line search on wres4
            acc = None
            for _ in range(12):
                cand = [u[k] - lam * s[k] for k in range(4)]
                Fc, rc = self.F(cand)
                if Fc is not None:
                    rr = wres4(Fc, self.x0)
                    if rr < res:
                        acc = (cand, rr)
                        break
                lam /= 3
            if verbose:
                print("   it%-2d res=%.3e step=%.3e lam=%.3g a20=%s"
                      % (it, float(res), float(sc), float(lam), mp.nstr(u[3], 12)))
            if acc is None:
                break                       # cannot decrease: at the noise floor
            u = acc[0]
            if abs(u[3]) > mp.mpf(a20max) or max(abs(v) for v in u[:3]) > mp.mpf(parmax):
                return None, "DIVERGED(a20=%s)" % mp.nstr(u[3], 8), best
            if stall > 4:
                break
        F, r = self.F(best)
        if F is None:
            return None, r["status"], best
        res = wres4(F, self.x0)
        if res > mp.mpf(want):
            return None, "no-converge(res=%.2e)" % float(res), best
        return best, r, res


def perko43(sw, u):
    """Perko 1995 Thm 4.3 nondegeneracy quantities, mu1 = a11 (index 0),
    mu_j = a01(1), a10(2), a20(3)."""
    J = sw.jac(u)             # J[i][j] = d^(i)_{mu_j}, i = (d, d_r, d_rr, d_rrr)
    if J is None:
        return None
    out = {"d_mu1": J[0][0], "d_r_mu1": J[1][0], "d_rr_mu1": J[2][0]}
    for j, nm in ((1, "a01"), (2, "a10"), (3, "a20")):
        out["J_d_dr_%s" % nm] = J[0][0] * J[1][j] - J[0][j] * J[1][0]
        out["J_d_drr_%s" % nm] = J[0][0] * J[2][j] - J[0][j] * J[2][0]
        out["J_dr_drr_%s" % nm] = J[1][0] * J[2][j] - J[1][j] * J[2][0]
    out["d_rrr_x"] = J[3][0]          # not a Thm 4.3 condition; logged for context
    vals = [v for k, v in out.items() if k != "d_rrr_x"]
    out["min_abs"] = min(abs(v) for v in vals)
    return out


def try_from_cusp_point(eng, a, a20, mu, x0, side=1, verbose=False):
    """Seed the swallow-tail Newton from a point already ON the cusp manifold."""
    sw = Swallow(eng, a, x0, side=side)
    u0 = [mp.mpf(mu[0]), mp.mpf(mu[1]), mp.mpf(mu[2]), mp.mpf(a20)]
    F0, r0 = sw.F(u0)
    u, r, res = sw.newton(u0, verbose=verbose)
    out = {"a": mp.nstr(a, 30), "x0": mp.nstr(x0, 30), "side": side,
           "seed_a20": mp.nstr(a20, 20),
           "seed_Dxxx": mp.nstr(F0[3], 12) if F0 else None}
    if u is None:
        out["status"] = "FAIL:%s" % r
        if res is not None:
            out["best"] = [mp.nstr(v, 24) for v in res]
        return out
    out.update({"status": "OK",
                "a11": mp.nstr(u[0], 34), "a01": mp.nstr(u[1], 34),
                "a10": mp.nstr(u[2], 34), "a20": mp.nstr(u[3], 34),
                "res": mp.nstr(res, 8),
                "D": mp.nstr(r["D"], 10), "Dx": mp.nstr(r["Dx"], 10),
                "Dxx": mp.nstr(r["Dxx"], 10), "Dxxx": mp.nstr(r["Dxxx"], 10),
                "T": mp.nstr(r["T"], 16),
                "V1": mp.nstr(V1_of(a, u[0], u[1]), 10),
                "L": mp.nstr(L_of(a, u[3], u[0], u[1], u[2]), 10),
                "a20_moved": mp.nstr(u[3] - mp.mpf(a20), 10)})
    pk = perko43(sw, u)
    if pk:
        out["perko43"] = {k: mp.nstr(v, 10) for k, v in pk.items()}
    return out
