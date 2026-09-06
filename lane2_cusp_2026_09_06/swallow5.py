"""Minimum-norm swallow-tail Newton over ALL FIVE Cherkas parameters.

The swallow-tail set {D = D_x = D_xx = D_xxx = 0} is 4 equations; the moduli space
(a, a20, a11, a01, a10) is 5-dimensional, so at a fixed section point x0 the set is
generically a CURVE, not a point.  Forcing the Newton to use one designated fourth
parameter (as swallow.py does with a20) can therefore send it off to infinity even
when a swallow-tail exists in some other direction -- exactly what happens on the
a = 3 curves, where the centre curve a20_c(3) = +0.49 lies outside the admissible
region a20 < a-3 = 0, so no a20 whatever can make d7 vanish.

Here the step is the MINIMUM-NORM solution of J s = F (J is 4x5), measured in
RELATIVE parameter size, so the Newton moves whichever parameters are cheapest and
never has to blow one of them up.

Guards: the antisaddle condition L = 2a - a01 - a10 - 2a20 > 0 (otherwise A is not
a focus and there is no nest), a != 2, a != 1/3, and a bound on |parameters|.
"""
import mpmath as mp
from engine import Engine, L_of, V1_of

mp.mp.dps = 50

FD = mp.mpf("1e-13")
NAMES = ("a11", "a01", "a10", "a20", "a")


def wres4(F, x0):
    r = abs(mp.mpf(x0) - 1)
    if r == 0:
        r = mp.mpf("1e-300")
    return max(abs(F[0]), abs(F[1]) * r, abs(F[2]) * r * r / 2, abs(F[3]) * r ** 3 / 6)


def _solve(A, b):
    """Dense solve with partial pivoting (n up to 5)."""
    n = len(b)
    M = [list(A[i]) + [b[i]] for i in range(n)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[p][c]) == 0:
            return None
        M[c], M[p] = M[p], M[c]
        for r in range(n):
            if r == c:
                continue
            f = M[r][c] / M[c][c]
            for k in range(c, n + 1):
                M[r][k] -= f * M[c][k]
    return [M[i][n] / M[i][i] for i in range(n)]


class Swallow5:
    """v = (a11, a01, a10, a20, a) at fixed x0 and section side."""

    def __init__(self, eng, x0, side=1):
        self.eng = eng
        self.x0 = mp.mpf(x0)
        self.side = side

    def val(self, v):
        return self.eng.D(v[4], v[3], v[0], v[1], v[2], self.x0, side=self.side)

    def F(self, v):
        r = self.val(v)
        if r["status"] != "OK":
            return None, r
        return [r["D"], r["Dx"], r["Dxx"], r["Dxxx"]], r

    def jac(self, v):
        """4x5 Jacobian: J[i][j] = d^(i) / d v_j."""
        J = [[mp.mpf(0)] * 5 for _ in range(4)]
        for j in range(5):
            h = FD * max(mp.mpf(1), abs(v[j]))
            vp = list(v); vp[j] += h
            vm = list(v); vm[j] -= h
            rp = self.val(vp); rm = self.val(vm)
            if rp["status"] != "OK" or rm["status"] != "OK":
                return None
            for i, k in enumerate(("D", "Dx", "Dxx", "Dxxx")):
                J[i][j] = (rp[k] - rm[k]) / (2 * h)
        return J

    def ok(self, v, parmax):
        if max(abs(t) for t in v) > mp.mpf(parmax):
            return False
        if abs(v[4] - 2) < mp.mpf("1e-6") or abs(v[4] - mp.mpf(1) / 3) < mp.mpf("1e-6"):
            return False
        return L_of(v[4], v[3], v[0], v[1], v[2]) > 0

    def newton(self, v, itmax=50, verbose=False, want="1e-24", parmax="2e3"):
        v = list(v)
        best, bestres = list(v), mp.inf
        stall = 0
        for it in range(itmax):
            F, r = self.F(v)
            if F is None:
                return None, r["status"], best
            res = wres4(F, self.x0)
            if res < bestres:
                bestres, best, stall = res, list(v), 0
            else:
                stall += 1
            if res < mp.mpf(want):
                break
            J = self.jac(v)
            if J is None:
                return None, "jac-fail", best
            # minimum-norm step in RELATIVE parameter size:
            #   s = W J'^T (J' J'^T)^{-1} F,   J' = J W,  W = diag(max(1,|v_j|))
            W = [max(mp.mpf(1), abs(t)) for t in v]
            Jp = [[J[i][j] * W[j] for j in range(5)] for i in range(4)]
            G = [[sum(Jp[i][k] * Jp[j][k] for k in range(5)) for j in range(4)]
                 for i in range(4)]
            y = _solve(G, F)
            if y is None:
                return None, "singular", best
            s = [W[j] * sum(Jp[i][j] * y[i] for i in range(4)) for j in range(5)]
            sc = max(abs(s[j]) / W[j] for j in range(5))
            lam = mp.mpf(1)
            if sc > mp.mpf("0.3"):
                lam = mp.mpf("0.3") / sc
            acc = None
            for _ in range(14):
                cand = [v[j] - lam * s[j] for j in range(5)]
                if self.ok(cand, parmax):
                    Fc, rc = self.F(cand)
                    if Fc is not None and wres4(Fc, self.x0) < res:
                        acc = cand
                        break
                lam /= 3
            if verbose:
                print("   it%-2d res=%.3e step=%.2e lam=%.2g  a=%s a20=%s"
                      % (it, float(res), float(sc), float(lam),
                         mp.nstr(v[4], 10), mp.nstr(v[3], 10)), flush=True)
            if acc is None:
                break
            v = acc
            if stall > 5:
                break
        F, r = self.F(best)
        if F is None:
            return None, r["status"], best
        res = wres4(F, self.x0)
        if res > mp.mpf(want):
            return None, "no-converge(res=%.2e)" % float(res), best
        return best, r, res

    def perko43(self, v):
        """Perko 1995 Thm 4.3 quantities with mu1 = a11 (the rotating parameter)."""
        J = self.jac(v)
        if J is None:
            return None
        out = {"d_mu1": J[0][0], "d_r_mu1": J[1][0], "d_rr_mu1": J[2][0]}
        for j in range(1, 5):
            nm = NAMES[j]
            out["J_d_dr_%s" % nm] = J[0][0] * J[1][j] - J[0][j] * J[1][0]
            out["J_d_drr_%s" % nm] = J[0][0] * J[2][j] - J[0][j] * J[2][0]
            out["J_dr_drr_%s" % nm] = J[1][0] * J[2][j] - J[1][j] * J[2][0]
        out["min_abs"] = min(abs(t) for t in out.values())
        return out
