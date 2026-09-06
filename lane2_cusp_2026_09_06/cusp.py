"""Lane 2: the manifold of TRIPLE limit cycles (cusp manifold) of the Cherkas family,
entered from the Bautin small-amplitude region and continued out in amplitude,
watching D_xxx for a sign change (= swallowtail = multiplicity-four limit cycle).

Coordinates
-----------
shape parameters   (a, a20)        held fixed along one cusp curve
unfolding params   mu = (a11, a01, a10)     mu1 = a11 is Cherkas's ROTATING parameter,
                                            so Perko's d_{mu1} != 0 is the natural
                                            nondegeneracy anchor.
section coordinate x0 on { y = -1, x > 1 }

cusp equations     F(mu, x0) = ( D, D_x, D_xx ) = 0        (3 equations, 4 unknowns)
swallowtail        additionally D_xxx = 0

All x-derivatives come from the degree-3 jet (exact to integration accuracy);
only the PARAMETER derivatives are finite differences, and those are taken in
binary128 where a centred difference with h ~ 1e-11 is good to ~1e-21.
"""
import json, math, os, sys, time
import mpmath as mp
from engine import Engine, third_order, V7_of, L_of

mp.mp.dps = 50

# --------------------------------------------------------------------- utils

def wres(F, x0):
    """Amplitude-weighted residual of (D, D_x, D_xx): the size, in units of
    length, of each Taylor term of D over the scale r = x0 - 1."""
    r = abs(mp.mpf(x0) - 1)
    if r == 0:
        r = mp.mpf("1e-300")
    return max(abs(F[0]), abs(F[1]) * r, 0.5 * abs(F[2]) * r * r)


def solve3(A, b):
    """Gaussian elimination with partial pivoting on a 3x3 (or nxn) system."""
    n = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[p][c]) == 0.0:
            return None
        M[c], M[p] = M[p], M[c]
        for r in range(n):
            if r == c:
                continue
            f = M[r][c] / M[c][c]
            for k in range(c, n + 1):
                M[r][k] -= f * M[c][k]
    return [M[i][n] / M[i][i] for i in range(n)]


def nullvec(A):
    """Unit null vector of a (n-1) x n matrix A (n=4 here), via cofactor expansion."""
    n = len(A[0])
    v = []
    for j in range(n):
        sub = [[A[i][k] for k in range(n) if k != j] for i in range(len(A))]
        v.append(((-1) ** j) * det(sub))
    nrm = mp.sqrt(sum(x * x for x in v))
    if nrm == 0:
        return None
    return [x / nrm for x in v]


def det(M):
    n = len(M)
    if n == 1:
        return M[0][0]
    if n == 2:
        return M[0][0] * M[1][1] - M[0][1] * M[1][0]
    if n == 3:
        return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
                - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
    raise ValueError


# ------------------------------------------------------------ residual + jac

FD_H = [mp.mpf("1e-13")] * 3   # parameter FD steps for (a11, a01, a10)


class Cusp:
    def __init__(self, eng, a, a20, side=1):
        self.eng = eng
        self.a = mp.mpf(a)
        self.a20 = mp.mpf(a20)
        self.side = side

    def val(self, mu, x0):
        return self.eng.D(self.a, self.a20, mu[0], mu[1], mu[2], x0, side=self.side)

    def F(self, mu, x0):
        r = self.val(mu, x0)
        if r["status"] != "OK":
            return None, r
        return [r["D"], r["Dx"], r["Dxx"]], r

    def jac_mu(self, mu, x0):
        """3x3 matrix M[i][j] = d^(i) / d mu_j  for i in (d, d_r, d_rr)."""
        M = [[mp.mpf(0)] * 3 for _ in range(3)]
        for j in range(3):
            h = FD_H[j] * max(mp.mpf(1), abs(mu[j]))
            mup = list(mu); mup[j] += h
            mum = list(mu); mum[j] -= h
            rp = self.val(mup, x0); rm = self.val(mum, x0)
            if rp["status"] != "OK" or rm["status"] != "OK":
                return None
            M[0][j] = (rp["D"] - rm["D"]) / (2 * h)
            M[1][j] = (rp["Dx"] - rm["Dx"]) / (2 * h)
            M[2][j] = (rp["Dxx"] - rm["Dxx"]) / (2 * h)
        return M

    def jac_full(self, mu, x0, r=None):
        """3x4 Jacobian of (D,Dx,Dxx) w.r.t. (a11,a01,a10,x0)."""
        M = self.jac_mu(mu, x0)
        if M is None:
            return None
        if r is None:
            r = self.val(mu, x0)
            if r["status"] != "OK":
                return None
        col = [r["Dx"], r["Dxx"], r["Dxxx"]]        # d/dx0 of (D, D_x, D_xx)
        return [M[i] + [col[i]] for i in range(3)]

    # -------------------------------------------------------------- Newton
    def newton_mu(self, mu, x0, itmax=40, verbose=False, want="1e-26"):
        """Solve F(mu; x0)=0 in mu at fixed x0.  Iterate until the Newton step
        STAGNATES (we are at the engine's own noise floor), then report."""
        mu = list(mu)
        best, bestres, bestr = list(mu), mp.inf, None
        prev = mp.inf
        for it in range(itmax):
            F, r = self.F(mu, x0)
            if F is None:
                return None, r["status"]
            res = wres(F, x0)
            if res < bestres:
                bestres, best, bestr = res, list(mu), r
            M = self.jac_mu(mu, x0)
            if M is None:
                return None, "jac-fail"
            s = solve3(M, F)
            if s is None:
                return None, "singular"
            sc = max(abs(s[k]) / max(mp.mpf(1), abs(mu[k])) for k in range(3))
            if verbose:
                print("   it%-2d res=%.3e |F|=(%.2e,%.2e,%.2e) step=%.2e"
                      % (it, float(res), float(abs(F[0])), float(abs(F[1])), float(abs(F[2])), float(sc)))
            if sc >= prev / 2 and it >= 3:
                break                      # stagnated: at the noise floor
            prev = sc
            mu = [mu[k] - s[k] for k in range(3)]
            if sc < mp.mpf("1e-45"):
                break
        F, r = self.F(best, x0)
        if F is None:
            return None, r["status"]
        if wres(F, x0) > mp.mpf(want):
            return None, "no-converge(res=%.2e)" % wres(F, x0)
        return best, r

    # ------------------------------------------------- pseudo-arclength step
    def newton_arc(self, z, tang, zpred, tol="1e-26", itmax=40):
        """Solve F(z)=0 with tang.(z - zpred) = 0.  z = (a11,a01,a10,x0)."""
        z = list(z)
        for it in range(itmax):
            mu, x0 = z[:3], z[3]
            F, r = self.F(mu, x0)
            if F is None:
                return None, r["status"]
            J = self.jac_full(mu, x0, r)
            if J is None:
                return None, "jac-fail"
            A = J + [tang[:]]
            b = F + [sum(tang[k] * (z[k] - zpred[k]) for k in range(4))]
            s = solve4(A, b)
            if s is None:
                return None, "singular"
            z = [z[k] - s[k] for k in range(4)]
            sc = max(abs(s[k]) / max(mp.mpf(1), abs(z[k])) for k in range(4))
            if sc < mp.mpf(tol):
                F, r = self.F(z[:3], z[3])
                if F is None:
                    return None, r["status"]
                return z, r
        F, r = self.F(z[:3], z[3])
        if F is not None and wres(F, z[3]) < mp.mpf(tol) * 100:
            return z, r
        return None, "no-converge"


def solve4(A, b):
    n = 4
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[p][c]) == 0.0:
            return None
        M[c], M[p] = M[p], M[c]
        for r in range(n):
            if r == c:
                continue
            f = M[r][c] / M[c][c]
            for k in range(c, n + 1):
                M[r][k] -= f * M[c][k]
    return [M[i][n] / M[i][i] for i in range(n)]


# ------------------------------------------------------- Perko nondegeneracy

def perko_data(M):
    """Perko 1995 Thm 4.3 nondegeneracy quantities from M[i][j] = d^(i)_{mu_j}.
    mu_1 = a11 (index 0), mu_j for j = a01(1), a10(2)."""
    out = {"d_mu1": M[0][0], "d_r_mu1": M[1][0], "d_rr_mu1": M[2][0]}
    for j, nm in ((1, "a01"), (2, "a10")):
        out["J_d_dr_%s" % nm]   = M[0][0] * M[1][j] - M[0][j] * M[1][0]
        out["J_d_drr_%s" % nm]  = M[0][0] * M[2][j] - M[0][j] * M[2][0]
        out["J_dr_drr_%s" % nm] = M[1][0] * M[2][j] - M[1][j] * M[2][0]
    out["min_abs"] = min(abs(v) for v in out.values())
    return out
