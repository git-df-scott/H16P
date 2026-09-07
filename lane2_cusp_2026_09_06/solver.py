#!/usr/bin/env python3
"""Solver success contract.

No consumer may read a non-null parameter vector as convergence.  Every solve
returns a Result carrying an explicit status, the FINAL defining residual
vector (componentwise), its scaling, and the best diagnostic iterate SEPARATELY
from the converged root (which is None unless status is CONVERGED).

Residual scaling: component i of Phi=(D,D_x,D_xx,D_xxx) is compared against the
Taylor term it controls at radius r0 = x0-1, i.e. |Phi_i| * r0^i / i!.  A small
weighted norm may not conceal a large unresolved derivative: EVERY component
must pass, not just the norm.
"""
import mpmath as mp
from cusp import solve4, FD_H

CONVERGED, STALLED, SINGULAR = "CONVERGED", "STALLED", "SINGULAR"
RETURN_FAILED, GUARD_EXIT, ITERATION_LIMIT = "RETURN_FAILED", "GUARD_EXIT", "ITERATION_LIMIT"

class Result:
    def __init__(self, status, root=None, best=None, resid=None, scaled=None,
                 detail="", its=0, jac_det=None):
        self.status, self.root, self.best = status, root, best
        self.resid, self.scaled, self.detail, self.its = resid, scaled, detail, its
        self.jac_det = jac_det
    @property
    def converged(self):
        return self.status == CONVERGED
    def __repr__(self):
        r = "Result(%s" % self.status
        if self.scaled is not None:
            r += ", scaled=[%s]" % ", ".join(mp.nstr(v, 3) for v in self.scaled)
        return r + (", %s)" % self.detail if self.detail else ")")

def scaled_residual(F, x0):
    """Componentwise: |Phi_i| * r0^i / i!  -- the Taylor term each one controls."""
    r0 = abs(mp.mpf(x0) - 1)
    fact = [1, 1, 2, 6]
    return [abs(F[i])*r0**i/fact[i] for i in range(len(F))]

def newton_swallow(c, u0, x00, tol="1e-24", itmax=60, verbose=False):
    """Solve (D,D_x,D_xx,D_xxx)=0 in (a11,a01,a10,x0) at fixed shape.
    Needs only the degree-4 jet: det Jac(F,G)/(u,s) = det(F_u) * D_ssss."""
    u = [mp.mpf(v) for v in u0]; x0 = mp.mpf(x00)
    prev = mp.inf; best = (mp.inf, list(u), x0, None); its = 0
    for it in range(itmax):
        its = it + 1
        r = c.val(list(u), x0)
        if r["status"] != "OK":
            return Result(GUARD_EXIT if "GUARD" in r["status"] else RETURN_FAILED,
                          best=(best[1], best[2]), detail=r["status"], its=its)
        F = [r["D"], r["Dx"], r["Dxx"], r["Dxxx"]]
        sc = scaled_residual(F, x0); m = max(sc)
        if m < best[0]: best = (m, list(u), x0, r)
        if m < mp.mpf(tol):
            J = _jac(c, u, x0, r)
            det = None
            if J is not None:
                M = mp.matrix(4, 4)
                for i in range(4):
                    for j in range(4): M[i, j] = J[i][j]
                det = mp.det(M)
            return Result(CONVERGED, root=(list(u), x0), best=(list(u), x0),
                          resid=F, scaled=sc, its=its, jac_det=det)
        J = _jac(c, u, x0, r)
        if J is None:
            return Result(RETURN_FAILED, best=(best[1], best[2]), detail="jacobian", its=its)
        s = solve4(J, F)
        if s is None:
            return Result(SINGULAR, best=(best[1], best[2]), detail="singular Jacobian", its=its)
        step = max(abs(s[k])/max(mp.mpf(1), abs(([*u, x0])[k])) for k in range(4))
        if verbose:
            print("   it%-2d maxscaled=%.3e Dxxx=%+.4e step=%.2e"
                  % (it, float(m), float(F[3]), float(step)), flush=True)
        if step >= prev/2 and it >= 4:
            rb = best[3]
            Fb = [rb["D"], rb["Dx"], rb["Dxx"], rb["Dxxx"]] if rb else None
            return Result(STALLED, best=(best[1], best[2]), resid=Fb,
                          scaled=scaled_residual(Fb, best[2]) if Fb else None,
                          detail="step stopped decreasing", its=its)
        prev = step
        u = [u[k] - s[k] for k in range(3)]; x0 = x0 - s[3]
    rb = best[3]
    Fb = [rb["D"], rb["Dx"], rb["Dxx"], rb["Dxxx"]] if rb else None
    return Result(ITERATION_LIMIT, best=(best[1], best[2]), resid=Fb,
                  scaled=scaled_residual(Fb, best[2]) if Fb else None, its=its)

def _jac(c, u, x0, r):
    M = [[mp.mpf(0)]*4 for _ in range(4)]
    for j in range(3):
        h = FD_H[j]*max(mp.mpf(1), abs(u[j]))
        up, um = list(u), list(u); up[j] += h; um[j] -= h
        rp, rm = c.val(up, x0), c.val(um, x0)
        if rp["status"] != "OK" or rm["status"] != "OK": return None
        for i, kk in enumerate(("D", "Dx", "Dxx", "Dxxx")):
            M[i][j] = (rp[kk] - rm[kk])/(2*h)
    for i, kk in enumerate(("Dx", "Dxx", "Dxxx", "Dxxxx")):
        M[i][3] = r[kk]
    return M
