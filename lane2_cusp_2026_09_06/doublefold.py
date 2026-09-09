#!/usr/bin/env python3
"""SEPARATED DOUBLE-FOLD system.

   F(mu, s1, s2) = ( D(s1), D_s(s1), D(s2), D_s(s2) ) = 0

with s1 != s2 in the same regular return-domain component, D_ss(s1), D_ss(s2)
both nonzero.  Square slice: unknowns (p1, p2, s1, s2) for a chosen parameter
pair; the rest of mu is frozen for that slice.

Jacobian: the s-columns are free from the jet ((D_s,D_ss) at each section);
the parameter columns are central differences.

Separation guard: |s1-s2| must stay above SEP_MIN -- s1 -> s2 is a cusp /
higher-multiplicity boundary event, not two separated folds.
"""
import mpmath as mp
from cusp import solve4, FD_H
import solver as SV

SEP_MIN = mp.mpf("0.05")
# The system has a TRIVIAL branch at the focus: D(1)=0 identically and D_s(1)=V1,
# so s->1 with V1->0 solves it.  That is the Bautin/Hopf boundary (PROTOCOL
# mechanism (a), capped at 3 small cycles), not a separated double fold.
MIN_AMP = mp.mpf("0.08")
PARAMS = ("a", "a20", "a11", "a01", "a10")

def get(c, nm):
    return getattr(c, nm) if nm in ("a", "a20") else None

class DF:
    def __init__(self, c, mu):
        self.c = c
        self.mu = dict(mu)          # a11, a01, a10 (a, a20 live on c)
    def _u(self):
        return [self.mu["a11"], self.mu["a01"], self.mu["a10"]]
    def val(self, s):
        return self.c.val(self._u(), s)
    def F(self, s1, s2):
        r1, r2 = self.val(s1), self.val(s2)
        if r1["status"] != "OK" or r2["status"] != "OK":
            return None, (r1, r2)
        return [r1["D"], r1["Dx"], r2["D"], r2["Dx"]], (r1, r2)
    def setp(self, nm, v):
        if nm in ("a", "a20"): setattr(self.c, nm, v)
        else: self.mu[nm] = v
    def getp(self, nm):
        return getattr(self.c, nm) if nm in ("a", "a20") else self.mu[nm]

def jac(df, s1, s2, pair, r1, r2):
    """4x4 d(D(s1),Dx(s1),D(s2),Dx(s2)) / d(p1,p2,s1,s2)."""
    M = [[mp.mpf(0)]*4 for _ in range(4)]
    for j, nm in enumerate(pair):
        old = df.getp(nm)
        h = mp.mpf("1e-9")*max(mp.mpf(1), abs(old))
        df.setp(nm, old + h); a1, a2 = df.val(s1), df.val(s2)
        df.setp(nm, old - h); b1, b2 = df.val(s1), df.val(s2)
        df.setp(nm, old)
        if any(z["status"] != "OK" for z in (a1, a2, b1, b2)): return None
        M[0][j] = (a1["D"]  - b1["D"]) /(2*h); M[1][j] = (a1["Dx"] - b1["Dx"])/(2*h)
        M[2][j] = (a2["D"]  - b2["D"]) /(2*h); M[3][j] = (a2["Dx"] - b2["Dx"])/(2*h)
    M[0][2], M[1][2] = r1["Dx"], r1["Dxx"]
    M[2][3], M[3][3] = r2["Dx"], r2["Dxx"]
    return M

def solve(df, s1, s2, pair, tol="1e-22", itmax=50, verbose=True):
    s1, s2 = mp.mpf(s1), mp.mpf(s2); prev = mp.inf
    best = (mp.inf, None)
    for it in range(itmax):
        if abs(s1 - s2) < SEP_MIN:
            return SV.Result("SEPARATION_LOST", detail="|s1-s2|=%s" % mp.nstr(abs(s1-s2),4), its=it)
        if min(abs(s1-1), abs(s2-1)) < MIN_AMP:
            return SV.Result("AMPLITUDE_LOST",
                             detail="fold collapsing onto the focus: s-1 = (%s, %s)"
                                    % (mp.nstr(s1-1,4), mp.nstr(s2-1,4)), its=it)
        F, (r1, r2) = df.F(s1, s2)
        if F is None:
            st = r1["status"] if r1["status"] != "OK" else r2["status"]
            return SV.Result(SV.GUARD_EXIT if "GUARD" in st else SV.RETURN_FAILED,
                             detail=st, its=it)
        sc = [abs(F[0]), abs(F[1])*abs(s1-1), abs(F[2]), abs(F[3])*abs(s2-1)]
        m = max(sc)
        if m < best[0]: best = (m, (s1, s2, dict(df.mu), df.c.a, df.c.a20))
        if verbose and it < 12:
            print("   it%-2d max|F|scaled=%.3e  D(s1)=%+.3e D(s2)=%+.3e  s=(%.5f,%.5f)"
                  % (it, float(m), float(F[0]), float(F[2]), float(s1), float(s2)), flush=True)
        if m < mp.mpf(tol):
            return SV.Result(SV.CONVERGED, root=(s1, s2, dict(df.mu)), resid=F, scaled=sc,
                             its=it, best=best[1])
        J = jac(df, s1, s2, pair, r1, r2)
        if J is None: return SV.Result(SV.RETURN_FAILED, detail="jacobian", its=it, best=best[1])
        st = solve4(J, F)
        if st is None: return SV.Result(SV.SINGULAR, its=it, best=best[1])
        scale = max(abs(st[k]) for k in range(4))
        if scale >= prev/2 and it >= 12:
            return SV.Result(SV.STALLED, detail="step stopped decreasing",
                             resid=F, scaled=sc, its=it, best=best[1])
        prev = scale
        # DAMPING: shrink the step until the amplitude/separation guards hold,
        # so the solver cannot jump onto the trivial focus branch in one step.
        lam = mp.mpf(1)
        for _ in range(40):
            t1, t2 = s1 - lam*st[2], s2 - lam*st[3]
            if (min(abs(t1-1), abs(t2-1)) >= MIN_AMP and abs(t1-t2) >= SEP_MIN):
                break
            lam /= 2
        else:
            return SV.Result("AMPLITUDE_LOST",
                             detail="no damped step keeps the folds off the focus", its=it,
                             resid=F, scaled=sc, best=best[1])
        for j, nm in enumerate(pair): df.setp(nm, df.getp(nm) - lam*st[j])
        s1 -= lam*st[2]; s2 -= lam*st[3]
    return SV.Result(SV.ITERATION_LIMIT, resid=F, scaled=sc, its=itmax, best=best[1])
