"""THE (3+1) ROUTE: a cusp point whose focus is simultaneously weak.

Motivation.  At a cusp point the displacement has a TRIPLE root at x0.  On the
right section, just outside the focus D ~ (exp(int div) - 1)(x-1), so sgn D =
sgn V1; just inside the cusp point D ~ (1/6) D_xxx (x-x0)^3, so sgn D =
-sgn D_xxx.  Hence D has an odd number of extra zeros between the focus and the
triple cycle -- i.e. a FOURTH limit cycle in the same nest once the cusp is
unfolded into three -- exactly when

        sgn(V1) * sgn(D_xxx) > 0.

Over all 15289 non-degenerate cusp records computed in this lane that product is
NEGATIVE without a single exception.  It can only flip where D_xxx = 0 (the
swallow-tail, which §II.4 shows is not reachable from the Bautin end) or where

        V1 = a11 + a01 - 2a - 1 = 0,

a Hopf bifurcation at the focus -- and that is an ALGEBRAIC condition costing no
integration at all.  V1 = 0 can simply be solved: a01 = 2a + 1 - a11.  Imposing
it leaves the square system

        unknowns  (a11, a10, a20)          3
        equations (D, D_x, D_xx) = 0       3      at fixed (a, x0)

whose solutions are cusp points with a weak focus.  Crossing V1 = 0 along the
cusp manifold flips the indicator and delivers the fourth cycle.

Bautin is not violated: only ONE of the four cycles bifurcates from the focus;
the other three sit at x0, of normal size.
"""
import json, os, time
import mpmath as mp
from engine import Engine, L_of, V1_of
from cusp import solve3, wres

mp.mp.dps = 50

FD = mp.mpf("1e-13")


class CuspV1:
    """Unknowns v = (a11, a10, a20) at fixed (a, x0); a01 = 2a + 1 - a11 slaves
    V1 to zero identically, so the focus is weak at every iterate."""

    def __init__(self, eng, a, x0, side=1):
        self.eng = eng
        self.a = mp.mpf(a)
        self.x0 = mp.mpf(x0)
        self.side = side

    def a01_of(self, a11):
        return 2 * self.a + 1 - a11

    def val(self, v, x0=None):
        return self.eng.D(self.a, v[2], v[0], self.a01_of(v[0]), v[1],
                          self.x0 if x0 is None else x0, side=self.side)

    def F(self, v):
        r = self.val(v)
        if r["status"] != "OK":
            return None, r
        return [r["D"], r["Dx"], r["Dxx"]], r

    def jac(self, v):
        M = [[mp.mpf(0)] * 3 for _ in range(3)]
        for j in range(3):
            h = FD * max(mp.mpf(1), abs(v[j]))
            vp = list(v); vp[j] += h
            vm = list(v); vm[j] -= h
            rp = self.val(vp); rm = self.val(vm)
            if rp["status"] != "OK" or rm["status"] != "OK":
                return None
            for i, k in enumerate(("D", "Dx", "Dxx")):
                M[i][j] = (rp[k] - rm[k]) / (2 * h)
        return M

    def ok(self, v, parmax=mp.mpf("5e3")):
        a01 = self.a01_of(v[0])
        if max(abs(t) for t in (v[0], v[1], v[2], a01)) > parmax:
            return False
        return L_of(self.a, v[2], v[0], a01, v[1]) > 0

    def amplitude(self, v):
        r0 = (self.x0 - 1) if self.side > 0 else (1 - self.x0)
        best = None
        for fac in (mp.mpf("0.5"), mp.mpf("1.7")):
            xr = 1 + r0 * fac * (1 if self.side > 0 else -1)
            q = self.val(v, xr)
            if q["status"] == "OK":
                b = abs(q["D"])
                best = b if best is None else max(best, b)
        return best

    def newton(self, v, itmax=40, verbose=False, want="1e-24"):
        v = list(v)
        best, bestres = list(v), mp.inf
        stall = 0
        for it in range(itmax):
            F, r = self.F(v)
            if F is None:
                return None, r["status"], best
            res = wres(F, self.x0)
            if res < bestres:
                bestres, best, stall = res, list(v), 0
            else:
                stall += 1
            if res < mp.mpf(want):
                break
            M = self.jac(v)
            if M is None:
                return None, "jac-fail", best
            s = solve3(M, F)
            if s is None:
                return None, "singular", best
            sc = max(abs(s[k]) / max(mp.mpf(1), abs(v[k])) for k in range(3))
            lam = mp.mpf(1)
            if sc > mp.mpf("0.3"):
                lam = mp.mpf("0.3") / sc
            acc = None
            for _ in range(14):
                cand = [v[k] - lam * s[k] for k in range(3)]
                if self.ok(cand):
                    Fc, rc = self.F(cand)
                    if Fc is not None and wres(Fc, self.x0) < res:
                        acc = cand
                        break
                lam /= 3
            if verbose:
                print("   it%-2d res=%.3e step=%.2e lam=%.2g a11=%s a20=%s"
                      % (it, float(res), float(sc), float(lam),
                         mp.nstr(v[0], 10), mp.nstr(v[2], 10)), flush=True)
            if acc is None or stall > 4:
                break
            v = acc
        F, r = self.F(best)
        if F is None:
            return None, r["status"], best
        res = wres(F, self.x0)
        if res > mp.mpf(want):
            return None, "no-converge(res=%.2e)" % float(res), best
        return best, r, res


def from_cusp_record(eng, rec, verbose=False):
    """Project a cusp record onto V1 = 0 and re-solve."""
    a = mp.mpf(rec["a"]); x0 = mp.mpf(rec["x0"])
    c = CuspV1(eng, a, x0, side=rec["side"])
    v0 = [mp.mpf(rec["a11"]), mp.mpf(rec["a10"]), mp.mpf(rec["a20"])]
    v, r, res = c.newton(v0, verbose=verbose)
    out = {"a": mp.nstr(a, 30), "x0": mp.nstr(x0, 30), "side": rec["side"],
           "seed": {k: rec[k] for k in ("a11", "a01", "a10", "a20", "V1", "Dxxx")},
           "src": rec.get("_src"), "i": rec.get("_i")}
    if v is None:
        out["status"] = "FAIL:%s" % r
        return out
    a01 = c.a01_of(v[0])
    amp = c.amplitude(v)
    out.update({"status": "OK",
                "a11": mp.nstr(v[0], 34), "a01": mp.nstr(a01, 34),
                "a10": mp.nstr(v[1], 34), "a20": mp.nstr(v[2], 34),
                "res": mp.nstr(res, 8),
                "V1": mp.nstr(V1_of(a, v[0], a01), 8),
                "L": mp.nstr(L_of(a, v[2], v[0], a01, v[1]), 10),
                "D": mp.nstr(r["D"], 8), "Dx": mp.nstr(r["Dx"], 8),
                "Dxx": mp.nstr(r["Dxx"], 8), "Dxxx": mp.nstr(r["Dxxx"], 12),
                "Dxxxx": mp.nstr(r["Dxxxx"], 12) if "Dxxxx" in r else None,
                "T": mp.nstr(r["T"], 16),
                "amp": mp.nstr(amp, 8) if amp is not None else None,
                "degenerate": bool(amp is None or amp < mp.mpf("1e-24"))})
    return out
