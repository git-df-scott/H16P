"""Homotopy driving V1 to zero ALONG the cusp manifold.

Unknowns u = (a11, a01, a10, a20); equations

        D = D_x = D_xx = 0,     V1 - t = 0,        V1 = a11 + a01 - 2a - 1

at fixed (a, x0).  Square (4x4).  The fourth row of the Jacobian is exact,
dV1/du = (1, 1, 0, 0) -- no integration.  Start at t = V1(seed) and walk t to
zero; a solution at t = 0 is a cusp point whose focus is simultaneously weak, so
crossing t = 0 flips sgn(V1)*sgn(D_xxx) and puts a FOURTH cycle inside the
triple (see cusp_v1.py for the sign argument).

Everything is logged, including the centre-variety amplitude guard: a homotopy
that reaches V1 = 0 by degenerating to a centre has proved nothing.
"""
import json, sys, time
import mpmath as mp
from engine import Engine, L_of, V1_of
from cusp import solve4, wres

mp.mp.dps = 50
FD = mp.mpf("1e-13")
AMP_FLOOR = mp.mpf("1e-24")


class V1Homotopy:
    def __init__(self, eng, a, x0, side=1):
        self.eng = eng
        self.a = mp.mpf(a)
        self.x0 = mp.mpf(x0)
        self.side = side

    def val(self, u, x0=None):
        return self.eng.D(self.a, u[3], u[0], u[1], u[2],
                          self.x0 if x0 is None else x0, side=self.side)

    def F(self, u, t):
        r = self.val(u)
        if r["status"] != "OK":
            return None, r
        return [r["D"], r["Dx"], r["Dxx"], V1_of(self.a, u[0], u[1]) - t], r

    def jac(self, u):
        J = [[mp.mpf(0)] * 4 for _ in range(4)]
        for j in range(4):
            h = FD * max(mp.mpf(1), abs(u[j]))
            up = list(u); up[j] += h
            um = list(u); um[j] -= h
            rp = self.val(up); rm = self.val(um)
            if rp["status"] != "OK" or rm["status"] != "OK":
                return None
            for i, k in enumerate(("D", "Dx", "Dxx")):
                J[i][j] = (rp[k] - rm[k]) / (2 * h)
        J[3] = [mp.mpf(1), mp.mpf(1), mp.mpf(0), mp.mpf(0)]   # exact dV1/du
        return J

    def ok(self, u, parmax=mp.mpf("5e3")):
        if max(abs(v) for v in u) > parmax:
            return False
        return L_of(self.a, u[3], u[0], u[1], u[2]) > 0

    def amplitude(self, u):
        r0 = (self.x0 - 1) if self.side > 0 else (1 - self.x0)
        best = None
        for fac in (mp.mpf("0.5"), mp.mpf("1.7")):
            xr = 1 + r0 * fac * (1 if self.side > 0 else -1)
            q = self.val(u, xr)
            if q["status"] == "OK":
                b = abs(q["D"])
                best = b if best is None else max(best, b)
        return best

    def newton(self, u, t, itmax=25, want="1e-24"):
        u = list(u)
        best, bestres = list(u), mp.inf
        for it in range(itmax):
            F, r = self.F(u, t)
            if F is None:
                return None, r["status"]
            res = max(wres(F[:3], self.x0), abs(F[3]))
            if res < bestres:
                bestres, best = res, list(u)
            if res < mp.mpf(want):
                return best, r
            J = self.jac(u)
            if J is None:
                return None, "jac-fail"
            s = solve4(J, F)
            if s is None:
                return None, "singular"
            sc = max(abs(s[k]) / max(mp.mpf(1), abs(u[k])) for k in range(4))
            lam = mp.mpf(1)
            if sc > mp.mpf("0.25"):
                lam = mp.mpf("0.25") / sc
            acc = None
            for _ in range(12):
                cand = [u[k] - lam * s[k] for k in range(4)]
                if self.ok(cand):
                    Fc, rc = self.F(cand, t)
                    if Fc is not None and max(wres(Fc[:3], self.x0), abs(Fc[3])) < res:
                        acc = cand
                        break
                lam /= 3
            if acc is None:
                break
            u = acc
        F, r = self.F(best, t)
        if F is not None and max(wres(F[:3], self.x0), abs(F[3])) < mp.mpf(want) * 100:
            return best, r
        return None, "no-converge(res=%.2e)" % float(bestres)


def run(eng, rec, nstep=60, verbose=True):
    a = mp.mpf(rec["a"]); x0 = mp.mpf(rec["x0"])
    h = V1Homotopy(eng, a, x0, side=rec["side"])
    u = [mp.mpf(rec["a11"]), mp.mpf(rec["a01"]), mp.mpf(rec["a10"]), mp.mpf(rec["a20"])]
    t0 = V1_of(a, u[0], u[1])
    out = {"a": mp.nstr(a, 24), "x0": mp.nstr(x0, 24), "side": rec["side"],
           "V1_start": mp.nstr(t0, 10), "src": rec.get("_src"), "i": rec.get("_i"),
           "path": [], "status": None}
    t = t0
    frac = mp.mpf(1)
    step = mp.mpf(1) / nstep
    stalls = 0
    while frac > 0:
        tgt = t0 * max(mp.mpf(0), frac - step)
        un, r = h.newton(u, tgt)
        if un is None:
            step /= 2
            stalls += 1
            if step < mp.mpf("1e-6"):
                out["status"] = ("STUCK at V1 = %s (%.3g%% of the way): %s"
                                 % (mp.nstr(t, 8), float(100 * (1 - frac)), r))
                break
            continue
        u = un
        frac = max(mp.mpf(0), frac - step)
        t = V1_of(a, u[0], u[1])
        amp = h.amplitude(u)
        rec2 = {"V1": mp.nstr(t, 10), "a11": mp.nstr(u[0], 24), "a01": mp.nstr(u[1], 24),
                "a10": mp.nstr(u[2], 24), "a20": mp.nstr(u[3], 24),
                "Dxxx": mp.nstr(r["Dxxx"], 10),
                "Dxxxx": mp.nstr(r["Dxxxx"], 10) if "Dxxxx" in r else None,
                "amp": mp.nstr(amp, 8) if amp is not None else None,
                "L": mp.nstr(L_of(a, u[3], u[0], u[1], u[2]), 8)}
        out["path"].append(rec2)
        if amp is None or amp < AMP_FLOOR:
            out["status"] = "CENTRE_VARIETY at V1 = %s (amp %s): proves nothing" % (
                mp.nstr(t, 8), rec2["amp"])
            break
        if verbose and len(out["path"]) % 10 == 0:
            print("    V1=%-16s a20=%-14s Dxxx=%-14s amp=%s"
                  % (rec2["V1"], rec2["a20"][:12], rec2["Dxxx"][:12], rec2["amp"]), flush=True)
        step = min(mp.mpf(1) / nstep, step * mp.mpf("1.5"))
        if frac == 0:
            out["status"] = "REACHED_V1_ZERO"
            out["final"] = rec2
            break
    if out["status"] is None:
        out["status"] = "INCOMPLETE"
    return out


def main():
    src = sys.argv[1]
    targets = [float(v) for v in sys.argv[2].split(",")]
    outfile = sys.argv[3]
    rows = [json.loads(l) for l in open(src)]
    eng = Engine(quad=True)
    res = []
    for tgt in targets:
        s = min(rows, key=lambda r: abs(float(r["x0"]) - tgt))
        s["_src"] = src; s["_i"] = rows.index(s)
        print("=== seed x0=%s  V1=%s  Dxxx=%s" % (s["x0"][:12], s["V1"][:12], s["Dxxx"][:12]),
              flush=True)
        o = run(eng, s)
        res.append(o)
        json.dump(res, open(outfile, "w"), indent=1)
        print("  -> %s" % o["status"], flush=True)
        if o["status"] == "REACHED_V1_ZERO":
            print("     FINAL: %s" % json.dumps(o["final"], indent=1), flush=True)
    eng.close()


if __name__ == "__main__":
    main()
