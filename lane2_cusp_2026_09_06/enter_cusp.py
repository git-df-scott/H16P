"""TASK 2 — enter the cusp manifold from the Bautin small-amplitude region.

At a third-order weak focus (V1=V3=V5=0, V7!=0) the Bautin unfolding gives
    D(r) ~ V1 r + V3 r^3 + V5 r^5 + V7 r^7
and D = V7 r (r^2 - r0^2)^3 is a TRIPLE cycle at r0 when
    V5 = -3 r0^2 V7,  V3 = 3 r0^4 V7,  V1 = -r0^6 V7.
We use the weak-focus point itself as the Newton seed (the required shift is
O(r0^2) in V5, O(r0^4) in V3, O(r0^6) in V1) and Newton on the TRUE displacement
    D = D_x = D_xx = 0
in (a11, a01, a10) at fixed x0 = 1 + r0.
"""
import json, sys, time
from engine import Engine, third_order, V7_of, L_of
from cusp import Cusp, perko_data, wres

SHAPES = [
    ("row1", 3.0, -12.0),
    ("row2", 1.5, -15.0),
    ("row3", -2.0, 12.0),
    ("row4", -2.0, -1.0),
    ("row5", -4.0, -1.0),
    ("row6", 5.0, -50.0),
    ("row7", 8.0 / 11.0, -12.0),
    ("row8", 1.04, -120.0),
]


def enter(eng, a, a20, r0, verbose=True):
    c = Cusp(eng, a, a20)
    mu0 = list(third_order(a, a20))
    x0 = 1.0 + r0
    mu, r = c.newton_mu(mu0, x0, verbose=verbose)
    return c, mu, r


def main():
    eng = Engine(quad=True)
    print("engine:", eng.banner)
    out = []
    for (name, a, a20) in SHAPES:
        mu0 = third_order(a, a20)
        v7 = V7_of(a, a20)
        L = L_of(a, a20, *mu0)
        cond = (a - 3 - a20) / (1 - 3 * a)
        print("\n=== %s  a=%.6f a20=%.6f   weak-focus mu=(%.6f,%.6f,%.6f)  V7=%.6g  L=%.4f  cond=%.4f"
              % (name, a, a20, mu0[0], mu0[1], mu0[2], v7, L, cond))
        rec = {"name": name, "a": a, "a20": a20, "mu_weakfocus": list(mu0),
               "V7_stratum": v7, "L": L, "cond_neg": cond, "entries": []}
        for r0 in (0.01, 0.02, 0.05):
            t0 = time.time()
            c, mu, r = enter(eng, a, a20, r0, verbose=False)
            if mu is None:
                print("  r0=%-6.3f  FAILED: %s" % (r0, r))
                rec["entries"].append({"r0": r0, "status": str(r)})
                continue
            M = c.jac_mu(mu, 1.0 + r0)
            pk = perko_data(M) if M else None
            res = wres([r["D"], r["Dx"], r["Dxx"]], 1.0 + r0)
            print("  r0=%-6.3f  mu=(%.15f, %.15f, %.15f)" % (r0, mu[0], mu[1], mu[2]))
            print("            res=%.2e  D_xxx=%+.8e   (48 r0^4 V7 = %+.4e)  %.1fs"
                  % (res, r["Dxxx"], 48 * r0 ** 4 * v7, time.time() - t0))
            rec["entries"].append({"r0": r0, "status": "OK", "mu": mu, "res": res,
                                   "D": r["D"], "Dx": r["Dx"], "Dxx": r["Dxx"],
                                   "Dxxx": r["Dxxx"], "T": r["T"],
                                   "dmu_shift": [mu[k] - mu0[k] for k in range(3)],
                                   "perko": pk})
        out.append(rec)
    json.dump(out, open("entry.json", "w"), indent=1)
    eng.close()
    print("\ncalls:", eng.ncalls)


if __name__ == "__main__":
    main()
