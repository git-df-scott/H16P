"""TASK 5 — the Andronov-Hopf / beta* picture at cusp points, for Lane 1.

Cherkas's Andronov-Hopf function AH(x) = a11 is defined implicitly by
D(x; a11 = AH(x)) = 0, with a11 the rotating parameter, so D is strictly
monotone in a11 (Duff/Perko) and bisection is safe.

At a cusp point (D = D_x = D_xx = 0),
    AH'  = -D_x   / D_a11 = 0
    AH'' = -D_xx  / D_a11 = 0        (after using AH' = 0)
    AH'''= -D_xxx / D_a11 != 0
so AH has a DEGENERATE INFLECTION at x0: this is exactly the point at which two
interior extrema of AH merge, i.e. the boundary between the "AH has 2 extrema"
region (3 cycles in the nest) and the "AH has 0 extrema" region (1 cycle).

A swallow-tail would be AH''' = 0 as well: three extrema merging.  Four cycles in
the nest <=> AH has THREE interior extrema, which is what Lane 1 is searching for
directly -- so Lane 1's third extremum of beta* and this lane's swallow-tail are
the SAME event, and the (3+1) route of REPORT II.7 is its other realisation.
"""
import json, sys, time
import mpmath as mp
from engine import Engine
from cusp import Cusp
from probe import AH

mp.mp.dps = 50


def sweep(c, mu, lo, hi, n, halfwidth="0.5"):
    xs, vs = [], []
    a11 = mu[0]
    for i in range(n + 1):
        x = mp.mpf(lo) + (mp.mpf(hi) - mp.mpf(lo)) * i / n
        v = AH(c, mu, x, lo=a11 - mp.mpf(halfwidth), hi=a11 + mp.mpf(halfwidth))
        if v is None:
            v = AH(c, mu, x)          # widen once
        xs.append(x); vs.append(v)
    ext = []
    for i in range(1, len(vs) - 1):
        if vs[i - 1] is None or vs[i] is None or vs[i + 1] is None:
            continue
        if (vs[i] - vs[i - 1] > 0) != (vs[i + 1] - vs[i] > 0):
            ext.append({"x": mp.nstr(xs[i], 12), "a11": mp.nstr(vs[i], 16),
                        "type": "max" if vs[i] > vs[i - 1] else "min"})
    return xs, vs, ext


def main():
    src = sys.argv[1]
    targets = [float(v) for v in sys.argv[2].split(",")]
    outfile = sys.argv[3]
    rows = [json.loads(l) for l in open(src)]
    eng = Engine(quad=True)
    out = []
    for tgt in targets:
        s = min(rows, key=lambda r: abs(float(r["x0"]) - tgt))
        c = Cusp(eng, mp.mpf(s["a"]), mp.mpf(s["a20"]), side=s["side"])
        mu = [mp.mpf(s["a11"]), mp.mpf(s["a01"]), mp.mpf(s["a10"])]
        x0 = mp.mpf(s["x0"]); r0 = x0 - 1
        lo, hi = 1 + r0 / 12, x0 + r0 * mp.mpf("0.9")
        t0 = time.time()
        xs, vs, ext = sweep(c, mu, lo, hi, 40)
        rec = {"src": src, "a": s["a"], "a20": s["a20"], "x0": s["x0"],
               "a11_cusp": s["a11"], "Dxxx": s["Dxxx"], "V1": s["V1"],
               "window": [mp.nstr(lo, 12), mp.nstr(hi, 12)],
               "n_defined": sum(1 for v in vs if v is not None),
               "n_interior_extrema": len(ext), "extrema": ext,
               "x": [mp.nstr(v, 12) for v in xs],
               "AH": [None if v is None else mp.nstr(v, 18) for v in vs],
               "wall_s": time.time() - t0}
        out.append(rec)
        json.dump(out, open(outfile, "w"), indent=1)
        print("cusp x0=%-12s : AH defined at %d/%d samples, %d interior extrema  (%.0fs)"
              % (s["x0"][:10], rec["n_defined"], len(xs), len(ext), rec["wall_s"]), flush=True)
        for e in ext:
            print("      %s at x=%s  a11=%s" % (e["type"], e["x"], e["a11"]), flush=True)
        # the cusp signature: AH should be flat to second order at x0
        i = min(range(len(xs)), key=lambda k: abs(xs[k] - x0))
        if 0 < i < len(xs) - 1 and None not in (vs[i - 1], vs[i], vs[i + 1]):
            h = xs[i + 1] - xs[i]
            d1 = (vs[i + 1] - vs[i - 1]) / (2 * h)
            d2 = (vs[i + 1] - 2 * vs[i] + vs[i - 1]) / (h * h)
            print("      near x0: AH' ~ %s   AH'' ~ %s  (both -> 0 at a cusp)"
                  % (mp.nstr(d1, 6), mp.nstr(d2, 6)), flush=True)
            rec["AH_prime_near_x0"] = mp.nstr(d1, 8)
            rec["AH_second_near_x0"] = mp.nstr(d2, 8)
            json.dump(out, open(outfile, "w"), indent=1)
    eng.close()


if __name__ == "__main__":
    main()
