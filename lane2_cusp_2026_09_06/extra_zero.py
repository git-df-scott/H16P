"""THE CHEAP ROUTE TO FOUR IN ONE NEST.

Perko's swallow-tail is a *sufficient* local mechanism for four simple limit
cycles in one nest, not the only one.  In Cherkas's Andronov-Hopf picture, four
cycles <=> AH(x) = a11 has THREE interior extrema.  Extrema of AH are folds
(D_x = 0) and they are created in PAIRS at cusps, so parity is preserved: three
extrema can arise either by three merging at a swallow-tail, OR by a cusp
creating a pair alongside a pre-existing extremum somewhere else in the nest.

The second route is testable directly and cheaply on every cusp point already
computed:

    at a cusp point D has a TRIPLE root at x0, which is ONE sign change.
    If D has any FURTHER zero elsewhere in the nest, then perturbing into the
    cuspidal region splits the triple root into three simple cycles and the
    extra zero survives => FOUR limit cycles in one nest.

So: scan D over the whole nest at each cusp point and count certified sign
changes (PROTOCOL rule 1, two-tolerance noise estimate).  Anything with >= 2
sign changes at the cusp point is a candidate; unfold it and re-count.
"""
import glob, json, os, sys, time
import mpmath as mp
from engine import Engine
from cusp import Cusp, wres
from probe import Noise, sign_changes, refine_root, triple_confirm

mp.mp.dps = 50

AMP_FLOOR = mp.mpf("1e-24")


def nest_end(c, mu, x0, cap="200"):
    """Outer end of the nest: the last section point from which the orbit still
    returns.  Found by DOUBLING the offset from the cusp point to bracket the
    failure, then bisecting -- a linear walk costs thousands of engine calls on
    the large-amplitude curves and is the wrong tool.
    """
    sgn = 1 if c.side > 0 else -1
    r0 = (mp.mpf(x0) - 1) if c.side > 0 else (1 - mp.mpf(x0))
    good = mp.mpf(x0)
    off = r0 / 40
    bad = None
    for _ in range(60):
        x = mp.mpf(x0) + sgn * off
        if c.side < 0 and x <= mp.mpf("1e-4"):
            bad = mp.mpf("1e-4")
            break
        if c.side > 0 and x > mp.mpf(cap):
            return mp.mpf(cap)          # nest extends past the cap; report the cap
        q = c.val(mu, x)
        if q["status"] != "OK":
            bad = x
            break
        good = x
        off *= 2
    if bad is None:
        return good
    for _ in range(60):
        mid = (good + bad) / 2
        q = c.val(mu, mid)
        if q["status"] == "OK":
            good = mid
        else:
            bad = mid
        if abs(bad - good) < abs(good) * mp.mpf("1e-12"):
            break
    return good


def scan_point(c, mu, x0, noise, nsample=90):
    """Count certified sign changes of D over the WHOLE nest at a cusp point."""
    amp = c.amplitude(mu, x0)
    if amp is None or amp < AMP_FLOOR:
        return {"status": "CENTRE_VARIETY", "amp": None if amp is None else mp.nstr(amp, 6)}
    smax = nest_end(c, mu, x0)
    r0 = (mp.mpf(x0) - 1) if c.side > 0 else (1 - mp.mpf(x0))
    inner = 1 + r0 / 400 * (1 if c.side > 0 else -1)
    lo, hi = (inner, smax) if c.side > 0 else (smax, inner)
    # First pass WITHOUT the noise engine (half the calls).  Only a point that
    # already shows >= 2 sign changes is worth certifying under PROTOCOL rule 1,
    # and only then do we pay for the finer grid and the two-tolerance check.
    xs, ds, good, weak, fails = sign_changes(c, mu, lo, hi, nsample, None)
    if len(good) >= 2 and noise is not None:
        xs, ds, good, weak, fails = sign_changes(c, mu, lo, hi, nsample * 4, noise)
    roots = ([refine_root(c, mu, b["lo"], b["hi"]) for b in good]
             if len(good) >= 2 else [])
    return {"status": "OK",
            "amp": mp.nstr(amp, 6),
            "nest_end": mp.nstr(smax, 16),
            "window": [mp.nstr(lo, 14), mp.nstr(hi, 14)],
            "n_certified": len(good), "n_uncertified": len(weak), "fails": fails,
            "roots": [mp.nstr(r, 20) for r in roots if r is not None],
            "brackets": [{"lo": mp.nstr(b["lo"], 16), "hi": mp.nstr(b["hi"], 16),
                          "min_abs": mp.nstr(b["min_abs"], 8),
                          "noise": mp.nstr(b["noise"], 8)} for b in good]}


def main():
    pats = sys.argv[1:-1] or ["ledger/cusp_*.jsonl", "ledger_grid/cusp_*.jsonl"]
    outfile = sys.argv[-1] if len(sys.argv) > 1 else "extra_zero.json"
    stride = int(os.environ.get("STRIDE", "20"))
    files = []
    for p in pats:
        files.extend(sorted(glob.glob(p)))
    seeds = []
    for p in files:
        rows = [json.loads(l) for l in open(p)]
        for i, r in enumerate(rows):
            if i % stride:
                continue
            r["_src"] = os.path.basename(p)
            r["_i"] = i
            seeds.append(r)
    print("%d cusp points from %d ledgers" % (len(seeds), len(files)), flush=True)
    eng = Engine(quad=True)
    noise = Noise()
    res, hits = [], 0
    t0 = time.time()
    for k, s in enumerate(seeds):
        c = Cusp(eng, mp.mpf(s["a"]), mp.mpf(s["a20"]), side=s["side"])
        mu = [mp.mpf(s["a11"]), mp.mpf(s["a01"]), mp.mpf(s["a10"])]
        try:
            out = scan_point(c, mu, mp.mpf(s["x0"]), noise)
        except Exception as e:
            out = {"status": "EXC:%r" % e}
        out.update({"src": s["_src"], "i": s["_i"], "a": s["a"], "a20": s["a20"],
                    "x0": s["x0"], "a11": s["a11"], "a01": s["a01"], "a10": s["a10"],
                    "side": s["side"], "Dxxx": s["Dxxx"]})
        res.append(out)
        json.dump(res, open(outfile, "w"), indent=1)
        n = out.get("n_certified", -1)
        flag = ""
        if n is not None and n >= 2:
            hits += 1
            flag = "   <<<<< EXTRA ZERO IN THE NEST"
        print("%4d/%4d %-30s i=%-4d x0=%-12s nest_end=%-12s certified sign changes=%s%s"
              % (k + 1, len(seeds), out["src"][:30], out["i"], out["x0"][:10],
                 (out.get("nest_end") or "-")[:10], n, flag), flush=True)
    noise.close(); eng.close()
    print("\n%d of %d cusp points have >= 2 certified sign changes of D in the nest"
          " (%.0f s, %d engine calls)" % (hits, len(res), time.time() - t0, eng.ncalls))


if __name__ == "__main__":
    main()
