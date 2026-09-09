"""PROTOCOL rule 7 validation of the Lane-2 engine.

(1) Cherkas rows 1-8: reproduce the cycle COUNT (3 in the nest) and positions.
    Rows 3,4,5 have their published crossings at x < 1, i.e. on the OTHER ray of
    the line y = -1 from the focus A = (1,-1); those rows are scanned with the
    left section.  The nest domain end s_max is found by walking out until the
    return fails, and never counted beyond (PROTOCOL rule 4).
(2) Row 4's Andronov-Hopf function must have 2 interior extrema on [0.6,0.9].
(3) At a fold (D = D_x = 0), D_xx must be nonzero and must agree between the
    jet value and a centred finite difference of D_x.
(4) Every bracket is certified with the two-tolerance noise estimate of rule 1.
"""
import json, sys, time
import mpmath as mp
from engine import Engine, third_order
from cusp import Cusp, solve3
from probe import Noise, sign_changes, refine_root, AH, ah_sweep

mp.mp.dps = 50

ROWS = [
    (1, "3", "-12", "-1.398", "8.4", "15.28", [1.26, 1.98, 3.95], 1),
    (2, "1.5", "-15", "0.79993", "3.2", "9.17", [1.4, 1.9, 3.1], 1),
    (3, "-2", "12", "10.999", "-14", "-26.1", [0.32, 0.66, 0.8], -1),
    (4, "-2", "-1", "9.49965", "-12.5", "6.955", [0.56, 0.75, 0.87], -1),
    (5, "-4", "-1", "13.9987", "-21", "12.4", [0.63, 0.80, 0.88], -1),
    (6, "5", "-50", "-5.49995", "16.5", "76.45", [1.05, 1.16, 1.5], 1),
    (7, "0.72727272727272727272727272727273", "-12", "2.1502",
     "0.30454545454545454545454545454545", "-26.5", [1.28, 2.15, 4.43], 1),
    (8, "1.04", "-120", "1.51997", "1.56", "-79.6", [1.29, 2.22, 4.63], 1),
]


def domain_end(c, mu, start, step, side):
    """Walk outward from the focus until the return fails; return the last good x."""
    x = mp.mpf(start)
    last = None
    for _ in range(400):
        r = c.val(mu, x)
        if r["status"] != "OK":
            break
        last = x
        x = x + step if side > 0 else x - step
        if side < 0 and x <= mp.mpf("1e-3"):
            break
    return last


def main():
    eng = Engine(quad=True)
    noise = Noise()
    print("engine:", eng.banner)
    print("noise engine:", noise.loose.banner, "(looser tolerance, rule-1 two-tolerance test)")
    rep = {"engine": eng.banner, "rows": []}
    counts_ok = 0
    for (rid, a, a20, a11, a01, a10, published, side) in ROWS:
        t0 = time.time()
        c = Cusp(eng, mp.mpf(a), mp.mpf(a20), side=side)
        mu = [mp.mpf(a11), mp.mpf(a01), mp.mpf(a10)]
        base = mp.mpf(1) + mp.mpf("0.0005") * side
        smax = domain_end(c, mu, base, mp.mpf("0.02"), side)
        if side > 0:
            lo, hi = base, smax
        else:
            lo, hi = smax, base
        xs, ds, good, weak, fails = sign_changes(c, mu, lo, hi, 500, noise)
        roots = [refine_root(c, mu, b["lo"], b["hi"]) for b in good]
        roots = sorted(float(x) for x in roots if x is not None)
        cnt_ok = len(good) == len(published)
        counts_ok += cnt_ok
        dev = max((abs(r - p) for r, p in zip(roots, sorted(published)))
                  if len(roots) == len(published) else [float("nan")])
        print("row %d (side %+d): %d certified cycles (%d uncertified) roots=%s"
              % (rid, side, len(good), len(weak), ["%.5f" % r for r in roots]))
        print("            published=%s  count_ok=%s  max|dev|=%.4f  s_max=%.4f  fails=%d  (%.0fs)"
              % (published, cnt_ok, dev, float(smax), fails, time.time() - t0))
        rep["rows"].append({"id": rid, "side": side, "roots": roots,
                            "published": published, "count_ok": bool(cnt_ok),
                            "max_dev": dev, "s_max": mp.nstr(smax, 12),
                            "n_certified": len(good), "n_uncertified": len(weak),
                            "failed_samples": fails,
                            "brackets": [{"lo": mp.nstr(b["lo"], 16), "hi": mp.nstr(b["hi"], 16),
                                          "min_abs": mp.nstr(b["min_abs"], 6),
                                          "noise": mp.nstr(b["noise"], 6)} for b in good]})
    print("\ncounts reproduced: %d/8" % counts_ok)

    # ---- (2) row 4 Andronov-Hopf function -------------------------------
    print("\n--- row 4 Andronov-Hopf function on [0.6,0.9] (published: 2 extrema) ---")
    c4 = Cusp(eng, mp.mpf("-2"), mp.mpf("-1"), side=-1)
    mu4 = [mp.mpf("9.49965"), mp.mpf("-12.5"), mp.mpf("6.955")]
    sw = ah_sweep(c4, mu4, "0.6", "0.9", 60)
    print("  AH defined at %d/61 samples, interior extrema = %d"
          % (sw["n_defined"], sw["n_interior_extrema"]))
    for e in sw["extrema"]:
        print("    %s at x=%s  a11=%s" % (e["type"], e["x"], e["a11"]))
    rep["row4_AH"] = sw

    # ---- (3) fold check --------------------------------------------------
    print("\n--- fold check: D_xx from the jet vs a centred difference of D_x ---")
    rep["fold_check"] = fold_check(eng)

    noise.close(); eng.close()
    json.dump(rep, open("validation.json", "w"), indent=1)
    print("\nengine calls:", eng.ncalls)


def fold_check(eng):
    """Newton to a fold (D = D_x = 0) in (x0, a11) on the row-4 shape."""
    c = Cusp(eng, mp.mpf("-2"), mp.mpf("-1"), side=-1)
    a01, a10 = mp.mpf("-12.5"), mp.mpf("6.955")
    x0, a11 = mp.mpf("0.80"), mp.mpf("9.49965")
    for it in range(60):
        mu = [a11, a01, a10]
        r = c.val(mu, x0)
        if r["status"] != "OK":
            return {"status": r["status"]}
        h = mp.mpf("1e-13")
        rp = c.val([a11 + h, a01, a10], x0)
        rm = c.val([a11 - h, a01, a10], x0)
        if rp["status"] != "OK" or rm["status"] != "OK":
            return {"status": "param-fd-fail"}
        Da = (rp["D"] - rm["D"]) / (2 * h)
        Dxa = (rp["Dx"] - rm["Dx"]) / (2 * h)
        dd = r["Dx"] * Dxa - Da * r["Dxx"]
        if dd == 0:
            return {"status": "singular"}
        sx = (r["D"] * Dxa - Da * r["Dx"]) / dd
        sa = (r["Dx"] * r["Dx"] - r["Dxx"] * r["D"]) / dd
        x0 -= sx; a11 -= sa
        if max(abs(sx), abs(sa)) < mp.mpf("1e-40"):
            break
    mu = [a11, a01, a10]
    r = c.val(mu, x0)
    hh = mp.mpf("1e-13")
    rp = c.val(mu, x0 + hh); rm = c.val(mu, x0 - hh)
    Dxx_fd = (rp["Dx"] - rm["Dx"]) / (2 * hh)
    rel = abs(Dxx_fd - r["Dxx"]) / abs(r["Dxx"])
    print("  fold at a11=%s x0=%s" % (mp.nstr(a11, 22), mp.nstr(x0, 22)))
    print("    D    = %s   D_x = %s" % (mp.nstr(r["D"], 6), mp.nstr(r["Dx"], 6)))
    print("    D_xx jet = %s" % mp.nstr(r["Dxx"], 22))
    print("    D_xx  fd = %s   rel diff = %.2e" % (mp.nstr(Dxx_fd, 22), float(rel)))
    print("    D_xxx    = %s   (nonzero => the fold is a simple fold, not a cusp)"
          % mp.nstr(r["Dxxx"], 12))
    return {"a11": mp.nstr(a11, 30), "x0": mp.nstr(x0, 30),
            "D": mp.nstr(r["D"], 8), "Dx": mp.nstr(r["Dx"], 8),
            "Dxx_jet": mp.nstr(r["Dxx"], 22), "Dxx_fd": mp.nstr(Dxx_fd, 22),
            "rel_diff": mp.nstr(rel, 6), "Dxxx": mp.nstr(r["Dxxx"], 12)}


if __name__ == "__main__":
    main()
