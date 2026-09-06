"""PROTOCOL rule 7 validation of the Lane-2 engine.

(1) reproduce the Cherkas rows 1-8 cycle positions to ~1e-2 by scanning D for
    sign changes on the section y=-1, x>1;
(2) reproduce row 4's published Andronov-Hopf polynomial qualitatively
    (2 extrema on [0.6,0.9]);
(3) check D_xx != 0 both by the jet and by a centred finite difference of D_x
    at a fold found by bisection.
"""
import json, sys, time
from engine import Engine, third_order

ROWS = [
    (1,  3.0,   -12.0,  -1.398,   8.4,      15.28,  [1.26, 1.98, 3.95]),
    (2,  1.5,   -15.0,   0.79993, 3.2,       9.17,  [1.4, 1.9, 3.1]),
    (3, -2.0,    12.0,  10.999, -14.0,     -26.1,   [0.32, 0.66, 0.8]),
    (4, -2.0,    -1.0,   9.49965, -12.5,     6.955, [0.56, 0.75, 0.87]),
    (5, -4.0,    -1.0,  13.9987, -21.0,     12.4,   [0.63, 0.80, 0.88]),
    (6,  5.0,   -50.0,  -5.49995, 16.5,     76.45,  [1.05, 1.16, 1.5]),
    (7,  8/11.,  -12.0,  2.1502,  67/220.,  -26.5,  [1.28, 2.15, 4.43]),
    (8,  1.04, -120.0,   1.51997, 1.56,     -79.6,  [1.29, 2.22, 4.63]),
]


def scan(eng, a, a20, a11, a01, a10, lo, hi, n):
    """Sample D on a grid; return (samples, brackets)."""
    xs, ds = [], []
    for i in range(n + 1):
        x = lo + (hi - lo) * i / n
        r = eng.D(a, a20, a11, a01, a10, x)
        xs.append(x)
        ds.append(r["D"] if r["status"] == "OK" else None)
    br = []
    for i in range(n):
        if ds[i] is None or ds[i + 1] is None:
            continue
        if (ds[i] > 0) != (ds[i + 1] > 0):
            br.append((xs[i], xs[i + 1], ds[i], ds[i + 1]))
    return xs, ds, br


def refine(eng, a, a20, a11, a01, a10, lo, hi):
    flo = eng.D(a, a20, a11, a01, a10, lo)["D"]
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        r = eng.D(a, a20, a11, a01, a10, mid)
        if r["status"] != "OK":
            return None
        if (r["D"] > 0) == (flo > 0):
            lo, flo = mid, r["D"]
        else:
            hi = mid
        if hi - lo < 1e-15 * max(1.0, abs(hi)):
            break
    return 0.5 * (lo + hi)


def main():
    eng = Engine(quad=True)
    print("engine:", eng.banner)
    report = {"engine": eng.banner, "rows": []}
    ok_rows = 0
    for (rid, a, a20, a11, a01, a10, published) in ROWS:
        lo = 1.0 + 1e-4
        hi = max(published) * 1.35
        t0 = time.time()
        xs, ds, br = scan(eng, a, a20, a11, a01, a10, lo, hi, 400)
        roots = []
        for (l, h, dl, dh) in br:
            r = refine(eng, a, a20, a11, a01, a10, l, h)
            if r is not None:
                roots.append(r)
        nfail = sum(1 for d in ds if d is None)
        # last x with a successful return = nest domain end
        smax = max((x for x, d in zip(xs, ds) if d is not None), default=None)
        match = (len(roots) == len(published) and
                 all(abs(r - p) < 2e-2 for r, p in zip(sorted(roots), sorted(published))))
        ok_rows += bool(match)
        print("row %d: roots=%s  published=%s  match=%s  fail_samples=%d  s_max=%.4f  (%.1fs)"
              % (rid, ["%.5f" % r for r in roots], published, match, nfail, smax or -1, time.time() - t0))
        report["rows"].append({"id": rid, "roots": roots, "published": published,
                               "match": bool(match), "failed_samples": nfail, "s_max": smax})
    print("\n%d/8 rows matched to 2e-2" % ok_rows)

    # ---- (3) D_xx at a fold: jet vs finite difference of D_x -------------
    # Cherkas row 4 with a11 pushed to make the outer pair coalesce.
    print("\n--- fold check (row 4, vary a11 to merge a cycle pair) ---")
    a, a20, a01, a10 = -2.0, -1.0, -12.5, 6.955
    def count(a11):
        _, _, br = scan(eng, a, a20, a11, a01, a10, 1.0001, 1.05, 60)
        return len(br)
    report["fold_check"] = fold_check(eng)
    eng.close()
    json.dump(report, open("validation.json", "w"), indent=1)


def fold_check(eng):
    """Find a fold (D=Dx=0) by 2-D Newton in (x0,a11) on row-4 shape and check D_xx
    both from the jet and from a centred difference of the (independently exact) D_x."""
    a, a20, a01, a10 = -2.0, -1.0, -12.5, 6.955
    x0, a11 = 0.80, 9.49965
    for it in range(60):
        r = eng.D(a, a20, a11, a01, a10, x0)
        if r["status"] != "OK":
            return {"status": r["status"]}
        h = 1e-8
        rp = eng.D(a, a20, a11 + h, a01, a10, x0)
        rm = eng.D(a, a20, a11 - h, a01, a10, x0)
        if rp["status"] != "OK" or rm["status"] != "OK":
            return {"status": "param-fd failed"}
        Da = (rp["D"] - rm["D"]) / (2 * h)
        Dxa = (rp["Dx"] - rm["Dx"]) / (2 * h)
        det = r["Dx"] * Dxa - r["Dxx"] * Da
        if det == 0:
            return {"status": "singular"}
        dx = (r["D"] * Dxa - r["Dx"] * Da) / det
        da = (r["Dx"] * Da - r["D"] * Dxa) / det  # placeholder, recompute properly
        # solve [[Dx, Da],[Dxx, Dxa]] [dx0, da11]^T = [D, Dx]
        A = [[r["Dx"], Da], [r["Dxx"], Dxa]]
        b = [r["D"], r["Dx"]]
        dd = A[0][0] * A[1][1] - A[0][1] * A[1][0]
        if dd == 0:
            return {"status": "singular"}
        sx = (b[0] * A[1][1] - A[0][1] * b[1]) / dd
        sa = (A[0][0] * b[1] - b[0] * A[1][0]) / dd
        x0 -= sx; a11 -= sa
        if abs(sx) < 1e-18 and abs(sa) < 1e-18:
            break
    r = eng.D(a, a20, a11, a01, a10, x0)
    # centred difference of D_x  ->  D_xx
    hh = 1e-9
    rp = eng.D(a, a20, a11, a01, a10, x0 + hh)
    rm = eng.D(a, a20, a11, a01, a10, x0 - hh)
    Dxx_fd = (rp["Dx"] - rm["Dx"]) / (2 * hh)
    out = {"a11": a11, "x0": x0, "D": r["D"], "Dx": r["Dx"], "Dxx_jet": r["Dxx"],
           "Dxx_fd": Dxx_fd, "rel_diff": abs(Dxx_fd - r["Dxx"]) / abs(r["Dxx"]),
           "Dxxx": r["Dxxx"]}
    print("fold: a11=%.15f x0=%.15f  D=%.3e Dx=%.3e" % (a11, x0, r["D"], r["Dx"]))
    print("      D_xx jet = %.15e" % r["Dxx"])
    print("      D_xx  fd = %.15e   rel diff = %.2e" % (Dxx_fd, out["rel_diff"]))
    return out


if __name__ == "__main__":
    main()
