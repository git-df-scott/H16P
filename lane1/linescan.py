"""One-dimensional scans of the interior-extremum count along straight lines
through the live coefficient space.

The (1+1)-ES climbs of the first campaign stalled on step-size collapse, and a
Newton solve for the cusp only finds what it is started near.  A line scan
cannot do either: for each random direction v in the 8 live directions it walks
lambda across a fixed grid and counts the interior extrema of beta* at every
point, so the whole segment is seen, including the far end where the field no
longer resembles the seed.

Recorded per direction:
  * lam_break        -- the smallest |lambda| at which the count drops below 2,
                        i.e. how far the seed's three-cycle configuration
                        survives along that line;
  * max n_extrema    -- and every lambda that reaches 3 or more (a TRIGGER);
  * the best fold margin seen and where.
"""
import argparse, json, math, os, time
import numpy as np
import engine as E
import ahsweep as A
import sweep as W

HERE = os.path.dirname(os.path.abspath(__file__))


def scan(L0, phi, v, lams, n=220):
    out = []
    for lam in lams:
        n0 = float(np.linalg.norm(L0))
        L = L0 + lam * n0 * v
        nn = float(np.linalg.norm(L))
        if nn == 0:
            continue
        L = np.ascontiguousarray(L * (n0 / nn))
        r = A.evaluate(L, phi, n=n)
        f = r[0] if isinstance(r, tuple) else r
        out.append((float(lam), L, f))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", default="cherkas4")
    ap.add_argument("--tag", default="s1")
    ap.add_argument("--dirs", type=int, default=200)
    ap.add_argument("--rngseed", type=int, default=41001)
    ap.add_argument("--n", type=int, default=220)
    ap.add_argument("--lam-max", type=float, default=0.6)
    ap.add_argument("--nlam", type=int, default=25)
    a = ap.parse_args()

    rng = np.random.default_rng(a.rngseed)
    led = A.Ledger(os.path.join(HERE, "ledger", f"line_{a.tag}.jsonl"))
    tab = W.seed_table()

    # geometric spacing on each side: the interesting scale spans 1e-4 to 1
    half = np.geomspace(1e-4, a.lam_max, a.nlam)
    lams = np.concatenate([-half[::-1], [0.0], half])

    for name in a.seed.split(","):
        L0 = np.ascontiguousarray(tab[name]["L"], float)
        phi = W.best_phi(L0)
        t0 = time.time()
        led.write(dict(kind="line_seed", seed=name, phi=float(phi),
                       lam_grid=[float(x) for x in lams],
                       local10=A.coeff_strings(L0)))
        best_max, n3 = 2, 0
        for k in range(a.dirs):
            v = A.project_live(L0, rng.standard_normal(10))
            nv = np.linalg.norm(v)
            if nv == 0:
                continue
            v = v / nv
            rows = scan(L0, phi, v, lams, n=a.n)
            counts = [(lam, f.get("n_extrema", -1), f.get("n_extrema_robust", -1),
                       f.get("score"), f.get("fold_margin"), f.get("status"))
                      for (lam, L, f) in rows]
            ok2 = [lam for (lam, c, _cr, _s, _fm, st) in counts
                   if st == "ok" and c >= 2]
            pos = [lam for lam in ok2 if lam > 0]
            neg = [lam for lam in ok2 if lam < 0]
            mx = max([c for (_l, c, _cr, _s, _fm, _st) in counts] or [-1])
            best_max = max(best_max, mx)
            rec = dict(kind="line", seed=name, dir=k,
                       v=[repr(float(x)) for x in v],
                       counts=[[l, c, cr, s, fm, st]
                               for (l, c, cr, s, fm, st) in counts],
                       max_extrema=mx,
                       lam_break_pos=(max(pos) if pos else 0.0),
                       lam_break_neg=(min(neg) if neg else 0.0))
            led.write(rec)
            if mx >= 3:
                n3 += 1
                for (lam, L, f) in rows:
                    if f.get("n_extrema", 0) >= 3:
                        f2 = A.evaluate(L, phi, n=450)
                        f2 = f2[0] if isinstance(f2, tuple) else f2
                        if f2.get("n_extrema", 0) >= 3:
                            W.check_trigger(L, phi, f2,
                                            f"line:{name}:dir{k}:lam{lam:.6g}", led)
            if (k + 1) % 20 == 0:
                print(f"{name} dirs={k+1}/{a.dirs} max_extrema={best_max} "
                      f"n3={n3} {time.time()-t0:.0f}s ledger={led.n}", flush=True)
        print(f"{name}: {a.dirs} directions x {len(lams)} lambda, "
              f"max extrema {best_max}, {n3} directions reaching 3, "
              f"{time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
