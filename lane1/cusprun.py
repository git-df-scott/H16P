"""Drive the cusp continuation over random directions from each seed."""
import argparse, json, math, os, time
import numpy as np
import engine as E
import seeds as S
import ahsweep as A
import cusp as C
import sweep as W

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", default="cherkas4")
    ap.add_argument("--tag", default="c1")
    ap.add_argument("--dirs", type=int, default=80)
    ap.add_argument("--rngseed", type=int, default=99001)
    ap.add_argument("--n", type=int, default=220)
    a = ap.parse_args()

    rng = np.random.default_rng(a.rngseed)
    led = A.Ledger(os.path.join(HERE, "ledger", f"cusp_{a.tag}.jsonl"))
    tab = W.seed_table()

    for name in a.seed.split(","):
        L0 = np.ascontiguousarray(tab[name]["L"], float)
        phi = W.best_phi(L0)
        r = A.evaluate(L0, phi, n=a.n)
        feat = r[0] if isinstance(r, tuple) else r
        if feat.get("status") != "ok":
            led.write(dict(kind="cusp_seed_failed", seed=name, feat=feat))
            continue
        s_lo, s_hi = feat["s_lo"], feat["s_hi"]
        dirh = E.rotation_direction(L0, phi, math.sqrt(s_lo * s_hi))
        regs = C.regions(feat)
        led.write(dict(kind="cusp_seed", seed=name, phi=float(phi),
                       n_extrema=feat["n_extrema"], s_lo=s_lo, s_hi=s_hi,
                       regions=[dict(u_lo=r[0], u_hi=r[1], tag=r[2],
                                     starts=r[3]) for r in regs],
                       dirhint=int(dirh),
                       local10=A.coeff_strings(L0)))
        t0 = time.time()
        nsol = 0
        for k in range(a.dirs):
            v = A.project_live(L0, rng.standard_normal(10))
            nv = np.linalg.norm(v)
            if nv == 0:
                continue
            v = v / nv
            for (u_lo, u_hi, tag, starts) in regs:
              for u0 in starts:
                  sol = C.solve_cusp(L0, phi, v, u0, dirh, u_lo=u_lo, u_hi=u_hi)
                  rec = dict(kind="cusp_try", seed=name, dir=k, u0=float(u0),
                             region=tag,
                             ok=bool(sol.get("ok")), why=sol.get("why"),
                             lam=sol.get("lam"), u=sol.get("u"),
                             s=(math.exp(sol["u"]) if sol.get("u") is not None else None),
                             au=sol.get("au"), auu=sol.get("auu"),
                             auuu=sol.get("auuu"),
                             v=[repr(float(x)) for x in v])
                  if sol.get("ok"):
                      nsol += 1
                      cross = C.cross_cusp(L0, phi, v, sol["lam"], n=a.n)
                      rec["cross"] = cross
                      mx = max([c["n_extrema"] for c in cross
                                if c["n_extrema"] is not None] or [-1])
                      rec["max_extrema_near_cusp"] = mx
                      if mx >= 3:
                          for c in cross:
                              if (c["n_extrema"] or 0) >= 3:
                                  L = C.field(L0, v, c["lam"])
                                  f2 = A.evaluate(L, phi, n=400)
                                  f2 = f2[0] if isinstance(f2, tuple) else f2
                                  W.check_trigger(L, phi, f2,
                                                  f"cusp:{name}:dir{k}:lam{c['lam']:.6g}",
                                                  led)
                  led.write(rec)
            if (k + 1) % 10 == 0:
                print(f"{name} dirs={k+1}/{a.dirs} solutions={nsol} "
                      f"{time.time()-t0:.0f}s ledger={led.n}", flush=True)
        print(f"{name}: {nsol} cusp solutions from {a.dirs} directions "
              f"x {sum(len(r[3]) for r in regs)} starts, {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
