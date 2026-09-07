"""Sweep the cusp manifold of triple limit cycles and read beta* on it.

For each (a, a20) in a grid and each target amplitude r0, `bautin.py` solves the
three Bautin conditions exactly, so every field produced here carries a limit
cycle of multiplicity THREE by construction -- a cusp of the Andronov-Hopf
curve -- rather than merely being near one.

What is being looked for is not the cusp itself.  A cusp at the focus end gives
two extrema of beta*, i.e. the three small cycles Bautin allows and no more.
The question is whether a field ON that manifold ALSO carries structure further
out: an extremum of beta* beyond the cusp is a third extremum, and a third
extremum is four cycles in the nest.  Bautin caps what comes out of the focus,
not what the rest of the nest does.
"""
import argparse, math, os, time
import numpy as np
import bautin as B
import engine as E
import ahsweep as A
import sweep as W

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="b1")
    ap.add_argument("--na", type=int, default=25)
    ap.add_argument("--na20", type=int, default=25)
    ap.add_argument("--a-lo", type=float, default=-6.0)
    ap.add_argument("--a-hi", type=float, default=6.0)
    ap.add_argument("--a20-lo", type=float, default=-130.0)
    ap.add_argument("--a20-hi", type=float, default=30.0)
    ap.add_argument("--r0", default="0.02,0.05,0.1,0.2,0.4,0.7,1.0,1.5")
    ap.add_argument("--n", type=int, default=220)
    ap.add_argument("--unfold", default="0,1e-6,-1e-6,1e-4,-1e-4,1e-2,-1e-2,0.3,-0.3")
    a = ap.parse_args()

    led = A.Ledger(os.path.join(HERE, "ledger", f"bautin_{a.tag}.jsonl"))
    r0s = [float(x) for x in a.r0.split(",")]
    eps = [float(x) for x in a.unfold.split(",")]
    avals = np.linspace(a.a_lo, a.a_hi, a.na)
    a20vals = np.linspace(a.a20_lo, a.a20_hi, a.na20)

    t0 = time.time()
    nsol = nfocus = nok = 0
    best = (2, None)
    for av in avals:
        if abs(av - 2.0) < 1e-9 or abs(1 - 3 * av) < 1e-9:
            continue                      # the family's own excluded values
        for a20v in a20vals:
            for r0 in r0s:
                try:
                    sols = B.triple_cycle_field(av, a20v, r0)
                except Exception:
                    continue
                for p0 in sols:
                  for (e1, e3) in [(x, y) for x in eps for y in eps]:
                    p = B.unfold(p0, e1, e3) if (e1 or e3) else p0
                    nsol += 1
                    v = B.vec12(p)
                    J = E.jac(v, (1.0, -1.0))
                    tr, dt = J[0, 0] + J[1, 1], float(np.linalg.det(J))
                    if not (dt > 0 and tr * tr < 4 * dt):
                        continue          # A is a saddle or a node, not a focus
                    nfocus += 1
                    L = E.local10(v, (1.0, -1.0))
                    phi = W.best_phi(L, n=110, Tmax=100.0, nstep=400_000)
                    f = W.evaluate_multi(L, [phi, phi + 0.5 * math.pi], n=a.n,
                                         Tmax=100.0, nstep=400_000)
                    if f.get("status") == "ok":
                        nok += 1
                    rec = dict(kind="bautin", a=float(av), a20=float(a20v),
                               r0=float(r0), a11=p["a11"], a01=p["a01"],
                               a10=p["a10"], V=p["V"],
                               e1=p.get("e1", 0.0), e3=p.get("e3", 0.0),
                               cusp_residual=p0["residual"],
                               trace=float(tr), det=dt, phi=float(phi),
                               local10=A.coeff_strings(L),
                               vec12=[repr(float(x)) for x in v],
                               status=f.get("status"),
                               n_extrema=f.get("n_extrema"),
                               n_extrema_robust=f.get("n_extrema_robust"),
                               n_run=f.get("n_run"),
                               wiggle_range=f.get("wiggle_range"),
                               fold_margin=f.get("fold_margin"),
                               fold_margin_kind=f.get("fold_margin_kind"),
                               extrema=f.get("extrema"))
                    led.write(rec)
                    ne = f.get("n_extrema", -1)
                    if ne > best[0]:
                        best = (ne, rec)
                    if ne >= 3:
                        f2 = W.evaluate_multi(L, [phi, phi + 0.5 * math.pi], n=450)
                        if f2.get("n_extrema", 0) >= 3:
                            W.check_trigger(L, f2.get("phi", phi), f2,
                                            f"bautin:a{av:.4g}:a20{a20v:.4g}:r0{r0:g}"
                                            f":e{p.get('e1',0):g},{p.get('e3',0):g}",
                                            led)
        print(f"a={av:+.3f}  solutions={nsol} foci={nfocus} resolved={nok} "
              f"best n_extrema={best[0]}  {time.time()-t0:.0f}s ledger={led.n}",
              flush=True)
    print(f"done: {nsol} cusp fields, {nfocus} with a focus at A, {nok} resolved, "
          f"max interior extrema {best[0]}, {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
