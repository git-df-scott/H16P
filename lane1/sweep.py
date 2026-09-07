"""Lane 1 sweep driver.

  python3 sweep.py perturb --reps 300
  python3 sweep.py climb --seed cherkas4 --iters 600
  python3 sweep.py climb-all --iters 400
"""
import argparse, json, math, os, time
import numpy as np
import engine as E
import seeds as S
import ahsweep as A

HERE = os.path.dirname(os.path.abspath(__file__))
LEDGER_DIR = os.path.join(HERE, "ledger")
os.makedirs(LEDGER_DIR, exist_ok=True)

PHI_CHOICES = (0.0, 0.5 * np.pi, np.pi, 1.5 * np.pi)


# --------------------------------------------------------------------- seeds
def seed_table():
    out = {}
    for rid in range(1, 9):
        sd = S.cherkas_seed(rid)
        out[f"cherkas{rid}"] = dict(L=E.local10(sd["vec12"], sd["focus"]),
                                    vec12=sd["vec12"], focus=sd["focus"])
    sd = S.perko_p3_seed()
    out["perkoP3"] = dict(L=E.local10(sd["vec12"], sd["focus"]),
                          vec12=sd["vec12"], focus=sd["focus"])
    return out


def best_phi(L, n=140, **kw):
    """Pick the ray on which the AH curve shows the most interior extrema.

    Extrema of beta* are multiple limit cycles, so the count is a property of
    the FIELD, not of the section -- but the s-interval on which the return map
    survives is very much a property of the section, and a ray whose domain is
    truncated early simply loses the extrema that lie beyond it.  Ranking by
    resolved length instead of by count therefore UNDERCOUNTS: on Cherkas row 1
    it picks a ray that reports one extremum where another ray reports two.
    Since a section can lose extrema by truncation but cannot invent them
    (spurious ones are excluded by the absolute prominence gate), the maximum
    over rays is the right estimator, with resolved length as the tie-break.
    """
    best, bphi = None, 0.0
    for phi in PHI_CHOICES:
        r = A.evaluate(L, phi, n=n, refine=False, **kw)
        f = r[0] if isinstance(r, tuple) else r
        if f.get("status") != "ok":
            continue
        key = (f["n_extrema"], f["n_run"])
        if best is None or key > best:
            best, bphi = key, phi
    return bphi


def evaluate_multi(L, phis, n=200, **kw):
    """Evaluate on several rays and keep the one that sees the most extrema."""
    best = None
    for phi in phis:
        r = A.evaluate(L, phi, n=n, **kw)
        f = r[0] if isinstance(r, tuple) else r
        if f.get("status") != "ok":
            if best is None:
                best = f
            continue
        if best is None or best.get("status") != "ok" or \
                (f["n_extrema"], f["n_run"]) > (best["n_extrema"], best["n_run"]):
            best = f
    return best


# ------------------------------------------------------------------ remote
def second_focus(L):
    """Foci of the local field other than the origin."""
    v = np.array([0.0, L[0], L[1], L[2], L[3], L[4],
                  0.0, L[5], L[6], L[7], L[8], L[9]])
    out = []
    for (p, tr, dt) in E.foci(v):
        if math.hypot(p[0], p[1]) > 1e-9:
            out.append((p, tr, dt))
    return v, out


def remote_features(L, n=120):
    v, fs = second_focus(L)
    if not fs:
        return None
    p, tr, dt = fs[0]
    LB = E.local10(v, p)
    phi = best_phi(LB, n=70)
    r = A.evaluate(LB, phi, n=n, refine=False,
                   Tmax=3000.0, Rmax=1e6)
    f = r[0] if isinstance(r, tuple) else r
    return dict(B=[float(p[0]), float(p[1])], trace=float(tr), det=float(dt),
                phi=float(phi), status=f.get("status"),
                n_extrema=f.get("n_extrema"), n_run=f.get("n_run"),
                height_range=f.get("height_range"),
                fold_margin=f.get("fold_margin"),
                score=f.get("score"), extrema=f.get("extrema"))


# ----------------------------------------------------------------- trigger
def check_trigger(L, phi, feat, tag, ledger, do_remote=True):
    """PROTOCOL rule 3.  Fires on >=3 interior extrema of beta*, or on >=4
    sign changes of D(.,b) in one nest."""
    if feat.get("n_extrema", 0) < 3:
        return None
    win = A.overlap_window(feat["extrema"])
    ts = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    rec = dict(kind="TRIGGER", tag=tag, when=ts,
               engine=E.ENGINE_NAME, engine_sha256=E.ENGINE_HASH,
               local10=A.coeff_strings(L), phi=float(phi),
               features=feat, overlap_window=win)
    if win:
        blev = 0.5 * (win[0] + win[1])
        br, s, D, st, noise = A.count_at_level(
            L, phi, blev, feat["s_lo"], feat["s_hi"], n=600)
        rec["b_level"] = float(blev)
        rec["n_sign_changes"] = len(br)
        rec["brackets"] = [[float(x) for x in q] for q in br]
        # second section, second engine
        import refengine as RF
        phi2 = phi + 0.5 * np.pi
        br2, *_ = A.count_at_level(L, phi2, blev, feat["s_lo"], feat["s_hi"], n=600)
        rec["phi2"] = float(phi2)
        rec["n_sign_changes_phi2"] = len(br2)
        ref = []
        for (s1, s2, D1, D2) in br:
            r1, k1 = RF.ret_once(L, phi, s1, blev)
            r2, k2 = RF.ret_once(L, phi, s2, blev)
            ref.append(dict(s1=float(s1), s2=float(s2),
                            D1_c=float(D1), D2_c=float(D2),
                            D1_ref=float(r1 - s1) if k1 == 0 else None,
                            D2_ref=float(r2 - s2) if k2 == 0 else None))
        rec["scipy_engine"] = ref
        try:
            import hiprec
            rec["binary128"] = hiprec.recheck(L, phi, blev, br)
        except Exception as exc:                       # never lose the trigger
            rec["binary128_error"] = repr(exc)
    if do_remote:
        rec["remote"] = remote_features(L)
    path = os.path.join(HERE, "..", f"TRIGGER_lane1_{ts}.json")
    json.dump(rec, open(path, "w"), indent=1, default=A._jd)
    ledger.write(dict(kind="trigger_written", path=os.path.basename(path), tag=tag))
    print("TRIGGER written:", path, flush=True)
    return path


# ------------------------------------------------------------------ record
def record(ledger, kind, name, L, phi, feat, extra=None):
    r = dict(kind=kind, seed=name, phi=float(phi),
             local10=A.coeff_strings(L),
             status=feat.get("status"), n_extrema=feat.get("n_extrema", -1),
             n_run=feat.get("n_run"), n=feat.get("n"),
             height_range=feat.get("height_range"),
             fold_margin=feat.get("fold_margin"),
             fold_margin_s=feat.get("fold_margin_s"),
             fold_margin_kind=feat.get("fold_margin_kind"),
             outer_turn=feat.get("outer_turn"), score=feat.get("score"),
             shoulder_curv=feat.get("shoulder_curv"),
             n_extrema_robust=feat.get("n_extrema_robust"),
             wiggle_range=feat.get("wiggle_range"),
             fold_margin_wiggle=feat.get("fold_margin_wiggle"),
             features_version=feat.get("features_version", 2),
             s_lo=feat.get("s_lo"), s_hi=feat.get("s_hi"),
             trace=feat.get("trace"), det=feat.get("det"),
             extrema=feat.get("extrema"),
             rtol=E.DEFAULTS["rtol"], btol=1e-10)
    if extra:
        r.update(extra)
    ledger.write(r)
    return r


def ev(L, phi, n):
    r = A.evaluate(L, phi, n=n)
    return r[0] if isinstance(r, tuple) else r


# ---------------------------------------------------------------- perturb
def run_perturb(args):
    rng = np.random.default_rng(args.rngseed)
    led = A.Ledger(os.path.join(LEDGER_DIR, f"perturb_{args.tag}.jsonl"))
    tab = seed_table()
    names = list(tab) if args.seed in (None, "all") else args.seed.split(",")
    t0 = time.time()
    best = {}
    for name in names:
        L0 = np.ascontiguousarray(tab[name]["L"], float)
        phi = best_phi(L0)
        f0 = ev(L0, phi, args.n)
        record(led, "seed", name, L0, phi, f0, dict(eps=0.0))
        rem = remote_features(L0)
        if rem:
            led.write(dict(kind="remote", seed=name, **rem))
        best[name] = (f0.get("n_extrema", -1), -f0.get("score", np.inf))
        norm = float(np.linalg.norm(L0))
        for eps in (1e-3, 1e-2, 1e-1):
            for k in range(args.reps):
                v = A.project_live(L0, rng.standard_normal(10))
                nv = np.linalg.norm(v)
                if nv == 0:
                    continue
                L = L0 + (eps * norm / nv) * v
                L = np.ascontiguousarray(L * (norm / np.linalg.norm(L)))
                f = ev(L, phi, args.n)
                record(led, "perturb", name, L, phi, f, dict(eps=eps, rep=k))
                if f.get("n_extrema", 0) >= 3:
                    check_trigger(L, phi, f, f"perturb:{name}:{eps}:{k}", led)
                key = (f.get("n_extrema", -1), -f.get("score", np.inf))
                if key > best[name]:
                    best[name] = key
        print(f"{name}: best (n_extrema, -fold) = {best[name]}  "
              f"ledger={led.n}  {time.time()-t0:.0f}s", flush=True)
    print("done", led.n, "rows", f"{time.time()-t0:.0f}s")


# ------------------------------------------------------------------ climb
def run_climb(args):
    rng = np.random.default_rng(args.rngseed)
    led = A.Ledger(os.path.join(LEDGER_DIR, f"climb_{args.tag}.jsonl"))
    tab = seed_table()
    names = list(tab) if args.seed in (None, "all") else args.seed.split(",")
    for name in names:
        L = np.ascontiguousarray(tab[name]["L"], float)
        norm = float(np.linalg.norm(L))
        phi = best_phi(L)
        f = ev(L, phi, args.n)
        cur = (f.get("n_extrema", -1), -f.get("score", np.inf))
        record(led, "climb_start", name, L, phi, f, dict(run=args.tag))
        sigma = args.sigma
        acc = 0
        t0 = time.time()
        for it in range(args.iters):
            v = A.project_live(L, rng.standard_normal(10))
            nv = np.linalg.norm(v)
            if nv == 0:
                continue
            Lp = L + (sigma * norm / nv) * v
            Lp = np.ascontiguousarray(Lp * (norm / np.linalg.norm(Lp)))
            fp = ev(Lp, phi, args.n)
            key = (fp.get("n_extrema", -1), -fp.get("score", np.inf))
            take = (fp.get("status") == "ok") and key > cur
            record(led, "climb", name, Lp, phi, fp,
                   dict(it=it, sigma=sigma, accepted=bool(take), run=args.tag))
            if fp.get("n_extrema", 0) >= 3:
                check_trigger(Lp, phi, fp, f"climb:{name}:{it}", led)
            if take:
                L, cur = Lp, key
                sigma = min(sigma * 1.4, 0.35)
                acc += 1
            else:
                sigma = max(sigma * 0.88, 1e-6)
            if (it + 1) % 50 == 0:
                print(f"{name} it={it+1} acc={acc} sigma={sigma:.3e} "
                      f"n_ext={cur[0]} score={-cur[1]:.5f} "
                      f"{time.time()-t0:.0f}s ledger={led.n}", flush=True)
        record(led, "climb_end", name, L, phi, ev(L, phi, args.n),
               dict(run=args.tag, accepted_total=acc))
        print(f"{name}: final n_extrema={cur[0]} score={-cur[1]:.6f}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["perturb", "climb"])
    ap.add_argument("--seed", default="all")
    ap.add_argument("--reps", type=int, default=60)
    ap.add_argument("--iters", type=int, default=400)
    ap.add_argument("--n", type=int, default=220)
    ap.add_argument("--sigma", type=float, default=0.02)
    ap.add_argument("--rngseed", type=int, default=20260906)
    ap.add_argument("--phi", type=float, default=None)
    ap.add_argument("--tag", default="r1")
    a = ap.parse_args()
    if a.mode == "perturb":
        run_perturb(a)
    else:
        run_climb(a)


if __name__ == "__main__":
    main()
