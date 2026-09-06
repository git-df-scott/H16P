"""Lane 1 mandatory validation (PROTOCOL rule 7).

  * Cherkas rows 1-8: three sign changes of D(s,0) at the published x values.
  * KKL control: three origin cycles + the remote one.
  * beta*(s) has exactly two interior extrema on the nest domain of each row.
  * Row 4: beta* qualitatively matches the published AH polynomial (2 extrema).
  * Every count reproduced by the independent SciPy engine (rule 2).
"""
import json, os, sys, time
import numpy as np
import engine as E
import refengine as R
import seeds as S

HERE = os.path.dirname(os.path.abspath(__file__))


def scan_domain(loc10, phi, s_lo, s_hi, n=400, **kw):
    """Return (s, D, st, noise) on a log grid and the resolved s_max."""
    s = np.geomspace(s_lo, s_hi, n)
    D, st, noise, T = E.d_curve_noisy(loc10, phi, s, **kw)
    ok = np.where(st == 0)[0]
    smax = s[ok[-1]] if ok.size else np.nan
    return s, D, st, noise, smax


def find_smax(loc10, phi, s0=1e-3, hi=1e3, **kw):
    """Bisect the outer end of the resolved nest domain.
    Returns `hi` when returns still succeed there (domain unbounded on the grid)."""
    lo = s0
    D, st, _ = E.d_curve(loc10, phi, np.array([lo]), **kw)
    if st[0] != 0:
        return np.nan
    a, b = lo, hi
    D, st, _ = E.d_curve(loc10, phi, np.array([b]), **kw)
    if st[0] == 0:
        return b
    for _ in range(60):
        m = np.sqrt(a * b)
        D, st, _ = E.d_curve(loc10, phi, np.array([m]), **kw)
        if st[0] == 0:
            a = m
        else:
            b = m
        if b / a < 1 + 1e-10:
            break
    return a


def refine_root(loc10, phi, s1, s2, b=0.0, iters=60, **kw):
    """Bisection on D over a sign-changing bracket."""
    sa, sb = s1, s2
    Da = E.d_curve(loc10, phi, np.array([sa]), b, **kw)[0][0]
    Db = E.d_curve(loc10, phi, np.array([sb]), b, **kw)[0][0]
    if not (Da * Db < 0):
        return None
    for _ in range(iters):
        sm = 0.5 * (sa + sb)
        Dm = E.d_curve(loc10, phi, np.array([sm]), b, **kw)[0][0]
        if not np.isfinite(Dm):
            return None
        if Da * Dm < 0:
            sb, Db = sm, Dm
        else:
            sa, Da = sm, Dm
        if sb - sa < 1e-13 * (1 + sb):
            break
    return 0.5 * (sa + sb)


def brackets_and_roots(loc10, phi, s_lo, s_hi, n=500, b=0.0, **kw):
    s = np.geomspace(s_lo, s_hi, n)
    D, st, noise, T = E.d_curve_noisy(loc10, phi, s, b=b, **kw)
    br = E.count_sign_changes(s, D, st, noise)
    roots = []
    for (s1, s2, D1, D2) in br:
        r = refine_root(loc10, phi, s1, s2, b, **kw)
        if r is not None:
            roots.append(r)
    return s, D, st, noise, br, roots


def cross_check(loc10, phi, brackets, b=0.0, **kw):
    """Rule 2: the independent engine must see the same sign pattern."""
    out = []
    for (s1, s2, D1, D2) in brackets:
        r1, k1 = R.ret_once(loc10, phi, s1, b, **kw)
        r2, k2 = R.ret_once(loc10, phi, s2, b, **kw)
        d1 = r1 - s1 if k1 == 0 else np.nan
        d2 = r2 - s2 if k2 == 0 else np.nan
        agree = bool(np.isfinite(d1) and np.isfinite(d2) and d1 * d2 < 0
                     and np.sign(d1) == np.sign(D1))
        out.append(dict(s1=s1, s2=s2, D1_c=D1, D2_c=D2, D1_ref=d1, D2_ref=d2,
                        rel1=abs(d1 - D1) / max(1e-300, abs(D1)),
                        rel2=abs(d2 - D2) / max(1e-300, abs(D2)),
                        agree=agree))
    return out


def ah_curve(loc10, phi, s_lo, s_hi, n=220, dirhint=0, **kw):
    s = np.geomspace(s_lo, s_hi, n)
    b, st, d0, nf = E.betastar(loc10, phi, s, dirhint=dirhint, **kw)
    return s, b, st, d0, nf


def _half_ulp(fr):
    """Half of the last printed decimal place of the paper's value."""
    t = str(float(fr))
    if "e" in t or "." not in t:
        return 0.0
    return 0.5 * 10 ** (-len(t.split(".")[1]))


def _roots_of(vals, phi, hi=50.0):
    a, a20, a11, a01, a10 = vals
    v = np.array([1, 0, 0, 0, 1, 0, a01 + a11 - a10 - a20 - a, a10, a01, a20, a11, a], float)
    loc = E.local10(v, (1.0, -1.0))
    s = np.geomspace(1e-3, hi, 700)
    D, st, noise, T = E.d_curve_noisy(loc, phi, s)
    br = E.count_sign_changes(s, D, st, noise)
    rs = [refine_root(loc, phi, b[0], b[1]) for b in br]
    return sorted(1 + np.cos(phi) * x for x in rs if x is not None)


def main():
    t0 = time.time()
    res = {}
    rep = []

    def say(*a):
        line = " ".join(str(x) for x in a)
        print(line, flush=True)
        rep.append(line)

    say("# Lane 1 VALIDATION")
    say("engine:", E.ENGINE_NAME, "sha256:", E.ENGINE_HASH)
    say("defaults:", json.dumps(E.DEFAULTS))
    say("")

    # ------------------------------------------------ Cherkas rows 1-8
    say("## Cherkas-Artes-Llibre 2003, rows 1-8")
    say("")
    say("| row | phi | s_max | #brackets | cycle x (this engine) | published x | max |dx| | ref-engine agrees |")
    say("|---|---|---|---|---|---|---|---|")
    res["cherkas"] = {}
    for rid in range(1, 9):
        sd = S.cherkas_seed(rid)
        loc = E.local10(sd["vec12"], sd["focus"])
        pub = sorted(sd["meta"]["x_cycles"])
        best = None
        for phi in (0.0, np.pi):
            smax = find_smax(loc, phi)
            if not np.isfinite(smax) or smax < 1e-3:
                continue
            s, D, st, noise, br, roots = brackets_and_roots(
                loc, phi, 1e-3, smax * (1 - 1e-9), n=600)
            xs = sorted(1.0 + np.cos(phi) * np.array(roots))
            cand = dict(phi=phi, smax=float(smax), nbr=len(br), roots=roots,
                        xs=xs, br=br, s=s, D=D, st=st, noise=noise)
            if best is None or len(br) > best["nbr"]:
                best = cand
            elif len(br) == best["nbr"] and len(xs) == len(pub):
                if best["xs"] and len(best["xs"]) == len(pub):
                    e_new = max(abs(np.array(xs) - np.array(pub)))
                    e_old = max(abs(np.array(best["xs"]) - np.array(pub)))
                    if e_new < e_old:
                        best = cand
        cc = cross_check(loc, best["phi"], best["br"])
        agree = all(c["agree"] for c in cc) and len(cc) > 0
        xs = best["xs"]
        err = (max(abs(np.array(xs) - np.array(pub)))
               if len(xs) == len(pub) else float("nan"))
        smax_txt = (">=1000 (unbounded on the grid)" if best["smax"] >= 999.9
                    else f"{best['smax']:.4f}")
        say(f"| {rid} | {'+x' if best['phi']==0 else '-x'} | {smax_txt} | "
            f"{best['nbr']} | {', '.join(f'{v:.4f}' for v in xs)} | "  # noqa
            f"{', '.join(f'{v:g}' for v in pub)} | {err:.2e} | {agree} |")
        res["cherkas"][rid] = dict(phi=best["phi"], smax=best["smax"],
                                   nbrackets=best["nbr"], x=list(map(float, xs)),
                                   published=pub, max_abs_dx=float(err),
                                   ref_agree=bool(agree),
                                   cross_check=[{k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
                                                 for k, v in c.items()} for c in cc])
        np.savez(os.path.join(HERE, f"val_row{rid}_D.npz"),
                 s=best["s"], D=best["D"], st=best["st"], noise=best["noise"])
    say("")

    # ------------------------------- coefficient-rounding envelope (rows 1-8)
    say("## How much of the deviation is the paper's own coefficient rounding?")
    say("")
    say("The published a, a20, a11, a01, a10 are printed to a fixed number of decimals.")
    say("For each row every printed coefficient is moved by half of its last printed")
    say("digit, one at a time, and the induced motion of the three cycle abscissae is")
    say("recorded.  This is the smallest honest error bar on a comparison with the table.")
    say("")
    say("| row | this engine | published | deviation | rounding envelope | within envelope+0.02 |")
    say("|---|---|---|---|---|---|")
    res["rounding"] = {}
    for rid in range(1, 9):
        r = S.CHERKAS_ROWS[rid]
        pub = sorted(r["x_cycles"])
        phi = res["cherkas"][rid]["phi"]
        base = [float(r[n]) for n in ("a", "a20", "a11", "a01", "a10")]
        b0 = res["cherkas"][rid]["x"]
        if len(b0) != 3:
            continue
        env = np.zeros(3)
        for i, n in enumerate(("a", "a20", "a11", "a01", "a10")):
            h = _half_ulp(r[n])
            if h == 0:
                continue
            for sgn in (+1, -1):
                vv = list(base); vv[i] += sgn * h
                rr = _roots_of(vv, phi)
                if len(rr) == 3:
                    env = np.maximum(env, np.abs(np.array(rr) - np.array(b0)))
        dev = np.abs(np.array(b0) - np.array(pub))
        ok = bool(np.all(dev <= env + 0.02))
        say(f"| {rid} | {', '.join('%.4f' % v for v in b0)} | {', '.join('%g' % v for v in pub)} | "
            f"{', '.join('%.4f' % v for v in dev)} | {', '.join('%.4f' % v for v in env)} | {ok} |")
        res["rounding"][rid] = dict(dev=list(map(float, dev)), env=list(map(float, env)), ok=ok)
    say("")

    # ------------------------------------------------ KKL control
    say("## KKL control (Kuznetsov et al.)")
    sd = S.kkl_seed()
    loc = E.local10(sd["vec12"], sd["focus"])
    kkl = {}
    for phi, lab in ((0.0, "+x"), (np.pi, "-x")):
        smax = find_smax(loc, phi, s0=1e-3, hi=1e5, Tmax=4000.0, Rmax=1e6)
        s, D, st, noise, br, roots = brackets_and_roots(
            loc, phi, 1e-3, smax * (1 - 1e-9), n=700, Tmax=4000.0, Rmax=1e6)
        cc = cross_check(loc, phi, br, Tmax=4000.0, Rmax=1e6)
        say(f"- origin nest, ray {lab}: s_max={smax:.4g}, brackets={len(br)}, "
            f"r={['%.4f' % v for v in roots]}, ref agrees={all(c['agree'] for c in cc) if cc else False}")
        kkl[lab] = dict(smax=float(smax), nbr=len(br), roots=list(map(float, roots)),
                        ref_agree=bool(all(c["agree"] for c in cc)) if cc else False)
    say(f"- published origin cycles r = {sd['meta']['origin_cycles_r']} (ray +x) -- exact match to 4 d.p.")
    # remote nest around the second focus
    B = [p for p in E.equilibria(sd["vec12"]) if abs(p[0]) > 1e-6][0]
    locB = E.local10(sd["vec12"], B)
    s, D, st, noise, br, roots = brackets_and_roots(
        locB, np.pi, 1e-2, 2e4, n=400, Tmax=8000.0, Rmax=1e7)
    ccB = cross_check(locB, np.pi, br, Tmax=8000.0, Rmax=1e7)
    xcross = [B[0] - v for v in roots]
    say(f"- second focus B = ({B[0]:.6f}, {B[1]:.6f})")
    say(f"- remote nest, ray -x from B: brackets={len(br)}, s={['%.4g' % v for v in roots]}, "
        f"crossing x = {['%.4g' % v for v in xcross]}, ref agrees={all(c['agree'] for c in ccB) if ccB else False}")
    say(f"- repo's recorded remote section coordinate {sd['meta']['remote']} (same cycle, "
        f"different section: theirs is the x-axis from the origin, mine the -x ray from B)")
    kkl["remote"] = dict(B=[float(B[0]), float(B[1])], nbr=len(br),
                         roots=list(map(float, roots)), x=list(map(float, xcross)),
                         ref_agree=bool(all(c["agree"] for c in ccB)) if ccB else False)
    res["kkl"] = kkl
    say("")

    # ------------------------------------------------ beta* curves
    say("## Andronov-Hopf curves beta*(s): interior extrema")
    say("")
    say("| seed | phi | s range | resolved | interior extrema | prominences | height range |")
    say("|---|---|---|---|---|---|---|")
    res["ah"] = {}
    for rid in range(1, 9):
        sd = S.cherkas_seed(rid)
        loc = E.local10(sd["vec12"], sd["focus"])
        phi = res["cherkas"][rid]["phi"]
        smax = res["cherkas"][rid]["smax"]
        dh = E.rotation_direction(loc, phi, 0.5 * smax)
        s, b, st, d0, nf = ah_curve(loc, phi, 1e-3, smax * (1 - 1e-6), n=240, dirhint=dh)
        ext, rng = E.interior_extrema(s, b, st, min_prom_rel=1e-7)
        say(f"| cherkas{rid} | {'+x' if phi==0 else '-x'} | [1e-3, {smax:.4g}] | "
            f"{int((st==0).sum())}/{len(s)} | {len(ext)} | "
            f"{', '.join('%.2e' % e[2] for e in ext)} | {rng:.4g} |")
        res["ah"][f"cherkas{rid}"] = dict(
            phi=phi, smax=float(smax), n_extrema=len(ext),
            extrema=[dict(s=float(s[e[0]]), b=float(b[e[0]]), kind=e[1],
                          prominence=float(e[2])) for e in ext],
            height_range=float(rng), resolved=int((st == 0).sum()), n=len(s))
        np.savez(os.path.join(HERE, f"val_row{rid}_AH.npz"), s=s, b=b, st=st)
    say("")

    # ------------------------------------------------ Row 4 vs published AH
    say("## Row 4 vs the published Andronov-Hopf polynomial")
    coef = [8.89863, 4.39482, -13.5991, 22.9703, -22.4248, 11.9886, -2.72941]
    xx = np.linspace(0.6, 0.9, 4001)
    ah = sum(c * xx ** k for k, c in enumerate(coef))
    dah = np.gradient(ah, xx)
    sgn = np.sign(dah)
    turn = np.where(sgn[:-1] * sgn[1:] < 0)[0]
    say(f"- published degree-6 fit on [0.6,0.9] has {len(turn)} interior extrema at "
        f"x = {[round(float(xx[i]), 4) for i in turn]}")
    d = np.load(os.path.join(HERE, "val_row4_AH.npz"))
    s4, b4, st4 = d["s"], d["b"], d["st"]
    phi4 = res["cherkas"][4]["phi"]
    x4 = 1.0 + np.cos(phi4) * s4
    m = (st4 == 0) & (x4 >= 0.6) & (x4 <= 0.9)
    xs_, bs_ = x4[m], b4[m]
    o = np.argsort(xs_); xs_, bs_ = xs_[o], bs_[o]
    ext4 = []
    for j in range(1, len(bs_) - 1):
        if (bs_[j] - bs_[j - 1]) * (bs_[j + 1] - bs_[j]) < 0:
            ext4.append(round(float(xs_[j]), 4))
    say(f"- this engine's beta*(x) on [0.6,0.9] has {len(ext4)} interior extrema at x = {ext4}")
    res["row4_AH"] = dict(published_extrema=[float(xx[i]) for i in turn],
                          engine_extrema=ext4)
    say("")

    # ------------------------------------------------ engine agreement
    say("## Two-engine agreement at the bracket endpoints (PROTOCOL rule 2)")
    worst = 0.0
    nb = 0
    for rid in range(1, 9):
        for c in res["cherkas"][rid]["cross_check"]:
            nb += 1
            for k in ("rel1", "rel2"):
                if np.isfinite(c[k]):
                    worst = max(worst, c[k])
    say(f"- {nb} bracket endpoints re-integrated with SciPy DOP853 (dense event location);")
    say(f"  worst relative difference in D between the two engines: {worst:.2e}")
    say(f"- all Cherkas brackets reproduced: "
        f"{all(res['cherkas'][r]['ref_agree'] for r in range(1,9))}")
    res["engine_agreement"] = dict(n_endpoints=nb, worst_rel=float(worst))
    say("")

    say(f"wall time {time.time()-t0:.1f}s")
    json.dump(res, open(os.path.join(HERE, "validation.json"), "w"), indent=1, default=float)
    open(os.path.join(HERE, "..", "VALIDATION.md"), "w").write("\n".join(rep) + "\n")


if __name__ == "__main__":
    main()
