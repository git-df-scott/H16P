"""lane1_ahcurve/validate.py -- mandatory engine validation (PROTOCOL rule 7).

Checks, for the nine fat seeds and the KKL control:
  V1  three sign changes of D in the primary nest at the published parameters,
      each clearing the two-tolerance noise floor (PROTOCOL rule 1);
  V2  cycle positions against the published values;
  V3  the same count on a second, rotated section (PROTOCOL rule 4);
  V4  every bracket endpoint reproduced by scipy DOP853 in global coordinates
      (PROTOCOL rule 2);
  V5  beta*(s) for the uniform rotation has exactly 2 interior extrema;
  V6  Cherkas rows only: the a11 Andronov-Hopf curve evaluated at the published
      x-values returns the published a11, and for row 4 the whole curve is
      compared with the published degree-6 polynomial;
  V7  remote nest count where a second focus exists (the "+1" of a (3,1)).
"""
import json, math, time, sys
import numpy as np
import engine, seeds, refengine, ledger

RTOL = 1e-12
RTOL_LOOSE = 1e-10
NGRID = 300
PTOL = 1e-10

def rot_dir(d, th):
    c, s = math.cos(th), math.sin(th)
    return (d[0]*c - d[1]*s, d[0]*s + d[1]*c)

def probe_smax(L, d, s_lo, s_hi, b):
    """Largest probed s with a successful return, on a 90-point geometric grid."""
    s = np.geomspace(s_lo, s_hi, 90)
    R, T, st = engine.returns(L, d, s, b=b, rtol=RTOL)
    ok = st == 0
    if not ok.any():
        return None, None
    i = int(np.max(np.nonzero(ok)[0]))
    return float(s[i]), (float(s[i+1]) if i+1 < s.size else None)

def refine(L, d, brs, b, tol=1e-12):
    out = []
    for (a, c) in brs:
        for _ in range(70):
            m = 0.5*(a+c)
            D, st = engine.displacement(L, d, [a, m], b=b, rtol=RTOL)
            if st[1] != 0:
                break
            if D[0]*D[1] < 0:
                c = m
            else:
                a = m
            if c-a < tol*max(1.0, c):
                break
        out.append(0.5*(a+c))
    return out

def nest_count(L, d, s_lo, s_hi, b):
    s = np.geomspace(s_lo, s_hi, NGRID)
    D, st, nz = engine.displacement_with_noise(L, d, s, b=b, rtol=RTOL, rtol_loose=RTOL_LOOSE)
    brs = engine.count_sign_changes(s, D, nz)
    return s, D, st, nz, brs

def beta_curve(L, d, s_lo, s_hi, n=NGRID, p0=None, span=1.2):
    s = np.geomspace(s_lo, s_hi, n)
    b, st, nz, nev = engine.curve_with_noise(L, d, s, mode="rot", p0=p0,
                                             ptol=PTOL, span=span,
                                             rtol=RTOL, rtol_loose=RTOL_LOOSE)
    return s, b, st, nz, nev

def extrema_of(b, nz, margin=8.0):
    good = np.isfinite(b) & np.isfinite(nz)
    if good.sum() < 5:
        return [], float("nan"), float("nan")
    tau = margin*float(np.max(nz[good]))
    ext = engine.turning_points(np.where(good, b, np.nan), tau)
    rng = float(np.nanmax(b[good]) - np.nanmin(b[good]))
    return ext, tau, rng

def run(seed, out):
    t0 = time.time()
    L = [float(v) for v in seed["local"]]
    G = [float(v) for v in seed["vec12"]]
    F = [float(v) for v in seed["focus"]]
    d = seed["direction"]
    b0 = seed["base_b"]
    sc = sorted(seed["s_cycles"])
    s_lo, s_hi_probe = 0.25*sc[0], 4.0*sc[-1]
    smax, sfail = probe_smax(L, d, s_lo, s_hi_probe, b0)
    s_hi = min(0.999*smax, s_hi_probe) if smax else s_hi_probe

    s, D, st, nz, brs = nest_count(L, d, s_lo, s_hi, b0)
    roots = refine(L, d, brs, b0)

    # V3 second section
    alt = {}
    for th in (0.7, -0.7):
        dd = rot_dir(d, th)
        sm2, _ = probe_smax(L, dd, s_lo, s_hi_probe, b0)
        hi2 = min(0.999*sm2, s_hi_probe) if sm2 else s_hi_probe
        s2, D2, st2, nz2, br2 = nest_count(L, dd, s_lo, hi2, b0)
        alt["theta=%+.2f" % th] = dict(n=len(br2), s_max=sm2,
                                       roots=refine(L, dd, br2, b0))

    # V4 second integrator at every bracket endpoint
    ref = []
    for (a, c) in brs:
        for sv in (a, c):
            De = float(engine.displacement(L, d, [sv], b=b0, rtol=RTOL)[0][0])
            Rr = refengine.ret_radius(G, F, d, sv, b=b0, rtol=RTOL)
            Dr = float("nan") if Rr is None else Rr - sv
            ref.append(dict(s=sv, D_engine=De, D_scipy=Dr,
                            absdiff=abs(De-Dr) if Rr is not None else None,
                            sign_agree=bool(Rr is not None and De*Dr > 0)))

    # V5 beta*
    sb, bb, bst, bnz, bnev = beta_curve(L, d, s_lo, s_hi)
    ext, tau, rng = extrema_of(bb, bnz)

    rec = dict(kind="validation", engine=engine.ENGINE, seed=seed["name"],
               family=seed["family"], row=seed["row"],
               vec12_exact=[str(v) for v in seed["vec12"]],
               focus=[str(v) for v in seed["focus"]],
               direction=list(d), base_b=b0, params=seed["params"],
               s_grid=dict(lo=s_lo, hi=s_hi, n=NGRID),
               s_max=smax, s_first_fail=sfail,
               n_sign_changes=len(brs), brackets=[list(x) for x in brs],
               cycle_s=roots, published_s=seed["s_cycles"],
               second_section=alt, second_integrator=ref,
               beta_extrema=[[int(i), int(g)] for i, g in ext],
               beta_extrema_s=[float(sb[i]) for i, _ in ext],
               beta_tau=tau, beta_range=rng,
               beta_status_ok=int((bst == 0).sum()), beta_n=int(bst.size),
               beta_returns=int(bnev),
               rtol=RTOL, rtol_loose=RTOL_LOOSE, ptol=PTOL,
               wall_s=round(time.time()-t0, 3))

    if seed["family"] == "cherkas":
        E = [float(v) for v in seed["evec_local"]]
        a11 = float(seed["params"]["a11"])
        xs = np.array(seed["x_cycles"], float)
        sv = np.abs(xs - 1.0)
        t, tst, tnz, _ = engine.curve_with_noise(L, d, sv, mode="lin", evec=E,
                                                 p0=0.0, ptol=PTOL, span=2.0,
                                                 rtol=RTOL, rtol_loose=RTOL_LOOSE)
        sg = np.geomspace(s_lo, s_hi, NGRID)
        tg, tgs, tgn, _ = engine.curve_with_noise(L, d, sg, mode="lin", evec=E,
                                                  p0=0.0, ptol=PTOL, span=2.0,
                                                  rtol=RTOL, rtol_loose=RTOL_LOOSE)
        gext, gtau, grng = extrema_of(tg, tgn)
        rec["ah_a11_at_published_x"] = [a11 + float(v) for v in t]
        rec["ah_a11_published"] = a11
        rec["ah_a11_max_dev"] = float(np.nanmax(np.abs(t)))
        rec["ah_a11_range_over_nest"] = grng
        rec["ah_a11_extrema_s"] = [float(sg[i]) for i, _ in gext]
        rec["ah_a11_extrema_sign"] = [int(g) for _, g in gext]
        if seed["row"] == 4:
            xg = np.linspace(0.6, 0.9, 61)
            sp = 1.0 - xg
            tp, _, npz, _ = engine.curve_with_noise(L, d, sp, mode="lin", evec=E,
                                                    p0=0.0, ptol=PTOL, span=2.0,
                                                    rtol=RTOL, rtol_loose=RTOL_LOOSE)
            AH = a11 + tp
            co = [8.89863, 4.39482, -13.5991, 22.9703, -22.4248, 11.9886, -2.72941]
            pub = sum(c*xg**k for k, c in enumerate(co))
            pe = engine.turning_points(pub, 1e-13)
            ee = engine.turning_points(AH, 8*float(np.nanmax(npz)))
            rec["row4"] = dict(
                x=[float(v) for v in xg], AH_engine=[float(v) for v in AH],
                AH_published=[float(v) for v in pub],
                max_abs_dev=float(np.nanmax(np.abs(AH-pub))),
                engine_range=float(AH.max()-AH.min()),
                published_range=float(pub.max()-pub.min()),
                engine_extrema_x=[float(xg[i]) for i, _ in ee],
                published_extrema_x=[float(xg[i]) for i, _ in pe])

    # V7 remote nest
    eqs = engine.equilibria(G)
    rem = []
    for (x, y) in eqs:
        if abs(x-F[0]) < 1e-7 and abs(y-F[1]) < 1e-7:
            continue
        Lr = engine.local_expand(G, x, y)
        if not engine.is_focus(Lr):
            continue
        best = None
        for dd in ((1.0, 0.0), (-1.0, 0.0), (0.0, 1.0), (0.0, -1.0)):
            sm, _ = probe_smax(Lr, dd, 1e-4, 1e5, b0)
            if sm is None:
                continue
            s2, D2, st2, nz2, br2 = nest_count(Lr, dd, 1e-4, 0.999*sm, b0)
            if best is None or len(br2) > best["n"]:
                best = dict(dir=list(dd), n=len(br2), s_max=sm,
                            roots=refine(Lr, dd, br2, b0))
        rem.append(dict(focus=[x, y], nest=best))
    rec["remote_foci"] = rem
    ledger.append("validation", rec)
    out.append(rec)
    return rec

def main():
    out = []
    allseeds = seeds.fat_seeds() + [seeds.kkl_control()]
    for sd in allseeds:
        r = run(sd, out)
        print("%-12s n=%d  s_max=%-10.4g  cycles=%s  beta*extrema=%d  wall=%.1fs"
              % (r["seed"], r["n_sign_changes"], r["s_max"] or float("nan"),
                 ["%.5f" % v for v in r["cycle_s"]], len(r["beta_extrema"]), r["wall_s"]),
              flush=True)
    with open("data/validation_summary.json", "w") as f:
        json.dump(out, f, indent=1)
    return out

if __name__ == "__main__":
    main()
