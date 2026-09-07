#!/usr/bin/env python3
"""Constrained variation of the swallow-tail indicator across the cusp locus.

  u = (a11, a01, a10),  q = (a, a20, x0),  F = (D, D_x, D_xx),  G = D_xxx.

On F = 0, du/dq = -F_u^{-1} F_q and

    dG/dq = G_q - G_u F_u^{-1} F_q          (linear solves, no explicit inverse)

Same for H = D_xxxx, then for the scale-free indicator nu = G/(H r0), r0 = x0-1:

    dnu/dq = (dG/dq)/(H r0) - G (dH/dq)/(H^2 r0) - [G/(H r0^2)] dr0/dq .

Every engine call is logged (PROTOCOL rule 5)."""
import json, os, sys, time
import mpmath as mp
from engine import Engine
from cusp import Cusp, solve3, wres
mp.mp.dps = 50

LEDGER = "ledger_opus/constrained.jsonl"
os.makedirs("ledger_opus", exist_ok=True)

def solveT(A, b):
    """Solve A^T w = b (3x3) by Gaussian elimination."""
    AT = [[A[j][i] for j in range(3)] for i in range(3)]
    return solve3(AT, b)

def constrained_derivs(eng, c, mu, x0, ha="1e-7", h20="1e-7"):
    r = c.val(mu, x0)
    if r["status"] != "OK": return None
    G, H = r["Dxxx"], r["Dxxxx"]
    Fu = c.jac_mu(mu, x0)                       # 3x3
    if Fu is None: return None
    # G_u, H_u by the same centred differences
    Gu, Hu = [mp.mpf(0)]*3, [mp.mpf(0)]*3
    Gu, Hu = list(Gu), list(Hu)
    from cusp import FD_H
    for j in range(3):
        h = FD_H[j]*max(mp.mpf(1), abs(mu[j]))
        mp_, mm_ = list(mu), list(mu); mp_[j] += h; mm_[j] -= h
        rp, rm = c.val(mp_, x0), c.val(mm_, x0)
        if rp["status"] != "OK" or rm["status"] != "OK": return None
        Gu[j] = (rp["Dxxx"] - rm["Dxxx"])/(2*h)
        Hu[j] = (rp["Dxxxx"] - rm["Dxxxx"])/(2*h)
    # shape columns of F_q, G_q, H_q
    Fq = [[mp.mpf(0)]*3 for _ in range(3)]; Gq = [mp.mpf(0)]*3; Hq = [mp.mpf(0)]*3
    for k, (attr, hstr) in enumerate((("a", ha), ("a20", h20))):
        h = mp.mpf(hstr)*max(mp.mpf(1), abs(getattr(c, attr)))
        old = getattr(c, attr)
        setattr(c, attr, old + h); rp = c.val(mu, x0)
        setattr(c, attr, old - h); rm = c.val(mu, x0)
        setattr(c, attr, old)
        if rp["status"] != "OK" or rm["status"] != "OK": return None
        Fq[0][k] = (rp["D"] - rm["D"])/(2*h)
        Fq[1][k] = (rp["Dx"] - rm["Dx"])/(2*h)
        Fq[2][k] = (rp["Dxx"] - rm["Dxx"])/(2*h)
        Gq[k] = (rp["Dxxx"] - rm["Dxxx"])/(2*h)
        Hq[k] = (rp["Dxxxx"] - rm["Dxxxx"])/(2*h)
    # x0 column is free
    Fq[0][2], Fq[1][2], Fq[2][2] = r["Dx"], r["Dxx"], r["Dxxx"]
    Gq[2] = H
    Hq[2] = mp.mpf(0)                            # D_xxxxx not available; flagged below
    wG = solveT(Fu, Gu); wH = solveT(Fu, Hu)
    if wG is None or wH is None: return None
    dG = [Gq[k] - sum(wG[i]*Fq[i][k] for i in range(3)) for k in range(3)]
    dH = [Hq[k] - sum(wH[i]*Fq[i][k] for i in range(3)) for k in range(3)]
    r0 = mp.mpf(x0) - 1
    nu = G/(H*r0)
    dnu = []
    for k in range(3):
        dr0 = mp.mpf(1) if k == 2 else mp.mpf(0)
        dnu.append(dG[k]/(H*r0) - G*dH[k]/(H*H*r0) - G*dr0/(H*r0*r0))
    return dict(G=G, H=H, nu=nu, dG=dG, dH=dH, dnu=dnu, Fu=Fu, r=r)

if __name__ == "__main__":
    src = sys.argv[1] if len(sys.argv) > 1 else "ledger/cusp_row5.jsonl"
    picks = [int(v) for v in (sys.argv[2].split(",") if len(sys.argv) > 2 else ["0","40","80","119"])]
    recs = [json.loads(l) for l in open(src) if l.strip()]
    eng = Engine()
    print("engine:", eng.banner)
    lg = open(LEDGER, "a")
    print("\nsource: %s" % src)
    print(" idx   x0        nu          dnu/da        dnu/da20      dnu/dx0       cond(F_u)")
    for i in picks:
        if i >= len(recs): continue
        R = recs[i]
        c = Cusp(eng, R["a"], R["a20"], side=int(R.get("side", 1)))
        mu = [mp.mpf(R["a11"]), mp.mpf(R["a01"]), mp.mpf(R["a10"])]
        x0 = mp.mpf(R["x0"])
        t0 = time.time()
        out = constrained_derivs(eng, c, mu, x0)
        if out is None:
            print("  %-5d FAILED" % i); continue
        Fu = out["Fu"]
        nrm = max(sum(abs(Fu[r][k]) for k in range(3)) for r in range(3))
        print("  %-5d %-9.6f %-11.6g %-13.5g %-13.5g %-13.5g %.3g"
              % (i, float(x0), float(out["nu"]), float(out["dnu"][0]),
                 float(out["dnu"][1]), float(out["dnu"][2]), float(nrm)))
        lg.write(json.dumps(dict(src=src, idx=i, a=R["a"], a20=R["a20"], x0=str(x0),
                                 nu=mp.nstr(out["nu"], 20),
                                 dnu_da=mp.nstr(out["dnu"][0], 20),
                                 dnu_da20=mp.nstr(out["dnu"][1], 20),
                                 dnu_dx0=mp.nstr(out["dnu"][2], 20),
                                 engine=eng.banner, calls=eng.ncalls,
                                 wall=time.time()-t0)) + "\n")
        lg.flush()
    print("\ntotal engine calls this run:", eng.ncalls)
    eng.close()
