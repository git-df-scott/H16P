#!/usr/bin/env python3
"""Connected predictor-corrector continuation ALONG the fold surface.

State: (mu, s_f) with D(s_f)=0, D_s(s_f)=0, D_ss(s_f)!=0, plus a tracked second
stationary point s_i (D_s(s_i)=0) with critical value h_i = D(s_i).

Steering (fold-preserving to first order):
   b_f = d_mu D(s_f);  b_i = d_mu D(s_i)
   w  = b_i - (b_i.b_f)/(b_f.b_f) b_f;   dmu = -sign(h_i)*step*w/|w|
Predictor: ds_f = -(d_mu D_s . dmu)/D_ss(s_f);  ds_i = -(d_mu D_s . dmu)/D_ss(s_i)
Corrector: Newton in (s_f, mu[pidx]) back onto (D, D_s) = 0; s_i re-Newtoned on D_s.

s_i is tracked CONTINUOUSLY by Newton from the previous step, so its loss is an
identified event (D_ss -> 0, domain exit, or return failure), not an absence
from a sampled profile.  A full inventory is taken every INV steps and at every
event.  Adaptive step with real backtracking.
"""
import json, sys
import mpmath as mp
from engine import Engine
from cusp import Cusp
import foldcont as FC
mp.mp.dps = 40

ROW   = int(sys.argv[1]) if len(sys.argv) > 1 else 7
PIDX  = int(sys.argv[2]) if len(sys.argv) > 2 else 2
STEP0 = mp.mpf(sys.argv[3]) if len(sys.argv) > 3 else mp.mpf("5e-4")
NSTEP = int(sys.argv[4]) if len(sys.argv) > 4 else 400
INV   = int(sys.argv[5]) if len(sys.argv) > 5 else 20
TAG   = sys.argv[6] if len(sys.argv) > 6 else ""

def Tval(m): return m["a11"] + m["a01"] - 2*m["a"] - 1
def detJ(m): return 2*m["a"] - m["a01"] - m["a10"] - 2*m["a20"]

def newton_stat(c, mu, s, itmax=40):
    """Newton on D_s(s)=0.  Returns (s, r) or (None, reason)."""
    for _ in range(itmax):
        r = FC.val(c, mu, s)
        if r["status"] != "OK": return None, r["status"]
        if r["Dxx"] == 0: return None, "DSS_ZERO"
        if abs(r["Dx"])*abs(s-1) < mp.mpf("1e-24"):
            return s, r
        d = r["Dx"]/r["Dxx"]
        if abs(d) > mp.mpf("0.5")*abs(s-1): d = mp.sign(d)*mp.mpf("0.5")*abs(s-1)
        s = s - d
        if s <= 1: return None, "LEFT_DOMAIN"
        if abs(d) < mp.mpf("1e-20"):
            r = FC.val(c, mu, s)
            return (s, r) if r["status"] == "OK" else (None, r["status"])
    return None, "STAT_NO_CONVERGE"

rec = [json.loads(l) for l in open("ledger_opus/fold_start.jsonl")]
rec = [r for r in rec if r["row"] == ROW][-1]
mu  = {p: mp.mpf(rec[p]) for p in FC.PARAMS}
s_f = mp.mpf(rec["s_f"])
s_i = mp.mpf(rec["stationary"][0])

eng = Engine(); print("engine:", eng.banner, flush=True)
lg = open("ledger_opus/foldrun_row%d_p%d%s.jsonl" % (ROW, PIDX, TAG), "a")
def emit(d):
    lg.write(json.dumps(d)+"\n"); lg.flush()
emit(dict(event="START", row=ROW, pidx=PIDX, step0=str(STEP0),
          **{p: mp.nstr(mu[p],25) for p in FC.PARAMS},
          s_f=mp.nstr(s_f,20), s_i=mp.nstr(s_i,20)))

step = STEP0
event = None
k = 0
while k < NSTEP:
    c = Cusp(eng, mu["a"], mu["a20"], side=1)
    rf = FC.val(c, mu, s_f)
    si, ri = newton_stat(c, mu, s_i)
    if rf["status"] != "OK":
        event = ("RETURN_FAILED_AT_FOLD", rf["status"]); break
    if si is None:
        event = ("SECOND_STATIONARY_LOST", ri); break
    if abs(si - s_f) < mp.mpf("0.01"):
        event = ("STATIONARY_COLLISION", mp.nstr(abs(si-s_f), 8)); break
    s_i = si; h_i = ri["D"]
    line = dict(step=k, **{p: mp.nstr(mu[p],25) for p in FC.PARAMS},
                s_f=mp.nstr(s_f,20), s_i=mp.nstr(s_i,20), h_i=mp.nstr(h_i,15),
                Dss_f=mp.nstr(rf["Dxx"],12), Dss_i=mp.nstr(ri["Dxx"],12),
                T=mp.nstr(Tval(mu),15), detJ=mp.nstr(detJ(mu),12),
                h=mp.nstr(step,6), calls=eng.ncalls)
    if k % INV == 0:
        inv = FC.inventory(c, mu)
        line["inventory"] = dict(domain=[mp.nstr(v,10) for v in inv["domain"]],
                                 cycles=[mp.nstr(v,12) for v in inv["cycles"]],
                                 stationary=[mp.nstr(v,12) for v in inv["stationary"]])
        print("k=%-4d s_f=%s s_i=%s h_i=%-13s Dss_f=%-11s Dss_i=%-11s T=%-12s cyc=%d sta=%d dom=[%s,%s]"
              % (k, mp.nstr(s_f,8), mp.nstr(s_i,8), mp.nstr(h_i,6), mp.nstr(rf["Dxx"],5),
                 mp.nstr(ri["Dxx"],5), mp.nstr(Tval(mu),6), len(inv["cycles"]),
                 len(inv["stationary"]), mp.nstr(inv["domain"][0],4), mp.nstr(inv["domain"][1],4)),
              flush=True)
    elif k % 5 == 0:
        print("   k=%-4d s_f=%s s_i=%s h_i=%-13s Dss_i=%-11s h=%s"
              % (k, mp.nstr(s_f,8), mp.nstr(s_i,8), mp.nstr(h_i,6),
                 mp.nstr(ri["Dxx"],5), mp.nstr(step,4)), flush=True)
    emit(line)
    if abs(h_i) < mp.mpf("1e-14"):
        event = ("SECOND_FOLD_REACHED", mp.nstr(h_i, 8)); break

    g  = FC.grad_mu(c, mu, s_f, which=("D","Dx"))
    gi = FC.grad_mu(c, mu, s_i, which=("D",))
    if g is None or gi is None:
        event = ("GRADIENT_FAILED", None); break
    bf, bi = g["D"], gi["D"]
    bb = sum(x*x for x in bf)
    pr = sum(x*y for x, y in zip(bi, bf))/bb
    w  = [bi[j] - pr*bf[j] for j in range(5)]
    nw = mp.sqrt(sum(x*x for x in w))
    if nw == 0: event = ("STEERING_DEGENERATE", None); break
    sg = 1 if h_i > 0 else -1
    sc = FC.scales(mu)
    ok = False
    for _ in range(12):                      # real backtracking
        dmu = [-sg*step*w[j]/nw for j in range(5)]
        ds  = -sum(g["Dx"][j]*dmu[j] for j in range(5))/rf["Dxx"]
        mu2 = {p: mu[p] + dmu[j]*sc[j] for j, p in enumerate(FC.PARAMS)}
        c2  = Cusp(eng, mu2["a"], mu2["a20"], side=1)
        m3, s3, rr = FC.correct(c2, mu2, s_f + ds, pidx=PIDX)
        if m3 is not None:
            si2, _ = newton_stat(c2, m3, s_i)
            if si2 is not None and abs(si2 - s3) > mp.mpf("0.01"):
                mu, s_f, s_i, ok = m3, s3, si2, True; break
        step = step/2
        if step < mp.mpf("1e-9"): break
    if not ok:
        event = ("STEP_COLLAPSE", "corrector/tracker failed below step 1e-9"); break
    step = min(step*mp.mpf("1.3"), STEP0*8)
    k += 1

print("\nEVENT:", event, "at step", k, flush=True)
if event:
    c = Cusp(eng, mu["a"], mu["a20"], side=1)
    try:
        inv = FC.inventory(c, mu, step=mp.mpf("0.005"))
        print("   final inventory: domain=%s cycles=%s stationary=%s"
              % ([mp.nstr(v,8) for v in inv["domain"]] if inv["domain"] else None,
                 [mp.nstr(v,12) for v in inv["cycles"]], [mp.nstr(v,12) for v in inv["stationary"]]), flush=True)
        fin = dict(domain=[mp.nstr(v,10) for v in inv["domain"]] if inv["domain"] else None,
                   cycles=[mp.nstr(v,15) for v in inv["cycles"]],
                   stationary=[mp.nstr(v,15) for v in inv["stationary"]])
    except Exception as e:
        fin = {"error": str(e)}
    emit(dict(event=event[0], detail=str(event[1]), step=k,
              **{p: mp.nstr(mu[p],25) for p in FC.PARAMS},
              s_f=mp.nstr(s_f,20), final_inventory=fin, calls=eng.ncalls))
print("calls:", eng.ncalls); eng.close()
