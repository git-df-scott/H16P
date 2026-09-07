#!/usr/bin/env python3
"""Independent ODE cycle counter, resolved, with a positive control.

The first version was under-resolved: 40 samples over d in (0,3] while two of
the engine-certified cycles sit 0.015 apart, so their sign change was invisible.
That is a resolution failure, not a broken detector -- the direction and the
event condition were right (dirtest.py).  This version samples ds = 0.01 and is
validated against the engine on the FIRST nest before being trusted on the second.
"""
import json, sys
import numpy as np
from scipy.integrate import solve_ivp

def field(mu):
    a, a20 = float(mu["a"]), float(mu["a20"])
    a11, a01, a10 = float(mu["a11"]), float(mu["a01"]), float(mu["a10"])
    a00 = a01 + a11 - a10 - a20 - a
    def F(t, z):
        x, y = z
        return [1 + x*y, a00 + a10*x + a20*x*x + a01*y + a11*x*y + a*y*y]
    return F, (a, a20, a11, a01, a10, a00)

def count(mu, fx, fy, sgn, dmax=3.0, ds=0.01, tend=4.0, rtol=1e-10):
    """Cycles crossing the horizontal ray from (fx,fy) in direction sgn."""
    F, _ = field(mu)
    def disp(d):
        z0 = [fx + sgn*d, fy]
        def ev(t, z): return z[1] - fy
        ev.direction = -1
        sol = solve_ivp(F, (0, tend), z0, rtol=rtol, atol=rtol*1e-2, events=ev, max_step=0.05)
        for tt, p in zip(sol.t_events[0], sol.y_events[0]):
            if tt > 1e-6 and sgn*(p[0]-fx) > 0: return sgn*(p[0]-fx) - d
        def ev2(t, z): return z[1] - fy
        ev2.direction = +1
        sol = solve_ivp(F, (0, tend), z0, rtol=rtol, atol=rtol*1e-2, events=ev2, max_step=0.05)
        for tt, p in zip(sol.t_events[0], sol.y_events[0]):
            if tt > 1e-6 and sgn*(p[0]-fx) > 0: return sgn*(p[0]-fx) - d
        return None
    n = int(dmax/ds)
    vals = []
    for k in range(1, n+1):
        d = ds*k; v = disp(d)
        if v is None:
            if vals: break
            continue
        vals.append((d, v))
    roots = []
    for i in range(len(vals)-1):
        if vals[i][1]*vals[i+1][1] < 0:
            lo, hi, f0 = vals[i][0], vals[i+1][0], vals[i][1]
            for _ in range(26):
                mid = (lo+hi)/2; fm = disp(mid)
                if fm is None: break
                if fm*f0 > 0: lo = mid
                else: hi = mid
            roots.append((lo+hi)/2)
    return roots, len(vals), (vals[-1][0] if vals else None)

def second_focus(mu):
    F, (a, a20, a11, a01, a10, a00) = field(mu)
    co = [a20, a10, a00 - a11, -a01, a]
    rs = np.roots(co) if abs(a20) > 1e-14 else np.roots(co[1:])
    for rr in rs:
        if abs(rr.imag) > 1e-9 or abs(rr.real - 1.0) < 1e-6 or abs(rr.real) < 1e-12: continue
        x = rr.real; y = -1.0/x
        J = np.array([[y, x], [a10 + 2*a20*x + a11*y, a01 + a11*x + 2*a*y]])
        det, tr = np.linalg.det(J), np.trace(J)
        if det > 0 and tr*tr - 4*det < 0: return (x, y, tr, det)
    return None

SRC = sys.argv[1] if len(sys.argv) > 1 else "ledger_opus/unfold.json"
KEY = sys.argv[2] if len(sys.argv) > 2 else "mu"
ROWS = json.load(open(SRC))
out = []
for r in ROWS:
    mu = r[KEY] if KEY in r else {k: r[k] for k in ("a","a20","a11","a01","a10")}
    lab = r.get("etaT") or r.get("eps") or "?"
    eng = [float(x) - 1.0 for x in r.get("local", [])]
    roots, nv, reach = count(mu, 1.0, -1.0, +1.0)
    if eng:
        ok = (len(roots) == len(eng)) and all(
            min(abs(v - e) for v in roots) < 2e-3 for e in eng) if roots else False
    else:
        ok = (len(roots) == 0)
    print("\n[%s] FIRST NEST  engine s-1 = %s" % (lab, ["%.6f" % v for v in eng]))
    print("            detector s-1 = %s  (samples %d, reach %.2f)  -> %s"
          % (["%.6f" % v for v in roots], nv, reach or -1, "AGREES" if ok else "DISAGREES"), flush=True)
    sf = second_focus(mu)
    if sf is None:
        print("            SECOND NEST: no second focus"); rem = None
    else:
        fx, fy, tr, det = sf
        rroots, rnv, rreach = count(mu, fx, fy, -1.0)
        rem = len(rroots)
        print("            SECOND NEST focus (%.5f, %.5f) trace %.5f det %.5f" % (fx, fy, tr, det))
        print("            detector cycles at d = %s  (samples %d, reach %.2f)"
              % (["%.6f" % v for v in rroots], rnv, rreach or -1), flush=True)
    out.append(dict(label=lab, first_engine=eng, first_detector=roots, agrees=bool(ok),
                    second_focus=None if sf is None else list(sf), remote=rem))
json.dump(out, open(SRC.replace(".json", "_remote2.json"), "w"), indent=1)
