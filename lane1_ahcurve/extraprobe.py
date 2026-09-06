"""Two extra validation probes, written to data/ for VALIDATION.md."""
import json, math
import numpy as np
from fractions import Fraction as F
import engine, seeds, ledger

def ulp(s):
    return 0.5*10.0**(-len(s.split('.')[-1])) if '.' in s else 0.5

def rounding_box():
    out = []
    for sd in seeds.fat_seeds():
        if sd["family"] != "cherkas":
            continue
        p = dict(sd["params"])
        xs = np.array(sd["x_cycles"], float); s = np.abs(xs-1.0); d = sd["direction"]
        L0 = [float(v) for v in sd["local"]]; E = [float(v) for v in sd["evec_local"]]
        t0, _, _, _ = engine.curve_with_noise(L0, d, s, mode="lin", evec=E, p0=0.0, span=2.0)
        shifts = np.zeros_like(t0)
        for key in ("a", "a20", "a01", "a10"):
            if "/" in p[key]:
                continue
            h = ulp(p[key])
            for sg in (+1, -1):
                q = dict(p); q[key] = str(F(p[key]) + F(sg*h))
                v = seeds.cherkas_vec12(q["a"], q["a20"], q["a11"], q["a01"], q["a10"])
                Lp = [float(z) for z in engine.local_expand(v, F(1), F(-1))]
                t, st, _, _ = engine.curve_with_noise(Lp, d, s, mode="lin", evec=E,
                                                     p0=0.0, span=2.0)
                if np.all(np.isfinite(t)):
                    shifts = np.maximum(shifts, np.abs(t-t0))
        out.append(dict(row=sd["row"], residual=[abs(float(v)) for v in t0],
                        shift=[float(v) for v in shifts], a11_ulp=ulp(p["a11"])))
    json.dump(out, open("data/rounding_box.json", "w"), indent=1)
    ledger.append("validation", dict(kind="rounding_box", engine=engine.ENGINE, rows=out))
    return out

def remote_probe():
    out = []
    for nm in ("cherkas7", "cherkas8"):
        sd = [x for x in seeds.fat_seeds() if x["name"] == nm][0]
        G = [float(v) for v in sd["vec12"]]
        eq = [e for e in engine.equilibria(G) if abs(e[0]-1) > 1e-6][0]
        L = engine.local_expand(G, eq[0], eq[1])
        Eloc = [float(v) for v in engine.local_expand(
            [F(0)]*6 + [F(1), F(0), F(0), F(0), F(1), F(0)],
            F(eq[0]).limit_denominator(10**12), F(eq[1]).limit_denominator(10**12))]
        a11 = float(sd["params"]["a11"])
        lo, hi, rays = math.inf, -math.inf, []
        for k in range(12):
            th = 2*math.pi*k/12
            d = (math.cos(th), math.sin(th))
            s = np.geomspace(1e-3, 50, 80)
            t, st, nz, _ = engine.curve_with_noise(L, d, s, mode="lin", evec=Eloc,
                                                   p0=0.0, span=40.0)
            good = np.isfinite(t) & (st == 0)
            D, dst, dnz = engine.displacement_with_noise(L, d, s, b=0.0)
            br = engine.count_sign_changes(s, D, dnz)
            if good.any():
                lo = min(lo, a11 + float(t[good].min()))
                hi = max(hi, a11 + float(t[good].max()))
            rays.append(dict(theta=th, n_ok=int(good.sum()), n_sign_changes=len(br),
                             s_max=float(s[good][-1]) if good.any() else None))
        out.append(dict(row=sd["row"], seed=nm, focus=[float(eq[0]), float(eq[1])],
                        trace=float(L[1]+L[8]), a11=a11, lo=lo, hi=hi, rays=rays,
                        status="UNRESOLVED"))
    json.dump(out, open("data/remote_probe.json", "w"), indent=1)
    ledger.append("validation", dict(kind="remote_probe", engine=engine.ENGINE, rows=out))
    return out

if __name__ == "__main__":
    rounding_box(); remote_probe(); print("ok")
