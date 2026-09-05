"""F18c: scaled unfolding of the hemicycle. With D ~ pi e0 - pi delta y + c2 y^2 + (second order) da*delta*y*log y,
cycles near the graphic need e0 ~ delta^2 and da ~ delta. Scan those scalings at a0 in {-1, -1/2}, b in {0.5,1,1.5}."""
import sys, json, time, numpy as np
sys.argv = ['x', 'mvneutral', '/dev/null', '0', '0']
import importlib.util; spec = importlib.util.spec_from_file_location('sl', 'sweep_log.py'); sl = importlib.util.module_from_spec(spec); spec.loader.exec_module(sl)
out = open('data/F18c_scaled.jsonl', 'a'); t0 = time.time(); hist = {}
for a0 in (-1.0, -0.5):
    for b in (0.5, 1.0, 1.5):
        for e1 in (1e-2, -1e-2):
            for dfac in (1e-1, 1e-2, 1e-3, -1e-1, -1e-2, -1e-3):
                delta = dfac*e1; e2 = -2*e1 + delta
                for k0 in (0.0, 0.3, 1.0, 3.0, 10.0, -0.3, -1.0, -3.0, -10.0):
                    e0 = k0*delta*delta
                    for kd in (0.0, 0.3, 1.0, 3.0, -0.3, -1.0, -3.0):
                        da = kd*abs(delta); a = a0+da
                        c = np.array([(b-2)/4, e1, 1-b, a, e2, b, e0, 0, 0, 0, -2.0, 0], float)
                        try: nests = sl.evaluate(c)
                        except Exception: continue
                        tot = sum(len(n['roots']) for n in nests); mx = max([len(n['roots']) for n in nests]+[0]); hist[tot] = hist.get(tot, 0)+1
                        if mx >= 3 or tot >= 4:
                            rec = dict(a0=a0, b=b, e0=e0, e1=e1, e2=e2, da=da, coef=c.tolist(), total=tot, nests=nests); out.write(json.dumps(rec)+"\n"); out.flush()
                            print(f"HIT total={tot} a0={a0} b={b} e1={e1:+g} delta={delta:+.2e} e0={e0:+.2e} da={da:+.2e} nests={[(len(n['roots']), ''.join(n['stab']), [f'{r:.3g}' for r in n['roots']]) for n in nests]}", flush=True)
                        if tot >= 5: print("FIVE OR MORE", rec, flush=True)
                print(f"a0={a0} b={b} e1={e1:+g} delta={delta:+.1e} done {time.time()-t0:.0f}s hist {dict(sorted(hist.items()))}", flush=True)
print("DONE", dict(sorted(hist.items())))
