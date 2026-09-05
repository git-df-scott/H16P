"""F15 targeted: structured unfolding of the neutral two-center system (MV family, a = -1 + da).
For each (b, da) and a set of perturbation directions (e0,e1,e2) on a sphere with log magnitudes,
count cycles per nest with the compactified counter.  Records everything (small runs)."""
import sys, json, time, numpy as np
sys.argv = ['x', 'mvneutral', '/dev/null', '0', '0'] + sys.argv[1:]
import importlib.util; spec = importlib.util.spec_from_file_location('sl', 'sweep_log.py'); sl = importlib.util.module_from_spec(spec); spec.loader.exec_module(sl)
out = open(sys.argv[5] if len(sys.argv) > 5 else 'data/F15_targeted.jsonl', 'a')
rng = np.random.default_rng(7)
t0 = time.time(); hist = {}
for b in (0.5, 1.0, 1.5):
    for da in (-0.02, -0.005, 0.0, 0.005, 0.02):
        for mag in (1e-4, 1e-3, 1e-2):
            for _ in range(6):
                v = rng.standard_normal(3); v /= np.linalg.norm(v); e0, e1, e2 = mag*v
                a = -1+da
                c = np.array([(b-2)/4, e1, 1-b, a, e2, b, e0, 0, 0, 0, -2.0, 0], float)
                try: nests = sl.evaluate(c)
                except Exception as ex: print("err", ex); continue
                tot = sum(len(n['roots']) for n in nests); hist[tot] = hist.get(tot, 0)+1
                rec = dict(b=b, da=da, mag=mag, e=[e0, e1, e2], coef=c.tolist(), total=tot, nests=nests)
                out.write(json.dumps(rec)+"\n"); out.flush()
                print(f"b={b} da={da:+.3f} mag={mag:g} total={tot} per nest={[ (len(n['roots']), ''.join(n['stab'])) for n in nests]} {time.time()-t0:.0f}s", flush=True)
                if tot >= 5: print("FIVE OR MORE", rec, flush=True)
print("DONE", dict(sorted(hist.items())))
