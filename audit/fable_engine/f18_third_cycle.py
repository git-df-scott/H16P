"""F18: unfold the neutral hemicycle near its codim-2 point e0 = 0, 2e1 + e2 = 0 with the ratio parameter da
and small residuals; count cycles per nest with the compactified counter.  Looks for three or more cycles born
near the upper hemicycle (upper nest count >= 3) or total >= 5."""
import sys, json, time, numpy as np
sys.argv = ['x', 'mvneutral', '/dev/null', '0', '0']
import importlib.util; spec = importlib.util.spec_from_file_location('sl', 'sweep_log.py'); sl = importlib.util.module_from_spec(spec); spec.loader.exec_module(sl)
out = open('data/F18_third_cycle.jsonl', 'a'); t0 = time.time(); hist = {}; best = 0
rng = np.random.default_rng(11)
for b in (0.5, 1.0, 1.5):
    for e1 in (1e-2, 1e-3, 1e-4, -1e-2, -1e-3, -1e-4):
        for da in (1e-4, 1e-3, 1e-2, 3e-2, -1e-4, -1e-3, -1e-2, -3e-2):
            for res in (0.0, 0.1, -0.1, 0.5, -0.5):           # residual of 2e1+e2 relative to e1
                for e0f in (0.0, 1e-3, -1e-3, 1e-1, -1e-1):  # e0 relative to |e1|
                    e2 = -2*e1 + res*e1; e0 = e0f*abs(e1); a = -1+da
                    c = np.array([(b-2)/4, e1, 1-b, a, e2, b, e0, 0, 0, 0, -2.0, 0], float)
                    try: nests = sl.evaluate(c)
                    except Exception as ex: continue
                    tot = sum(len(n['roots']) for n in nests); mx = max([len(n['roots']) for n in nests]+[0])
                    hist[tot] = hist.get(tot, 0)+1
                    if mx >= 3 or tot >= 4:
                        rec = dict(b=b, da=da, e0=e0, e1=e1, e2=e2, coef=c.tolist(), total=tot, nests=nests); out.write(json.dumps(rec)+"\n"); out.flush()
                        print(f"HIT total={tot} b={b} da={da:+g} e0={e0:+g} e1={e1:+g} e2={e2:+g} nests={[(len(n['roots']), ''.join(n['stab']), [f'{r:.3g}' for r in n['roots']]) for n in nests]}", flush=True)
                    if tot >= 5: print("FIVE OR MORE", rec, flush=True)
            print(f"b={b} e1={e1:+g} da={da:+g} done {time.time()-t0:.0f}s hist {dict(sorted(hist.items()))}", flush=True)
print("DONE", dict(sorted(hist.items())))
