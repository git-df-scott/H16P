#!/usr/bin/env python3
"""F11: locate the neutral point (saddle trace = 0) on the focus-type loop branch at fixed a, and report eta2 there.
Usage: python3 fable_f11_neutral.py a b_lo b_hi nb l_lo l_hi nl OUT.jsonl
Coarse: for each b, find focus-type loops (|eta2|>1e-6) in [l_lo,l_hi]; bracket sign change of trace in b; brentq."""
import sys, os, json, numpy as np, importlib.util
spec = importlib.util.spec_from_file_location('f11', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fable_f11_neutral_loop.py'))
f11 = importlib.util.module_from_spec(spec); spec.loader.exec_module(f11)
from scipy.optimize import brentq
from multiprocessing import Pool
a, blo, bhi, nb, llo, lhi, nl = [float(v) for v in sys.argv[1:8]]; out = sys.argv[8]
def focus_loops(b):
    ls = np.linspace(llo, lhi, int(nl))
    res = f11.scan_line((a, b, ls))
    return [r for r in res if abs(r['eta2']) > 1e-6]
if __name__ == "__main__":
    bs = np.linspace(blo, bhi, int(nb))
    with Pool(3) as pool:
        coarse = pool.map(focus_loops, bs)
    recs = []
    for b, res in zip(bs, coarse):
        for r in res:
            print(f"coarse a={a} b={b:.4f} l={r['l']:.6f} trace={r['trace']:+.6f} eta2={r['eta2']:+.5e} saddle={np.round(r['saddle'],4).tolist()}", flush=True)
        recs.append((b, res))
    with open(out, 'a') as f:
        for (b0, r0), (b1, r1) in zip(recs[:-1], recs[1:]):
            for p in r0:
                # match the nearest loop in the next b by l
                q = min(r1, key=lambda s: abs(s['l']-p['l']), default=None)
                if q is None or abs(q['l']-p['l']) > 0.3 or p['trace']*q['trace'] >= 0: continue
                lc = [p['l']]
                def tr(b):
                    res = focus_loops(b)
                    if not res: raise RuntimeError("no loop")
                    s = min(res, key=lambda s: abs(s['l']-lc[-1])); lc.append(s['l']); tr.last = s
                    return s['trace']
                try:
                    bstar = brentq(tr, b0, b1, xtol=1e-9)
                except Exception as e:
                    print("brentq failed", e); continue
                s = tr.last
                print(f"NEUTRAL a={a} b*={bstar:.9f} l={s['l']:.9f} trace={s['trace']:+.3e} eta2={s['eta2']:+.6e} saddle={s['saddle']}", flush=True)
                f.write(json.dumps(dict(a=a, bstar=bstar, **s))+"\n")
    print("DONE")
