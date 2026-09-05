#!/usr/bin/env python3
"""F11 branch tracking: at fixed a, follow the focus-type loop branch in b and record saddle trace vs eta2."""
import sys, os, json, numpy as np, importlib.util
spec = importlib.util.spec_from_file_location('f11', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fable_f11_neutral_loop.py'))
f11 = importlib.util.module_from_spec(spec); spec.loader.exec_module(f11)
from multiprocessing import Pool
a = float(sys.argv[1]); bs = np.linspace(float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])); out = sys.argv[5]
def job(b):
    ls = np.linspace(-3.0, 0.5, 36)
    return f11.scan_line((a, b, ls))
if __name__ == "__main__":
    with Pool(3) as pool, open(out, 'a') as f:
        for res in pool.imap(job, bs):
            for r in res:
                f.write(json.dumps(r)+"\n"); f.flush()
                tag = "CENTER" if abs(r['eta2']) < 1e-6 else "FOCUS"
                print(f"a={r['a']:.2f} b={r['b']:.3f} l={r['l']:.6f} which={r['which']} saddle={np.round(r['saddle'],4).tolist()} trace={r['trace']:+.5f} eta2={r['eta2']:+.4e} {tag}", flush=True)
