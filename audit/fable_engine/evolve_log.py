"""Evolutionary descent with the compactified counter (radius to e^40). Seeds: jsonl records with 'coef' and 'total'.
Score = total + 0.5*(largest nest count). Usage: python3 evolve_log.py SEEDS.jsonl OUT.jsonl GEN POP SIGMA"""
import sys, json, time, numpy as np
seeds_f, out_f, GEN, POP, SIG = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
sys.argv = ['x', 'mvneutral', '/dev/null', '0', '0']
import importlib.util; spec = importlib.util.spec_from_file_location('sl', 'sweep_log.py'); sl = importlib.util.module_from_spec(spec); spec.loader.exec_module(sl)
recs = [json.loads(l) for l in open(seeds_f)]; recs = [r for r in recs if r.get('total', 0) >= 4]
rng = np.random.default_rng(int(time.time()) % 99991)
def score(nests): return sum(len(n['roots']) for n in nests) + 0.5*max([len(n['roots']) for n in nests]+[0])
elites = []
for r in recs:
    c = np.array(r['coef']); nests = sl.evaluate(c); elites.append((score(nests), c, nests))
elites.sort(key=lambda t: -t[0]); print("seeds", len(elites), "best", elites[0][0], flush=True)
DIMS = [i for i in range(12) if abs(elites[0][1][i]) > 0 or i in (1, 8)]
best_total = max(sum(len(n['roots']) for n in e[2]) for e in elites); sig = SIG
with open(out_f, 'a') as f:
    for g in range(GEN):
        t0 = time.time(); res = []
        for _ in range(POP):
            par = elites[rng.integers(min(len(elites), 8))][1].copy(); scale = np.where(np.abs(par) > 1e-6, np.abs(par), 1.0)
            step = np.zeros(12); step[DIMS] = rng.standard_normal(len(DIMS)); c = par + sig*scale*step
            try: nests = sl.evaluate(c)
            except Exception: continue
            tot = sum(len(n['roots']) for n in nests); res.append((score(nests), c, nests))
            if tot >= best_total:
                f.write(json.dumps(dict(gen=g, total=tot, coef=c.tolist(), nests=nests))+"\n"); f.flush()
            if tot > best_total:
                best_total = tot; print("NEW MAX TOTAL", tot, c.tolist(), [(''.join(n['stab']), [f'{x:.3g}' for x in n['roots']]) for n in nests], flush=True)
            if tot >= 5: print("FIVE OR MORE", c.tolist(), flush=True)
        elites = sorted(elites+res, key=lambda t: -t[0])[:12]
        tots = [sum(len(n['roots']) for n in r[2]) for r in res]
        print(f"gen {g}: sigma {sig:.3g} best {elites[0][0]:.2f} totals {dict(zip(*np.unique(tots, return_counts=True)))} {time.time()-t0:.0f}s", flush=True)
        if g % 2: sig *= 0.97
print("DONE best_total", best_total)
