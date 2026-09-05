"""F4/F8: evolutionary exploration around a multi-cycle seed.
Score = total cycles (all nests) + (1 - smallest normalized gap of a non-crossing extremum of D/r).
Usage: python3 evolve.py SEED_NAME OUT.jsonl GENERATIONS POP SIGMA0 [--dims=i,j,k]
"""
import sys, json, time, numpy as np
import retmap as rm

SEEDS = {
 'kkl': np.array([0, 0, 1, 1, 1, 0, 0, -363889/5000, 3/2000, -10, 11/5, 0.7], float),
 'yz':  np.array([0, 1/20000, 1, 0, -30/7, 0, 0, -1, -500001/1e10, 1, 49182857/96810000000, -671/210], float),
 'shi': np.array([0, -2e-8, -1, -10, 5-0.1, 1, 0, 1, 0, 1, -25+8*(-1e-3)-9*(-0.1), 0], float),
}
name, out, GEN, POP, SIG = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
dims = list(range(12))
for a in sys.argv:
    if a.startswith('--dims='): dims = [int(t) for t in a[7:].split(',')]
if name.startswith('file:'):
    recs = [json.loads(l) for l in open(name[5:])]
    recs = [r for r in recs if 'coef' in r]
    recs.sort(key=lambda r: -r['count'])
    FILE_SEEDS = [np.array(r['coef']) for r in recs[:24]]
    c0 = FILE_SEEDS[0]
else:
    FILE_SEEDS = None
    c0 = SEEDS[name]
NR = 48
rng = np.random.default_rng(int(time.time()) % 100000)

def evaluate(cs):
    """cs: list of coefficient vectors. Returns list of dicts with total count, nests, score."""
    sets = []
    for si, c in enumerate(cs):
        try:
            ants = rm.antisaddles(c); eq = rm.equilibria(c)
        except Exception:
            continue
        for k, an in enumerate(ants):
            x, y = an['pt']
            others = [np.hypot(x-p[0], y-p[1]) for p in eq if np.hypot(x-p[0], y-p[1]) > 1e-9]
            scale = min(others) if others else 1.0
            rad = np.geomspace(1e-4*scale, 1e4*scale, NR)
            sets.append((si, k, c, (x, y), rad, an['tr']))
    res = [dict(total=0, nests=[], score=-1.0, coef=[float(v) for v in c]) for c in cs]
    if not sets:
        return res
    coef = np.array([s[2] for s in sets]); foc = np.array([s[3] for s in sets])
    dr = np.tile([[1.0, 0.0]], (len(sets), 1)); rad = np.array([s[4] for s in sets])
    R, T, st = rm.returns(coef, foc, dr, rad, 1e-10, 1e7, 2e3, 300000)
    for i, s in enumerate(sets):
        ok = st[i] == 0; k = NR if ok.all() else int(np.argmin(ok))
        if k < 3:
            continue
        D = R[i, :k]-rad[i, :k]; q = D/rad[i, :k]; sc = np.sign(D)
        idx = np.nonzero(sc[:-1]*sc[1:] < 0)[0]
        roots = [float(np.sqrt(rad[i, j]*rad[i, j+1])) for j in idx]
        stab = ['S' if D[j] > 0 else 'U' for j in idx]
        qmax = np.max(np.abs(q)) + 1e-300
        gaps = []
        for j in range(1, k-1):
            if (q[j]-q[j-1])*(q[j+1]-q[j]) < 0 and sc[j-1] == sc[j] == sc[j+1] and rad[i, j] > 3e-4*rad[i, -1]/1e8*1e8:
                gaps.append((float(rad[i, j]), float(q[j]/qmax)))
        r = res[s[0]]
        r['total'] += len(roots)
        r['nests'].append(dict(pt=[float(s[3][0]), float(s[3][1])], tr=float(s[5]), roots=roots, stab=stab, kvalid=k, gaps=gaps))
    for r in res:
        gmin = min([abs(g[1]) for n in r['nests'] for g in n['gaps']] + [1.0])
        r['score'] = r['total'] + (1.0-gmin)
    return res

base = evaluate([c0])[0]
print("seed:", base['total'], [(n['roots'], n['stab']) for n in base['nests']], "score", round(base['score'], 3), flush=True)
elites = [(base['score'], c0.copy())]
if FILE_SEEDS is not None:
    for r0 in evaluate(FILE_SEEDS[1:]):
        elites.append((r0['score'], np.array(r0['coef'])))
    elites.sort(key=lambda t: -t[0])
    print('file seeds loaded:', len(elites), 'best', elites[0][0], flush=True)
scale = np.where(np.abs(c0) > 1e-6, np.abs(c0), 1.0)
sig = SIG; best_total = base['total']
with open(out, 'a') as f:
    f.write(json.dumps(dict(gen=-1, **base))+"\n")
    for g in range(GEN):
        t0 = time.time(); cs = []
        for _ in range(POP):
            parent = elites[rng.integers(len(elites))][1].copy()
            step = np.zeros(12); step[dims] = rng.standard_normal(len(dims))
            cs.append(parent + sig*scale*step)
        res = evaluate(cs)
        pool = elites + [(r['score'], np.array(r['coef'])) for r in res]
        pool.sort(key=lambda t: -t[0]); elites = pool[:max(4, POP//8)]
        for r in res:
            if r['total'] >= best_total or r['total'] >= 5:
                f.write(json.dumps(dict(gen=g, **r))+"\n"); f.flush()
        tot = [r['total'] for r in res]
        mx = max(tot) if tot else 0
        if mx > best_total:
            best_total = mx; print("  NEW MAX TOTAL", mx, flush=True)
        print(f"gen {g}: sigma {sig:.3g} best score {elites[0][0]:.4f} totals {dict(zip(*np.unique(tot, return_counts=True)))} {time.time()-t0:.0f}s", flush=True)
        if mx >= 5:
            print("FIVE OR MORE:", [r for r in res if r['total'] >= 5][:3], flush=True)
        sig *= 0.97 if g % 2 else 1.0
print("DONE best_total", best_total, "elite score", elites[0][0])
