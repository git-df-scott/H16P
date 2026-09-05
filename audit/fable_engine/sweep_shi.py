"""F3/F5 sweep: Shi chart x' = lam x - y + l x^2 + m xy + n y^2, y' = x + a x^2 + b xy.
Counts limit cycles in every nest by sign changes of the return displacement.
Usage: python3 sweep_shi.py OUT.jsonl NSETS SEED [--lam0] [--L=4]
"""
import sys, json, time, numpy as np
import retmap as rm
from scipy.stats import qmc

out = sys.argv[1]; N = int(sys.argv[2]); seed = int(sys.argv[3])
lam0 = '--lam0' in sys.argv
store1 = '--store1' in sys.argv
L = 4.0
for a in sys.argv:
    if a.startswith('--L='): L = float(a[4:])
NR = 40; RMIN, RMAX = 1e-3, 3e3
dim = 5 if lam0 else 6
sob = qmc.Sobol(d=dim, scramble=True, seed=seed)
B = 256
t0 = time.time(); done = 0; hist = {}
with open(out, 'a') as f:
    while done < N:
        u = sob.random(B)
        pars = (2*u-1)*L
        sets = []   # (coef, focus, dir, radii, meta)
        for row in pars:
            if lam0:
                l, m, n, a, b = row; lam = 0.0
            else:
                lam, l, m, n, a, b = row; lam *= 0.1
            c = rm.shi_coef(lam, l, m, a, b, n)
            try:
                ants = rm.antisaddles(c)
            except Exception:
                continue
            eq = rm.equilibria(c)
            for k, an in enumerate(ants):
                x, y = an['pt']
                others = [np.hypot(x-p[0], y-p[1]) for p in eq if np.hypot(x-p[0], y-p[1]) > 1e-9]
                scale = min(others) if others else 1.0
                rad = np.geomspace(RMIN*scale, RMAX*max(scale, 1.0), NR)
                sets.append((c, (x, y), (1.0, 0.0), rad, dict(par=[float(v) for v in row], lam=lam, nest=k, tr=an['tr'], pt=[x, y])))
        if not sets:
            continue
        coef = np.array([s[0] for s in sets]); foc = np.array([s[1] for s in sets]); dr = np.array([s[2] for s in sets]); rad = np.array([s[3] for s in sets])
        R, T, st = rm.returns(coef, foc, dr, rad, 1e-9, 1e6, 2e3, 300000)
        for i, s in enumerate(sets):
            ok = st[i] == 0
            k = NR if ok.all() else int(np.argmin(ok))
            if k < 2:
                cnt = 0; roots = []
            else:
                D = R[i, :k]-rad[i, :k]
                sc = np.sign(D)
                idx = np.nonzero(sc[:-1]*sc[1:] < 0)[0]
                cnt = len(idx); roots = [float(rad[i, j]) for j in idx]
            hist[cnt] = hist.get(cnt, 0)+1
            if cnt >= 2 or (store1 and cnt == 1):
                rec = dict(s[4]); rec.update(coef=[float(v) for v in s[0]], count=cnt, roots=roots, kvalid=k, D=[float(v) for v in (R[i, :k]-rad[i, :k])], rad=[float(v) for v in rad[i, :k]])
                f.write(json.dumps(rec)+"\n"); f.flush()
        done += B
        if done % (B*8) == 0:
            print(f"{done} sets, {time.time()-t0:.0f}s, hist {dict(sorted(hist.items()))}", flush=True)
print("DONE", done, time.time()-t0, dict(sorted(hist.items())))
