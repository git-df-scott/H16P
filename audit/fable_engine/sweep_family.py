"""Generic nest-count sweep over a parametrized family of quadratic fields.
Usage: python3 sweep_family.py FAMILY OUT.jsonl NSETS SEED [--eps=1e-3] [--L=..]
Families:
  q4pert : Q4 center (rho=1,kappa=2) plus eps * random unit direction in the 12 coefficients
  mv     : Marin-Villadelprat family, a=-1: x'=(b-2)/4+e1 x+(1-b)y - x^2+e2 xy+b y^2, y'=e0-2xy
  kklx   : KKL with all quadratic coefficients free: x'=y+e x^2+g xy, y'=p x^2+q xy+c y^2+al x+be y
"""
import sys, json, time, numpy as np
import retmap as rm
from scipy.stats import qmc

fam, out, N, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
eps = 1e-3; L = 4.0
for a in sys.argv:
    if a.startswith('--eps='): eps = float(a[6:])
    if a.startswith('--L='): L = float(a[4:])

Q4 = np.array([0, 0, 1, 6, 4, -2, 0, -1, 0, 2, 8, -2], float)  # rho=1: b=0,c=2

def make(u):
    if fam == 'q4pert':
        v = rng.standard_normal(12); v /= np.linalg.norm(v)
        return Q4 + eps*v, dict(dir=[float(t) for t in v], eps=eps)
    if fam == 'mv':
        b = 0.1 + 3.9*u[0]; e0, e1, e2 = (2*u[1:4]-1)*L
        c = np.array([(b-2)/4, e1, 1-b, -1.0, e2, b, e0, 0, 0, 0, -2.0, 0], float)
        return c, dict(b=b, e0=e0, e1=e1, e2=e2)
    if fam == 'kklx':
        e, g = 1+(2*u[0]-1)*L, 1+(2*u[1]-1)*L
        p = -10+(2*u[2]-1)*4*L; q = 2.2+(2*u[3]-1)*L; c = 0.7+(2*u[4]-1)*0.5
        al = -73+(2*u[5]-1)*30; be = (2*u[6]-1)*0.05
        return np.array([0, 0, 1, e, g, 0, 0, al, be, p, q, c], float), dict(e=e, g=g, p=p, q=q, c=c, al=al, be=be)
    if fam == 'q3rpert':
        # reversible center at the Shi loop point (a=1), CLAUDE_ROUTES_4AB: X'=-Y(1+kX), Y'=X+pX^2+qY^2
        k, p, q = 5.54048179, -1.24519487, 0.22849752
        base = np.array([0, 0, -1, 0, -k, 0, 0, 1, 0, p, 0, q], float)
        v = rng.standard_normal(12); v /= np.linalg.norm(v)
        return base + eps*v, dict(dir=[float(t) for t in v], eps=eps)
    if fam == 'mvpert':
        # Marin-Villadelprat two-center Q3R at a0=-1, b0 in (0,2), plus eps * random direction
        b = 0.1 + 1.8*u[0]
        base = np.array([(b-2)/4, 0, 1-b, -1.0, 0, b, 0, 0, 0, 0, -2.0, 0], float)
        v = rng.standard_normal(12); v /= np.linalg.norm(v)
        return base + eps*v, dict(b=b, dir=[float(t) for t in v], eps=eps)
    raise SystemExit("unknown family")

dim = {'q4pert': 1, 'mv': 4, 'kklx': 7, 'q3rpert': 1, 'mvpert': 1}[fam]
rng = np.random.default_rng(seed)
sob = qmc.Sobol(d=dim, scramble=True, seed=seed)
NR = 64; RMIN, RMAX = 1e-4, 1e4
B = 256; t0 = time.time(); done = 0; hist = {}
with open(out, 'a') as f:
    while done < N:
        U = sob.random(B); sets = []
        for u in U:
            c, meta = make(u)
            try:
                ants = rm.antisaddles(c); eq = rm.equilibria(c)
            except Exception:
                continue
            for k, an in enumerate(ants):
                x, y = an['pt']
                others = [np.hypot(x-p[0], y-p[1]) for p in eq if np.hypot(x-p[0], y-p[1]) > 1e-9]
                scale = min(others) if others else 1.0
                rad = np.geomspace(RMIN*scale, RMAX*max(scale, 1.0), NR)
                m = dict(meta); m.update(nest=k, tr=float(an['tr']), pt=[float(x), float(y)])
                sets.append((c, (x, y), rm.away_dir((x, y), eq), rad, m))
        if not sets:
            done += B; continue
        coef = np.array([s[0] for s in sets]); foc = np.array([s[1] for s in sets]); dr = np.array([s[2] for s in sets]); rad = np.array([s[3] for s in sets])
        R, T, st = rm.returns(coef, foc, dr, rad, 1e-12, 1e7, 5e3, 500000)
        ks = [NR if (st[i] == 0).all() else int(np.argmin(st[i] == 0)) for i in range(len(sets))]
        edge = [i for i, k in enumerate(ks) if 1 <= k < NR]
        Dedge = {}
        if edge:
            re_, De_ = rm.edge_refine(coef[edge], foc[edge], dr[edge], [rad[i, ks[i]-1] for i in edge], [rad[i, ks[i]] for i in edge], 1e-12, 1e7, 5e3, 500000)
            for i, r_, d_ in zip(edge, re_, De_):
                if np.isfinite(d_): Dedge[i] = (r_, d_)
        for i, s in enumerate(sets):
            k = ks[i]
            if k < 2:
                cnt = 0; roots = []
            else:
                rr = rad[i, :k]; D = R[i, :k]-rr
                if i in Dedge:
                    rr = np.append(rr, Dedge[i][0]); D = np.append(D, Dedge[i][1])
                idx = rm.count_signs(rr, D); cnt = len(idx); roots = [float(rr[j]) for j in idx]
            hist[cnt] = hist.get(cnt, 0)+1
            if cnt >= 2:
                rec = dict(s[4]); rec.update(coef=[float(v) for v in s[0]], count=cnt, roots=roots, kvalid=k,
                                             D=[float(v) for v in (R[i, :k]-rad[i, :k])], rad=[float(v) for v in rad[i, :k]])
                f.write(json.dumps(rec)+"\n"); f.flush()
        done += B
        if done % (B*8) == 0:
            print(f"{done} sets, {time.time()-t0:.0f}s, hist {dict(sorted(hist.items()))}", flush=True)
print("DONE", done, time.time()-t0, dict(sorted(hist.items())))
