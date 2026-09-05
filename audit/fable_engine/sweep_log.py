"""F15: compactified cycle counting (log-polar returns, radius up to e^umax) over unfoldings of neutral-hemicycle
two-center systems.  Families:
  mvneutral : x' = (b-2)/4 + (1-b) y + a x^2 + b y^2 + e1 x + e2 xy ,  y' = -2xy + e0 ;  a = -1 + da
  kklstar   : KKL x'=y+x^2+xy, y'=-10x^2+(11/5)xy+c y^2+al x+be y near c* (J=0), free (c, K, be) and small changes of -10, 11/5
Usage: python3 sweep_log.py FAMILY OUT.jsonl NSETS SEED [--umax=40] [--store=3]
"""
import sys, json, time, numpy as np
import retmap as rm
from scipy.stats import qmc
fam, out, N, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
UMAX = 36.0; STORE = 3
for a in sys.argv:
    if a.startswith('--umax='): UMAX = float(a[7:])
    if a.startswith('--store='): STORE = int(a[8:])
rng = np.random.default_rng(seed)
def logu(u):   # log-uniform magnitude with random sign from u in [0,1)
    return (10**(-6+5*u[0]))*(1 if u[1] < 0.5 else -1)
def make(u):
    if fam == 'mvneutral':
        da = logu(np.array([u[0]*0.8+0.2, u[1]]))*10; b = 0.2+1.6*u[1]; a = -1+da   # |da| log-uniform in [1e-4, 1e-1]
        e0 = logu(u[2:4]); e1 = logu(u[4:6]); e2 = logu(u[6:8])
        c = np.array([(b-2)/4, e1, 1-b, a, e2, b, e0, 0, 0, 0, -2.0, 0], float)
        return c, dict(a=a, b=b, e0=e0, e1=e1, e2=e2)
    if fam == 'kklstar':
        cs = 0.968620633553494; c = cs + (2*u[0]-1)*0.02; K = 10**(-4+3*u[1]) * (1 if u[2] < 0.7 else -1)
        al = -(K+42)/(11*c/5-1); be = logu(u[3:5])*0.1; p = -10+(2*u[5]-1)*0.5; q = 2.2+(2*u[6]-1)*0.2
        cf = np.array([0, 0, 1, 1, 1, 0, 0, al, be, p, q, c], float)
        return cf, dict(c=c, K=K, al=al, be=be, p=p, q=q)
    if fam == 'shi':
        lam = (2*u[0]-1)*0.2; l, m, n, a, b = (2*u[1:6]-1)*4.0
        return rm.shi_coef(lam, l, m, a, b, n), dict(lam=lam, l=l, m=m, n=n, a=a, b=b)
    raise SystemExit
DIM = {'mvneutral': 8, 'kklstar': 7, 'shi': 6}[fam]
sob = qmc.Sobol(d=DIM, scramble=True, seed=seed)
def count_nest_log(c, pt, other, umax=UMAX):
    d = rm.away_dir(pt, other); th = float(np.arctan2(d[1], d[0]))
    scale = min([np.hypot(pt[0]-q[0], pt[1]-q[1]) for q in other if np.hypot(pt[0]-q[0], pt[1]-q[1]) > 1e-9] + [1.0])
    RT = dict(th0=th, umax=umax+4, Smax=2e4, maxsteps=600_000)
    NOISE = 5e-12
    def ret(us):
        """two-tolerance return: value at rtol 1e-12 and a per-point noise estimate from the 1e-11 run"""
        us = np.asarray(us, float)[None]
        u1, S, st = rm.returns_log(c[None], np.array([pt]), us, rtol=1e-12, **RT)
        u2, S2, st2 = rm.returns_log(c[None], np.array([pt]), us, rtol=1e-11, **RT)
        noise = np.where((st[0] == 0) & (st2[0] == 0), 10*np.abs(u1[0]-u2[0]) + NOISE, np.inf)
        return u1[0], st[0], noise
    u0 = np.arange(np.log(1e-3*scale), umax, 0.25)
    u1, st, nz = ret(u0)
    ok = st == 0; k = len(u0) if ok.all() else int(np.argmin(ok))
    uu = u0[:k]; D = (u1[:k]-uu); NZ = nz[:k]
    if 1 <= k < len(u0):   # edge bisection in u
        lo, hi = u0[k-1], u0[k]; De = None; Ne = np.inf
        for _ in range(8):
            mid = 0.5*(lo+hi); v1, s1, n1 = ret([mid])
            if s1[0] == 0: lo = mid; De = v1[0]-mid; Ne = n1[0]
            else: hi = mid
        if De is not None and lo > u0[k-1]: uu = np.append(uu, lo); D = np.append(D, De); NZ = np.append(NZ, Ne)
    # adaptive refinement: interior local minima of |D| without a sign change (near-fold pairs) and
    # any interval whose endpoints are both below 20x the noise floor (small displacements)
    for _ in range(2):
        add = []
        for i in range(1, len(D)-1):
            same = D[i-1]*D[i] > 0 and D[i]*D[i+1] > 0
            if same and abs(D[i]) < abs(D[i-1]) and abs(D[i]) < abs(D[i+1]) and abs(D[i]) < 0.2*max(abs(D[i-1]), abs(D[i+1])):
                add.extend(np.linspace(uu[i-1], uu[i+1], 9)[1:-1])
        add = [a for a in add if not np.any(np.abs(uu-a) < 1e-9)]
        if not add: break
        v1, s1, n1 = ret(add); good = s1 == 0
        uu = np.concatenate([uu, np.asarray(add)[good]]); D = np.concatenate([D, (v1-np.asarray(add))[good]]); NZ = np.concatenate([NZ, n1[good]])
        o = np.argsort(uu); uu, D, NZ = uu[o], D[o], NZ[o]
    idx = [i for i in range(len(D)-1) if D[i]*D[i+1] < 0 and min(abs(D[i]), abs(D[i+1])) > max(NZ[i], NZ[i+1])]
    roots = [float(np.exp(0.5*(uu[i]+uu[i+1]))) for i in idx]; stab = ['S' if D[i] > 0 else 'U' for i in idx]
    return roots, stab, int(k), float(np.exp(uu[-1])) if len(uu) else None
def evaluate(c):
    eq = rm.equilibria(c); res = []
    for a in rm.antisaddles(c):
        roots, stab, k, redge = count_nest_log(c, a['pt'], [tuple(p) for p in eq])
        res.append(dict(pt=[float(a['pt'][0]), float(a['pt'][1])], tr=float(a['tr']), roots=roots, stab=stab, k=k, redge=redge))
    return res
if __name__ == "__main__":
    t0 = time.time(); hist = {}; done = 0
    with open(out, 'a') as f:
        while done < N:
            for u in sob.random(64):
                c, meta = make(u)
                try: nests = evaluate(c)
                except Exception as ex: continue
                tot = sum(len(n['roots']) for n in nests); mx = max([len(n['roots']) for n in nests] + [0])
                hist[tot] = hist.get(tot, 0)+1
                if mx >= STORE or tot >= STORE+1:
                    rec = dict(meta); rec.update(coef=[float(v) for v in c], total=tot, nests=nests); f.write(json.dumps(rec)+"\n"); f.flush()
                    if tot >= 5: print("FIVE OR MORE:", rec, flush=True)
                done += 1
            print(f"{done} sets {time.time()-t0:.0f}s hist {dict(sorted(hist.items()))}", flush=True)
