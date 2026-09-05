"""D2 lemmas (A) and (B): for fixed a, follow the focus-type loop branch in (b,l), bisect for sigma = 0 in b,
report eta_2 there (lemma A: should be a center, eta_2 = 0) and the slopes d sigma/db, d eta_2/db along the branch
(lemma B: opposite signs)."""
import sys, numpy as np, importlib.util
spec = importlib.util.spec_from_file_location('f11', '/home/user/H16P/audit/fable_f11_neutral_loop.py'); f11 = importlib.util.module_from_spec(spec); spec.loader.exec_module(f11)
from multiprocessing import Pool
def loops_at(args):
    a, b = args
    ls = np.linspace(-2.6, -0.4, 45)
    return b, [r for r in f11.scan_line((a, b, ls)) if abs(r['eta2']) > 1e-7]   # focus-type only
def branch_at(a, b, lguess):
    ls = np.linspace(lguess-0.08, lguess+0.08, 17)
    res = [r for r in f11.scan_line((a, b, ls))]
    if not res: return None
    return min(res, key=lambda r: abs(r['l']-lguess))
if __name__ == "__main__":
    a = float(sys.argv[1])
    bs = np.arange(-4.0, 4.01, 0.25)
    with Pool(4) as p: rows = p.map(loops_at, [(a, b) for b in bs])
    pts = sorted([(b, r['l'], r['trace'], r['eta2']) for b, rs in rows for r in rs])
    print(f"a={a}: focus-type loops found at {len(pts)} (b,l) points")
    for b, l, s, e in pts: print(f"   b={b:+.2f} l={l:+.5f} sigma={s:+.5f} eta2={e:+.4e} product={s*e:+.3e}")
    # find sign changes of sigma along the branch (consecutive b with nearby l)
    for (b1, l1, s1, e1), (b2, l2, s2, e2) in zip(pts[:-1], pts[1:]):
        if s1*s2 < 0 and abs(l2-l1) < 0.3:
            lo, hi, llo = b1, b2, l1
            for _ in range(30):
                mid = 0.5*(lo+hi); r = branch_at(a, mid, llo + (mid-lo)/(hi-lo)*(l2-l1) if hi > lo else llo)
                if r is None: break
                if r['trace']*s1 > 0: lo, s1, llo = mid, r['trace'], r['l']
                else: hi = mid
                if hi-lo < 1e-7: break
            r = branch_at(a, 0.5*(lo+hi), llo)
            if r is None: continue
            # slopes along the branch
            rp = branch_at(a, 0.5*(lo+hi)+1e-3, r['l']); rm = branch_at(a, 0.5*(lo+hi)-1e-3, r['l'])
            ds = (rp['trace']-rm['trace'])/2e-3 if rp and rm else float('nan'); de = (rp['eta2']-rm['eta2'])/2e-3 if rp and rm else float('nan')
            l, b = r['l'], 0.5*(lo+hi)
            c1 = a*(b+2*l); c2 = b-3*l-5; c3 = a*a*(b+2*l+1)-(b+1)*(l+1)**2
            print(f"NEUTRAL LOOP a={a} b={b:.7f} l={l:.7f}: sigma={r['trace']:+.2e} eta2={r['eta2']:+.3e} | center factors a(b+2l)={c1:+.3e} (b-3l-5)={c2:+.3e} C3={c3:+.3e} | slopes dsigma/db={ds:+.4f} deta2/db={de:+.4f} product={ds*de:+.3e}")
