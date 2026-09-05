"""D2 lemma (B) at a crossing: bisect b in a bracket where sigma changes sign on the focus branch, with a robust
l-window, and report eta2, center factors and slopes at the neutral loop."""
import sys, numpy as np, importlib.util
spec = importlib.util.spec_from_file_location('f11', '/home/user/H16P/audit/fable_f11_neutral_loop.py'); f11 = importlib.util.module_from_spec(spec); spec.loader.exec_module(f11)
a, b1, b2, l1, l2 = [float(v) for v in sys.argv[1:6]]
def branch_at(b, lg, w=0.2):
    res = [r for r in f11.scan_line((a, b, np.linspace(lg-w, lg+w, 41)))]
    return min(res, key=lambda r: abs(r['l']-lg)) if res else None
r1 = branch_at(b1, l1); r2 = branch_at(b2, l2)
print("bracket:", (b1, r1['l'], r1['trace'], r1['eta2']) if r1 else None, (b2, r2['l'], r2['trace'], r2['eta2']) if r2 else None, flush=True)
lo, hi, slo, llo, lhi = b1, b2, r1['trace'], r1['l'], r2['l']
for it in range(28):
    mid = 0.5*(lo+hi); lg = llo + (mid-lo)/(hi-lo)*(lhi-llo)
    r = branch_at(mid, lg)
    if r is None: print("lost branch at b=", mid); break
    print(f"  it {it}: b={mid:.8f} l={r['l']:.7f} sigma={r['trace']:+.3e} eta2={r['eta2']:+.3e}", flush=True)
    if r['trace']*slo > 0: lo, slo, llo = mid, r['trace'], r['l']
    else: hi, lhi = mid, r['l']
    if hi-lo < 2e-7: break
b = 0.5*(lo+hi); r = branch_at(b, 0.5*(llo+lhi)); l = r['l']
rp = branch_at(b+2e-3, l); rm = branch_at(b-2e-3, l)
ds = (rp['trace']-rm['trace'])/4e-3; de = (rp['eta2']-rm['eta2'])/4e-3
c1 = a*(b+2*l); c2 = b-3*l-5; c3 = a*a*(b+2*l+1)-(b+1)*(l+1)**2
print(f"NEUTRAL LOOP a={a} b={b:.8f} l={l:.8f}: sigma={r['trace']:+.2e} eta2={r['eta2']:+.3e} | center factors a(b+2l)={c1:+.3e} (b-3l-5)={c2:+.3e} C3={c3:+.3e} | slopes dsigma/db={ds:+.5f} deta2/db={de:+.5f} product={ds*de:+.3e}")
