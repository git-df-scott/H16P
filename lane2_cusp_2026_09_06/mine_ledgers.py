#!/usr/bin/env python3
"""Mine every saved Lane 2 cusp ledger for a sign change of G = D_xxx along the
cusp curve.  A sign change on a regular cusp path is a swallow-tail (multiplicity
four cycle).  Also record which side of the centre curve a20_c(a) each shape is on.

a20_c(a) = 4a(a+1)(a-2)^2 / [ (a-1)(2a+1)^2 ]      (V7 = 0; entry sign of D_xxx flips)
"""
import json, glob, os
import numpy as np

def a20c(a):
    den = (a - 1)*(2*a + 1)**2
    return np.nan if abs(den) < 1e-15 else 4*a*(a+1)*(a-2)**2/den

def V7sign(a, a20):
    return np.sign(-150*(a-2)*(-4*a*(a+1)*(a-2)**2 + a20*(a-1)*(2*a+1)**2))

rows = []
for f in sorted(glob.glob('ledger/*.jsonl') + glob.glob('ledger_grid/*.jsonl')
                + glob.glob('ledger_axis/*.jsonl')):
    recs = []
    for line in open(f):
        line = line.strip()
        if not line: continue
        try: recs.append(json.loads(line))
        except Exception: pass
    if not recs: continue
    def g(r, k, d=np.nan):
        try: return float(r[k])
        except Exception: return d
    a = g(recs[0], 'a'); a20 = g(recs[0], 'a20')
    G = np.array([g(r, 'Dxxx') for r in recs])
    x0 = np.array([g(r, 'x0') for r in recs])
    L = np.array([g(r, 'L', np.nan) for r in recs])
    ok = np.isfinite(G)
    G, x0, L = G[ok], x0[ok], L[ok]
    if G.size < 2: continue
    s = np.sign(G)
    changes = int(np.sum(s[:-1]*s[1:] < 0))
    rows.append(dict(f=os.path.basename(f), a=a, a20=a20, n=G.size,
                     x0lo=x0.min(), x0hi=x0.max(), Gstart=G[0], Gend=G[-1],
                     changes=changes, a20c=a20c(a), V7s=V7sign(a, a20),
                     Lmin=np.nanmin(L) if L.size else np.nan))

print(" ledger                              a      a20      a20_c    side  n    x0 range          D_xxx start     D_xxx end      sign changes")
for r in rows:
    side = "below" if (np.isfinite(r['a20c']) and r['a20'] < r['a20c']) else "above"
    print(" %-34s %-6.3g %-8.4g %-8.4g %-5s %-4d %6.3f..%-8.3f %-15.4g %-14.4g %s"
          % (r['f'], r['a'], r['a20'], r['a20c'], side, r['n'], r['x0lo'], r['x0hi'],
             r['Gstart'], r['Gend'],
             ("*** %d ***" % r['changes']) if r['changes'] else "0"))
print()
print("shapes whose ENTRY sign differs from another shape's (the swallow-tail lever):")
for r in rows:
    print("   %-34s sign(V7)=%+d   sign(D_xxx start)=%+d  sign(D_xxx end)=%+d %s"
          % (r['f'], r['V7s'], np.sign(r['Gstart']), np.sign(r['Gend']),
             "  <== ENDS DIFFER" if np.sign(r['Gstart'])*np.sign(r['Gend']) < 0 else ""))
