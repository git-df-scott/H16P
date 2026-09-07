#!/usr/bin/env python3
"""The correct swallow-tail detector.

A nondegenerate multiplicity-four cycle needs D_xxx = 0 AND D_xxxx != 0, i.e.
a zero of the scale-free indicator nu = D_xxx/(D_xxxx r0).  A sign change of
D_xxx alone is NOT sufficient: it can be the whole jet changing sign, which is
what ledger cusp_c_am2p0_om0p08 does.

Scan every saved cusp curve on nu."""
import json, glob, os, numpy as np
print(" ledger                              a      a20      n    nu range              min|nu|   sign chg of nu   min|D_xxxx|")
best = []
for f in sorted(glob.glob('ledger/*.jsonl') + glob.glob('ledger_grid/*.jsonl')
                + glob.glob('ledger_axis/*.jsonl')):
    recs = [json.loads(l) for l in open(f) if l.strip()]
    if not recs: continue
    def col(k):
        out = []
        for r in recs:
            try: out.append(float(r[k]))
            except Exception: out.append(np.nan)
        return np.array(out)
    nu = col('nu'); H = col('Dxxxx'); a = col('a')[0]; a20 = col('a20')[0]
    m = np.isfinite(nu)
    if m.sum() < 2: continue
    nu = nu[m]; H = H[m]
    s = np.sign(nu); chg = int(np.sum(s[:-1]*s[1:] < 0))
    print(" %-34s %-6.3g %-8.4g %-4d [%+.6f, %+.6f]  %.5f   %-16s %.3g"
          % (os.path.basename(f), a, a20, nu.size, nu.min(), nu.max(),
             np.abs(nu).min(), ("*** %d ***" % chg) if chg else "0", np.abs(H).min()))
    best.append((np.abs(nu).min(), os.path.basename(f)))
print()
best.sort()
print("closest approaches to a swallow-tail (min |nu|):")
for v, n in best[:5]:
    print("   %-34s min|nu| = %.6f" % (n, v))
print()
print("A swallow-tail requires nu = 0.  The global minimum over every saved cusp")
print("curve is %.6f, attained on %s." % (best[0][0], best[0][1]))
