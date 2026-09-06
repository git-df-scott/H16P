#!/usr/bin/env python3
"""Is the minimum of nu interior (a genuine near-swallow-tail that shape motion
could push through zero) or an endpoint/asymptotic decay?"""
import json, glob, os, numpy as np
for f in ['ledger/cusp_row5.jsonl', 'ledger_grid/cusp_c_am2p9_o1p0.jsonl',
          'ledger_grid/cusp_c_am2p5_o0p08.jsonl']:
    recs = [json.loads(l) for l in open(f) if l.strip()]
    def col(k):
        out = []
        for r in recs:
            try: out.append(float(r[k]))
            except Exception: out.append(np.nan)
        return np.array(out)
    nu, x0, H, G = col('nu'), col('x0'), col('Dxxxx'), col('Dxxx')
    L, tr, V1 = col('L'), col('transv'), col('V1')
    m = np.isfinite(nu); nu, x0, H, G, L, tr, V1 = (v[m] for v in (nu, x0, H, G, L, tr, V1))
    k = int(np.argmin(np.abs(nu)))
    print("="*96)
    print("%s   a=%g a20=%g   n=%d" % (os.path.basename(f), col('a')[0], col('a20')[0], nu.size))
    print("  min|nu| = %.6f at index %d of %d  ->  %s"
          % (abs(nu[k]), k, nu.size-1,
             "INTERIOR minimum" if 0 < k < nu.size-1 else "AT AN ENDPOINT (asymptotic decay)"))
    lo, hi = max(0, k-4), min(nu.size, k+5)
    print("   idx    x0         nu          D_xxx        D_xxxx      L         transv    V1")
    for i in range(lo, hi):
        print("   %-5d %-10.6f %-11.6f %+-12.4e %+-11.4e %-9.5f %-9.5f %+.3e"
              % (i, x0[i], nu[i], G[i], H[i], L[i], tr[i], V1[i]))
    print("   tail of nu:", np.array2string(nu[-6:], precision=6))
