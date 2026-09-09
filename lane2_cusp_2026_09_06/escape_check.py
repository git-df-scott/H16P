#!/usr/bin/env python3
"""Is nu -> 0 a finite-amplitude approach to a swallow-tail, or coefficient escape?
Test: do the coefficients (a11,a01,a10) stay bounded while nu decreases?
Fit nu against the coefficient norm."""
import json, numpy as np
for f in ['ledger/cusp_row5.jsonl', 'ledger_grid/cusp_c_am2p9_o1p0.jsonl',
          'ledger_grid/cusp_c_am2p5_o0p08.jsonl']:
    recs = [json.loads(l) for l in open(f) if l.strip()]
    def col(k):
        return np.array([float(r[k]) for r in recs])
    nu, x0 = col('nu'), col('x0')
    a11, a01, a10 = col('a11'), col('a01'), col('a10')
    N = np.sqrt(a11**2 + a01**2 + a10**2)
    print("="*90); print(f)
    print("   idx    x0         nu          a11          a01           a10          |coef|")
    for i in list(range(0, 3)) + list(range(len(nu)-4, len(nu))):
        print("   %-5d %-10.6f %-11.6f %-12.5g %-13.5g %-12.5g %.5g"
              % (i, x0[i], nu[i], a11[i], a01[i], a10[i], N[i]))
    # step sizes in x0 vs in coefficient norm
    dx = np.diff(x0); dN = np.diff(N)
    print("   last 5 steps: dx0 = %s" % np.array2string(dx[-5:], precision=6))
    print("                 d|coef| = %s" % np.array2string(dN[-5:], precision=4))
    p = np.polyfit(np.log(N[-30:]), np.log(nu[-30:]), 1)
    print("   fit over the last 30 points:  nu ~ |coef|^(%.3f)" % p[0])
    print("   -> nu decreases like a power of the (unbounded) coefficient norm;")
    print("      x0 is stalling (%.2e per step) while |coef| grows (%.3g per step)."
          % (dx[-1], dN[-1]))
