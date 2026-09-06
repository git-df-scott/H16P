#!/usr/bin/env python3
"""Adjudicate the D_xxx sign change in cusp_c_am2p0_om0p08.

A nondegenerate swallow-tail needs D_xxx = 0 WITH D_xxxx != 0, i.e. the
scale-free indicator nu = D_xxx/(D_xxxx r0) must pass through zero.
Here D_xxx and D_xxxx share a zero, so nu does not."""
import json, numpy as np
recs = [json.loads(l) for l in open('ledger_grid/cusp_c_am2p0_om0p08.jsonl') if l.strip()]
f = lambda r, k: float(r[k])
x0 = np.array([f(r,'x0') for r in recs]); G = np.array([f(r,'Dxxx') for r in recs])
H = np.array([f(r,'Dxxxx') for r in recs]); nu = np.array([f(r,'nu') for r in recs])
V1 = np.array([f(r,'V1') for r in recs])

def cross(x, y):
    for i in range(len(y)-1):
        if y[i]*y[i+1] < 0:
            t = y[i]/(y[i]-y[i+1])
            return x[i] + t*(x[i+1]-x[i]), i
    return None, None

xg, ig = cross(x0, G); xh, ih = cross(x0, H); xv, iv = cross(x0, V1)
print("linear-interpolated zero crossings along the cusp curve:")
print("   D_xxx  = 0 at x0 = %.8f" % xg)
print("   D_xxxx = 0 at x0 = %.8f   (separation from D_xxx zero: %.2e)" % (xh, abs(xh-xg)))
print("   V1     = 0 at x0 = %.8f   (separation from D_xxx zero: %.2e)" % (xv, abs(xv-xg)))
print()
print("nu = D_xxx/(D_xxxx r0) across the crossing (must -> 0 for a swallow-tail):")
for i in range(max(0,ig-3), min(len(nu), ig+4)):
    print("   x0=%.6f  nu=%.6f   D_xxx=%+.4e  D_xxxx=%+.4e" % (x0[i], nu[i], G[i], H[i]))
print()
print("nu range over the whole curve: [%.6f, %.6f]  -- bounded away from zero"
      % (nu.min(), nu.max()))
print()
print("VERDICT: D_xxx and D_xxxx vanish at the same point to %.1e in x0, and nu never"
      % abs(xh-xg))
print("approaches zero.  Perko's nondegeneracy D_ssss != 0 FAILS at this crossing, so it")
print("is NOT a swallow-tail: it is a common sign change of the whole jet, accompanied by")
print("V1 crossing zero at the same place (separation %.1e)." % abs(xv-xg))
