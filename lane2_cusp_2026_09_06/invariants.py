#!/usr/bin/env python3
"""Is the escaping branch a genuine geometric degeneration, or motion along the
scaling orbit?

The Cherkas form xdot = 1+xy is preserved by (x,y,t) -> (alpha x, y/alpha, alpha t),
under which
    a -> a,  a01 -> alpha a01,  a11 -> alpha^2 a11,  a10 -> alpha^3 a10,
    a20 -> alpha^4 a20,  a00 -> alpha^2 a00.
Scale invariants:  a,  I1 = a11^2/a20,  I2 = a01^2/a11,  I3 = a10/(a01 a11).
The indicator nu = D_xxx/(D_xxxx r0) is itself invariant (D ~ alpha,
D_xxx ~ alpha^-2, D_xxxx ~ alpha^-3, r0 ~ alpha), so nu -> 0 cannot be a
coordinate artefact -- but it may still be attained only on the BOUNDARY of the
invariant shape space."""
import json, numpy as np

def invs(a, a20, a11, a01, a10):
    return (a11**2/a20 if a20 else np.nan,
            a01**2/a11 if a11 else np.nan,
            a10/(a01*a11) if (a01 and a11) else np.nan)

for f in ['ledger/cusp_row5.jsonl', 'ledger_grid/cusp_c_am2p9_o1p0.jsonl',
          'ledger_grid/cusp_c_am2p5_o0p08.jsonl']:
    recs = [json.loads(l) for l in open(f) if l.strip()]
    col = lambda k: np.array([float(r[k]) for r in recs])
    a, a20 = col('a')[0], col('a20')[0]
    a11, a01, a10, nu, x0 = col('a11'), col('a01'), col('a10'), col('nu'), col('x0')
    print("="*94); print("%s   a=%g a20=%g" % (f, a, a20))
    print("  idx    x0        nu          I1=a11^2/a20   I2=a01^2/a11   I3=a10/(a01 a11)")
    for i in [0, len(nu)//3, 2*len(nu)//3, len(nu)-1]:
        I1, I2, I3 = invs(a, a20, a11[i], a01[i], a10[i])
        print("  %-6d %-9.6f %-11.6f %-14.6g %-14.6g %.6g" % (i, x0[i], nu[i], I1, I2, I3))
    I1e, I2e, I3e = invs(a, a20, a11[-1], a01[-1], a10[-1])
    I1s, I2s, I3s = invs(a, a20, a11[0], a01[0], a10[0])
    print("  invariants move by factors: I1 x%.3g   I2 x%.3g   I3 x%.3g"
          % (I1e/I1s, I2e/I2s, I3e/I3s))
    print("  -> the escape is %s in invariant coordinates"
          % ("GENUINE (invariants move substantially)"
             if max(abs(np.log(abs(I1e/I1s))), abs(np.log(abs(I2e/I2s)))) > 0.3
             else "mostly along the scaling orbit"))
