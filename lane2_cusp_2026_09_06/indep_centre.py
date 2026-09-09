#!/usr/bin/env python3
"""PROTOCOL rule 2: confirm the centre on the INDEPENDENT engine.

indep_engine.py is a from-scratch mpmath Taylor integrator in Cartesian (x,y)
with its own step rule and event solver, sharing no code with cusp_engine.cpp.
If D is ~0 across the whole annulus on BOTH engines, the Newton limit is a
centre, not a quadruple cycle."""
import mpmath as mp
import indep_engine as IE
mp.mp.dps = 40

a   = mp.mpf(-2)
a20 = -mp.mpf(3254)/675
a11 = mp.mpf("8.067217612039481856440047")
a01 = -a11 - 3
a10 = mp.mpf("16.53363941409600843877444")
print("a=%s a20=%s (= -3254/675)" % (mp.nstr(a,4), mp.nstr(a20,20)))
print("a11=%s\na01=%s\na10=%s" % (mp.nstr(a11,22), mp.nstr(a01,22), mp.nstr(a10,22)))
print("\nindependent engine (mpmath Taylor, Cartesian, own event solver):")
print("   x0          D (independent)        |D|/x0-scale")
for s in ("1.01", "1.03", "1.0337", "1.06", "1.12", "1.20"):
    x0 = mp.mpf(s)
    try:
        d = IE.D(a, a20, a11, a01, a10, x0)
        dv = d[0] if isinstance(d, tuple) else (d["D"] if isinstance(d, dict) else d)
        print("   %-11s %-22s %s" % (s, mp.nstr(dv, 8), mp.nstr(abs(dv)/(x0-1), 6)))
    except Exception as e:
        print("   %-11s ERROR %s" % (s, repr(e)[:70]))
