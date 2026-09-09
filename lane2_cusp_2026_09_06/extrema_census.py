#!/usr/bin/env python3
"""Entry step for the separated double-fold system.

Two double cycles in ONE nest at the SAME parameter require the displacement D
to have TWO distinct stationary points with D = 0 at both.  Since D must be
one-signed between them, that needs D to have >= 3 interior stationary points
(min-max-min), i.e. the rotated-family curve beta* must have three interior
extrema.  With three cycles D has only TWO stationary points.

So: census the zeros and stationary points of D on the nest domain for the
(3,1) controls with a second focus (Cherkas rows 7 and 8) and for row 2."""
import json
import mpmath as mp
from engine import Engine
from cusp import Cusp
mp.mp.dps = 40

ROWS = {
 7: dict(a=mp.mpf(8)/11, a20=-12, a11="2.1502",  a01=mp.mpf(67)/220, a10="-26.5"),
 8: dict(a="1.04", a20=-120, a11="1.51997", a01="1.56", a10="-79.6"),
 2: dict(a="1.5", a20=-15, a11="0.79993", a01="3.2", a10="9.17"),
}
eng = Engine(); print("engine:", eng.banner)
for rid, P in ROWS.items():
    c = Cusp(eng, P["a"], P["a20"], side=1)
    u = [mp.mpf(P["a11"]), mp.mpf(P["a01"]), mp.mpf(P["a10"])]
    print("\n=== Cherkas row %d  (a=%s, a20=%s) ===" % (rid, mp.nstr(mp.mpf(P['a']),6), P["a20"]))
    xs, Ds, Dxs = [], [], []
    x = mp.mpf("1.02")
    while x < mp.mpf("6.5"):
        r = c.val(u, x)
        if r["status"] == "OK":
            xs.append(x); Ds.append(r["D"]); Dxs.append(r["Dx"])
        else:
            if xs: break
        x += mp.mpf("0.04")
    if not xs:
        print("   no return domain reached"); continue
    print("   return domain sampled: x0 in [%s, %s]  (%d points)"
          % (mp.nstr(xs[0], 6), mp.nstr(xs[-1], 6), len(xs)))
    zs = [(xs[i], xs[i+1]) for i in range(len(xs)-1) if Ds[i]*Ds[i+1] < 0]
    st = [(xs[i], xs[i+1]) for i in range(len(xs)-1) if Dxs[i]*Dxs[i+1] < 0]
    print("   D  sign changes (cycles):        %d  at %s"
          % (len(zs), [ "%.3f" % float((a+b)/2) for a, b in zs]))
    print("   D_x sign changes (stationary):   %d  at %s"
          % (len(st), [ "%.3f" % float((a+b)/2) for a, b in st]))
    for a_, b_ in st:
        m = (a_+b_)/2; r = c.val(u, m)
        print("      stationary near x0=%.4f : D = %-14s D_xx = %s"
              % (float(m), mp.nstr(r["D"], 6), mp.nstr(r["Dxx"], 6)))
print("\ncalls:", eng.ncalls); eng.close()
