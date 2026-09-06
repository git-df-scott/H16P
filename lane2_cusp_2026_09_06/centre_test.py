#!/usr/bin/env python3
"""Is the Newton limit a CENTRE (D identically zero on the annulus), rather than
a quadruple cycle?  Test D over a RANGE of x0 at the converged u, on both the
degree-4 engine and the independent binary128 engine (PROTOCOL rule 2)."""
import mpmath as mp
from engine import Engine
from cusp import Cusp
mp.mp.dps = 50

a, a20 = "-2.000000000000000000000000000000000", "-4.820740740740740740740740740740741"
u = [mp.mpf("8.067217612039481856440047"),
     mp.mpf("-11.06721761203948185644005"),
     mp.mpf("16.53363941409600843877444")]
print("converged u:  a11=%s\n              a01=%s\n              a10=%s"
      % (mp.nstr(u[0], 25), mp.nstr(u[1], 25), mp.nstr(u[2], 25)))
print("\nexact-looking relations at this point (a = -2):")
print("   a01 + a11 + 3   = %s" % mp.nstr(u[1] + u[0] + 3, 8))
print("   V1 = a11+a01-2a-1 = %s   <- trace at A vanishes" % mp.nstr(u[0]+u[1]-2*mp.mpf(a)-1, 8))

eng = Engine()
c = Cusp(eng, a, a20, side=1)
print("\nD over a RANGE of x0 at that fixed u (engine: %s):" % eng.banner)
print("   x0          D                    D_x                  D_xx                 D_xxx                D_xxxx")
for x0 in ("1.005", "1.02", "1.0337", "1.05", "1.10", "1.20", "1.35"):
    r = c.val(u, mp.mpf(x0))
    if r["status"] != "OK":
        print("   %-11s %s" % (x0, r["status"])); continue
    print("   %-11s %-20s %-20s %-20s %-20s %s"
          % (x0, mp.nstr(r["D"], 6), mp.nstr(r["Dx"], 6), mp.nstr(r["Dxx"], 6),
             mp.nstr(r["Dxxx"], 6), mp.nstr(r["Dxxxx"], 6)))
print("\ncalls:", eng.ncalls); eng.close()
