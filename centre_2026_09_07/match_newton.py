#!/usr/bin/env python3
"""Does k = sqrt(1627)/5 reproduce the numerical Newton limit exactly?"""
import mpmath as mp, sympy as S
mp.mp.dps = 40
k = mp.sqrt(1627)/5
newton = dict(a11=mp.mpf("8.067217612039481856440047"),
              a01=mp.mpf("-11.06721761203948185644005"),
              a10=mp.mpf("16.53363941409600843877444"),
              a20=-mp.mpf(3254)/675)
exact = dict(a11=k, a01=-k-3, a20=-2*k**2/27, a10=5*k*(k+3)/27)
print("k = sqrt(1627)/5 = %s\n" % mp.nstr(k, 30))
print("  coeff   exact (from k)                    Newton limit                      difference")
for nm in ("a11", "a01", "a20", "a10"):
    print("  %-7s %-34s %-34s %s"
          % (nm, mp.nstr(exact[nm], 25), mp.nstr(newton[nm], 25),
             mp.nstr(exact[nm] - newton[nm], 6)))
print()
kk = S.sqrt(1627)/5
print("exact a20 = -2k^2/27 =", S.nsimplify(S.simplify(-2*kk**2/27)), " = -3254/675 ->",
      S.simplify(-2*kk**2/27 + S.Rational(3254,675)) == 0)
print("exact a10 = 5k(k+3)/27 =", S.simplify(5*kk*(kk+3)/27))
print("   equals 1627/135 + sqrt(1627)/9 ->",
      S.simplify(5*kk*(kk+3)/27 - (S.Rational(1627,135) + S.sqrt(1627)/9)) == 0)
print()
print("3 < k < 9 ?  k = %s  ->  %s" % (mp.nstr(k, 12), 3 < k < 9))
print("det J(A) = (k-3)(9-k)/27 = %s  > 0" % mp.nstr((k-3)*(9-k)/27, 12))
