#!/usr/bin/env python3
"""Claude control map for Strike-4 viability.
Along gamma(1-eps) for a grid of (eps, kappa), compute by direct quadrature:
  P(tau1),P(tau2),P(tau3) [S1 needs P0>0, P(tau1)<0, P(tau2)>0, P(tau3)<0]
  Phi(tau1)=Y0+int_0^{tau1} Rcal*Omega*H  [necessary for S2: Phi(tau1)>0]
Also the limit-point quantities P_*(1;k), Phi_*(1;k), Z_*(1;k)."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import mpmath as mp
from claude_green_tools import *
mp.mp.dps = 40
print("=== limit point (94/77,-17/77,1): P_*(1;k), Phi_*(1;k) ===")
co_star = (mp.mpf(94)/77, -mp.mpf(17)/77, mp.mpf(1))
for k in (3, 5, 8, 8.5, 9, 12, 16, 30, 100, 1000, 10000):
    lift = Lift(k)
    P1, Phi1, Y0, P0 = P_Phi_at(lift, co_star, 1-mp.mpf('1e-14'), Hfun=Hstar)
    print(f"k={k:>7}: P0*={mp.nstr(P0,8)}  P*(1)={mp.nstr(P1,8)}  Phi*(1)={mp.nstr(Phi1,8)}  Phi*(1)/|Y0|={mp.nstr(Phi1/abs(Y0),8)}")
print("=== path gamma(1-eps): S1 signs and Phi(tau1) ===")
for ee in ("1e-2", "1e-3", "1e-4", "1e-5", "1e-6"):
    eps = mp.mpf(ee); r = 1-eps; co = coefficients_from_r(r); anchors = threshold_anchors(r)
    for k in (5, 8, 8.5, 9, 12, 20, 50, 200):
        lift = Lift(k)
        out = []
        for t in anchors:
            P, Phi, Y0, P0 = P_Phi_at(lift, co, t)
            out.append((P, Phi))
        s1 = P0 > 0 and out[0][0] < 0 and out[1][0] > 0 and out[2][0] < 0
        print(f"eps={ee} k={k:>5}: P0={mp.nstr(P0,6)} P(tau)=({mp.nstr(out[0][0],5)},{mp.nstr(out[1][0],5)},{mp.nstr(out[2][0],5)}) S1={s1}  Phi(tau1)={mp.nstr(out[0][1],7)} Phi/|Y0|={mp.nstr(out[0][1]/abs(Y0),6)}")
