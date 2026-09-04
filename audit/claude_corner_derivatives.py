#!/usr/bin/env python3
"""Claude: large-k corner linearization. At the limit point, compute
P_*(1;k), Phi_*(1;k) and the partial derivatives of P(1) and Phi(1) w.r.t.
(A,B,eta), then the derivative along the universal deviation direction
dir=(-643/462, 1105/462, 1) (eta-1 = v>0 along any late-root path)."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import mpmath as mp
from claude_green_tools import *
mp.mp.dps = 40
t_end = 1-mp.mpf('1e-14')
def basisH(j):
    return lambda u: primitive_basis_closed(u)[j]
for k in (8.5, 12, 20, 50, 100, 300, 1000, 3000, 10000):
    lift = Lift(k)
    pts = breakpoints(t_end)
    # partials: H = (A-1)K0 + B K1 - eta K2 + K3 ; Y0 partials from formula; P0 partial = dY1 - r dY0
    dY0 = [3*mp.mpf(c)/1361360 for c in (1326, 864, -2431)]
    dP0 = [-(mp.mpf(3)/2*(1+lift.a)+lift.r())*d for d in dY0]; dP0[2] -= mp.mpf(1)/192
    V = [mp.quad(lambda u: lift.Omega(u)*basisH(j)(u), pts) for j in range(3)]
    W = [mp.quad(lambda u: lift.Rcal(u)*lift.Omega(u)*basisH(j)(u), pts) for j in range(3)]
    sgn = (1, 1, -1)
    dP1 = [dP0[j]-sgn[j]*V[j] for j in range(3)]
    dPhi1 = [dY0[j]+sgn[j]*W[j] for j in range(3)]
    dirv = (-mp.mpf(643)/462, mp.mpf(1105)/462, mp.mpf(1))
    DP = sum(d*x for d, x in zip(dirv, dP1)); DPhi = sum(d*x for d, x in zip(dirv, dPhi1))
    P1, Phi1, Y0, P0 = P_Phi_at(lift, (mp.mpf(94)/77, -mp.mpf(17)/77, mp.mpf(1)), t_end, Hfun=Hstar)
    print(f"k={k:>6}: P*(1)={mp.nstr(P1,7)}  Phi*(1)={mp.nstr(Phi1,7)}  dP/dv={mp.nstr(DP,7)}  dPhi/dv={mp.nstr(DPhi,7)}  v needed for P(1)=0: {mp.nstr(-P1/DP,6)}  Phi there: {mp.nstr(Phi1-P1/DP*DPhi,7)}")
