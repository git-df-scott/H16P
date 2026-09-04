#!/usr/bin/env python3
"""Claude hostile numerical audit of the Strike-3 fixed-lambda obstruction.
Along gamma(1-eps), with k=1/(lambda*eps), compute P(1-eps) by direct
quadrature (no ODE) and compare 2304*pi*P/(eps^{5/6} L) with the claimed
limit  P_lambda(1) = D*B_lambda(1) - int_1^inf omega_lambda (e-Vv+Qv log v)."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import mpmath as mp
from claude_green_tools import *
mp.mp.dps = 45
D, e, V, Qc = mp.mpf(30)/77, -mp.mpf(15)/154, mp.mpf(45)/154, -mp.mpf(45)/(154*mp.log(2))
CO = mp.mpf(3)/(5*mp.cbrt(2))
def claimed_limit(lam, c=1):
    lam = mp.mpf(lam)
    def ytil(v): return lam**(mp.mpf(5)/6)/CO*Lift.O(mp.asinh(mp.sqrt(v/lam)))
    def omega(v): return ytil(v)/((v+lam)**mp.mpf(1.5)*v**mp.mpf(1.5))
    x = mp.asinh(mp.sqrt(c/lam))
    JO = mp.mpf(9)/50*mp.cosh(5*x/3)+mp.mpf(9)/2*mp.cosh(x/3)
    Bl = 2*lam**(mp.mpf(5)/6)/CO*(mp.mpf(23)/18*JO-Lift.O(x)*mp.tanh(x)+Lift.Ox(x)/2)
    # sanity: B' = omega c^2
    dB = mp.diff(lambda cc: 2*lam**(mp.mpf(5)/6)/CO*(mp.mpf(23)/18*(mp.mpf(9)/50*mp.cosh(5*mp.asinh(mp.sqrt(cc/lam))/3)+mp.mpf(9)/2*mp.cosh(mp.asinh(mp.sqrt(cc/lam))/3))-Lift.O(mp.asinh(mp.sqrt(cc/lam)))*mp.tanh(mp.asinh(mp.sqrt(cc/lam)))+Lift.Ox(mp.asinh(mp.sqrt(cc/lam)))/2), c)
    assert abs(dB-omega(c)*c*c) < mp.mpf('1e-25'), (dB, omega(c)*c*c)
    tail = mp.quad(lambda v: omega(v)*(e-V*v+Qc*v*mp.log(v)), [c, 10, 100, 1000, mp.inf])
    return D*Bl-tail, D*Bl, -tail
for lam in ("0.5", "2"):
    lim, b, tl = claimed_limit(lam)
    print(f"lambda={lam}: claimed limit P_lambda(1)={mp.nstr(lim,10)}  (D*B={mp.nstr(b,8)}, -tail={mp.nstr(tl,8)})")
    for ee in ("1e-2", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7"):
        eps = mp.mpf(ee); r = 1-eps; L = mp.log(432/eps)
        co = coefficients_from_r(r); k = 1/(mp.mpf(lam)*eps); lift = Lift(k)
        P, Phi, Y0, P0 = P_Phi_at(lift, co, r)
        scaled = 2304*mp.pi*P/(eps**(mp.mpf(5)/6)*L)
        print(f"  eps={ee} k={mp.nstr(k,6)} P(1-eps)={mp.nstr(P,8)} scaled={mp.nstr(scaled,8)}  Phi(1-eps)={mp.nstr(Phi,8)} Y0={mp.nstr(Y0,8)}")
