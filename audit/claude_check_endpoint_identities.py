#!/usr/bin/env python3
"""Claude checks of Strike-3 endpoint identities: coefficient limit,
H_* closed form (G1), beta moments (G2), variation-of-parameters identity
Y = Phi*y + P*y2 against the repository ODE solution, and Rcal closed form."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
import mpmath as mp
from claude_green_tools import *
mp.mp.dps = 40
# coefficient limit along the path
for ee in ("1e-4", "1e-8", "1e-12"):
    co = coefficients_from_r(1-mp.mpf(ee))
    print(f"eps={ee}: (A,B,eta)=({mp.nstr(co[0],10)},{mp.nstr(co[1],10)},{mp.nstr(co[2],10)})  target (1.220779,-0.220779,1)")
# (G1): H at the limit point equals H_*
co_star = (mp.mpf(94)/77, -mp.mpf(17)/77, mp.mpf(1))
for t in (mp.mpf('0.2'), mp.mpf('0.7'), mp.mpf('0.99')):
    assert abs(Hval(co_star, t)-Hstar(t)) < mp.mpf('1e-15'), t  # closed-moment evaluator loses ~20 digits to cancellation
print("(G1) H_* closed form: OK")
# (G2) beta moments with Omega_1=(1-t)^{-13/6}/(1152 t^2)
Om1 = lambda t: (1-t)**(-mp.mpf(13)/6)/(1152*t*t)
# substitute t=1-u^6 to remove the endpoint singularity (1-t)^{-5/6}
def sub(f): return mp.quad(lambda u: (f(1-u**6)*6*u**5 if u**6 > mp.mpf('1e-36') else mp.mpf(0)), [0, mp.mpf('0.3'), mp.mpf('0.6'), mp.mpf('0.9'), 1-mp.mpf('1e-12')])
m1 = sub(lambda t: Om1(t)*Hstar(t)); m2 = sub(lambda t: (1-t)**(-mp.mpf(2)/3)*Om1(t)*Hstar(t))
print("(G2) moments:", mp.nstr(m1*14784, 12), "(=1?)", mp.nstr(m2*14784, 12), "(=25?)")
assert abs(m1*14784-1) < 1e-8 and abs(m2*14784-25) < 2e-3  # second integrand ~(1-t)^(-5/6)*log: slow quadrature
# Rcal closed form vs integral, and Y = Phi*y + P*y2 vs repository ODE
from q4_reconstruction import reconstruct
import numpy as np
k = mp.mpf('2.5'); lift = Lift(k); a = float(lift.a)
t0 = mp.mpf('0.6')
R_int = mp.quad(lambda u: 1/(mp.sqrt((1-u)/(1-lift.a*u))*lift.y(u)**2), [0, t0])
print("Rcal closed", mp.nstr(lift.Rcal(t0), 15), "integral", mp.nstr(R_int, 15)); assert abs(R_int-lift.Rcal(t0)) < 1e-25
co = (mp.mpf('1.21'), mp.mpf('-0.13'), mp.mpf('1.21'))
sol = reconstruct(a, *map(float, co), t_end=0.95)
for t in (mp.mpf('0.3'), mp.mpf('0.6'), mp.mpf('0.9')):
    P, Phi, Y0, P0 = P_Phi_at(lift, co, t)
    Y_vp = Phi*lift.y(t)+P*lift.y(t)*lift.Rcal(t)
    Y_ode = sol.sol(float(t))[1]
    print(f"t={t}: Y(var.par.)={mp.nstr(Y_vp,12)}  Y(ODE)={Y_ode:.12e}")
    assert abs(Y_vp-Y_ode) < 1e-9
print("Variation-of-parameters decomposition Y=Phi*y+P*y2: OK")
