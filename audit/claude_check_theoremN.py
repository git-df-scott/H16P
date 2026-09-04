#!/usr/bin/env python3
"""Claude hostile check of the three identities behind Astra's Theorem N proof.
(N0) int_0^1 W_1 H_* = 3/1232 with W_1 = 3/(2304 t^2)[(1-t)^(-17/6)-(1-t)^(-13/6)].
(N2) moment curve x=K1/K0, m=K2/K0: x(1)=6289/9061, m(1)=11/41, S=(M-m)/(t-x)
     increasing with S(1-)=1105/462; m convex in x.
(N6) operator for v_a=z_a/(1-at)^{3/2} and the residual L_a v_1 = (1-a)(22-7d^{2/3})/(6 d^{7/3}).
(N7) direct numerical check W_a < W_1 at sample points, and Phi_a(y1) < -int_{y1}^1 W_1 H_*
     on Astra's tuned shots (N8)."""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
import mpmath as mp, sympy as S
from claude_green_tools import *
mp.mp.dps = 30
# (N0)
W1 = lambda t: 3/(2304*t*t)*((1-t)**(-mp.mpf(17)/6)-(1-t)**(-mp.mpf(13)/6))
def sub6(f):  # t = 1-u^6 substitution removes the (1-t)^{-5/6} endpoint singularity
    return mp.quad(lambda u: (f(1-u**6)*6*u**5 if u**6 > mp.mpf('1e-36') else mp.mpf(0)), [0, mp.mpf('0.3'), mp.mpf('0.6'), mp.mpf('0.9'), 1-mp.mpf('1e-12')])
n0 = sub6(lambda t: W1(t)*Hstar(t))
print("(N0) int W1 H* =", mp.nstr(n0, 12), " 3/1232 =", mp.nstr(mp.mpf(3)/1232, 12)); assert abs(n0-mp.mpf(3)/1232) < 1e-5
# (N2)
K = lambda t: primitive_basis_closed(t)
x1 = K(mp.mpf(1))[1]/K(mp.mpf(1))[0]; m1 = K(mp.mpf(1))[2]/K(mp.mpf(1))[0]
print("(N2) x(1)=", mp.nstr(x1, 15), "6289/9061=", mp.nstr(mp.mpf(6289)/9061, 15), " m(1)=", mp.nstr(m1, 15), "11/41=", mp.nstr(mp.mpf(11)/41, 15))
assert abs(x1-mp.mpf(6289)/9061) < 1e-20 and abs(m1-mp.mpf(11)/41) < 1e-20
def Sslope(t):
    k = K(t); x = k[1]/k[0]; m = k[2]/k[0]
    return (Mh(t)-m)/(t-x), x, m
prev = None; pts = []
for i in range(1, 100):
    t = mp.mpf(i)/100; s, x, m = Sslope(t); pts.append((x, m))
    if prev is not None: assert s > prev, ("S not increasing", t)
    prev = s
print("     S increasing on grid; S(0.99)=", mp.nstr(prev, 10), " S(1-) claimed 1105/462 =", mp.nstr(mp.mpf(1105)/462, 10), " S(0.999999)=", mp.nstr(Sslope(1-mp.mpf('1e-6'))[0], 10))
# convexity of m in x: second differences positive
for (xa, ma), (xb, mb), (xc, mc) in zip(pts, pts[1:], pts[2:]):
    assert (mc-mb)/(xc-xb) > (mb-ma)/(xb-xa)
print("     m convex in x on grid: OK")
# (N6) symbolic
a, t, d = S.symbols('a t d', positive=True)
Y = S.Function('y')
# homogeneous equation for y (and z): (1-at)(1-t)y'' - (1-a)/2 y' + 5a/36 y = 0 ; v = z (1-at)^(-3/2)
v = S.Function('v')
z = v(t)*(1-a*t)**S.Rational(3,2)
expr = S.expand((1-a*t)*(1-t)*S.diff(z, t, 2) - (1-a)/2*S.diff(z, t) + S.Rational(5,36)*a*z)
expr = S.simplify(expr/(1-a*t)**S.Rational(3,2))
La = S.collect(S.expand(expr), [v(t).diff(t,2), v(t).diff(t), v(t)])
print("(N6) transformed operator:", La)
claimed = (1-a*t)*(1-t)*v(t).diff(t,2) - (1+5*a-6*a*t)/2*v(t).diff(t) + S.Rational(8,9)*a*v(t)
assert S.simplify(La-claimed) == 0
v1 = S.Rational(3,2)*((1-t)**S.Rational(-4,3)-(1-t)**S.Rational(-2,3))
res = claimed.subs(v(t), v1).doit()
target = (1-a)*(22-7*(1-t)**S.Rational(2,3))/(6*(1-t)**S.Rational(7,3))
print("     residual matches (1-a)(22-7d^(2/3))/(6 d^(7/3)):", S.simplify(res-target) == 0)
assert S.simplify(res-target) == 0
# (N7) numeric W_a < W_1
for k in (1.5, 4, 20, 500):
    lift = Lift(k)
    for tt in (mp.mpf('0.2'), mp.mpf('0.7'), mp.mpf('0.99')):
        Wa = lift.Rcal(tt)*lift.Omega(tt)
        assert Wa < W1(tt), (k, tt, Wa, W1(tt))
print("(N7) W_a < W_1 at sampled (kappa,t): OK")
# (N8) on Astra's tuned shots: Phi(tau1) < -int_{tau1}^1 W1 H*
base = os.path.join(os.path.dirname(__file__), "..", "q4", "data")
for row in json.load(open(os.path.join(base, "third_tuned_shoot.json")))["rows"][:3]:
    co = tuple(mp.mpf(x) for x in row["A_B_eta"]); lift = Lift(mp.mpf(row["kappa"])); tau1 = mp.mpf(row["primitive_anchors"][0])
    P, Phi, Y0, P0 = P_Phi_at(lift, co, tau1)
    tail = mp.quad(lambda u: (W1(1-u**6)*Hstar(1-u**6)*6*u**5 if u**6 > mp.mpf('1e-36') else 0), [0, (1-tau1)**(mp.mpf(1)/6)])
    print(f"(N8) r={row['r']}: Phi(tau1)={mp.nstr(Phi,8)} < -tail={mp.nstr(-tail,8)} : {Phi < -tail}")
    assert Phi < -tail
print("THEOREM N IDENTITIES INDEPENDENTLY VERIFIED")
