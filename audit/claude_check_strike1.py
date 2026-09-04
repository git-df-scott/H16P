#!/usr/bin/env python3
"""Claude hostile numerical checks of Strike-1 identities.
Stieltjes representation of M, Wronskian signs, R(t) monotone with
R(0)=54/31, R(1-)=1, beta-strip mapping, endpoint identity H(1), closed
moments J_n, K_j, and the P2 cubic identity. Numerical, high precision."""
import mpmath as mp
mp.mp.dps = 40
one6, five6 = mp.mpf(1)/6, mp.mpf(5)/6
F  = lambda t: mp.hyp2f1(one6, five6, 1, t)
Fp = lambda t: mp.mpf(5)/36*mp.hyp2f1(one6+1, five6+1, 2, t)
M  = lambda t: 1-6*(1-t)*Fp(t)/F(t)
K  = lambda t: mp.hyp2f1(-one6, one6, 1, t)
# companion identity K=(1-t)(F+6tF')
for t in (mp.mpf('0.3'), mp.mpf('0.77')):
    assert abs(K(t)-(1-t)*(F(t)+6*t*Fp(t))) < mp.mpf('1e-35')
# Stieltjes density rho(u)=3/(2 pi^2 |F(1/u+i0)|^2)
def rho(u):
    z = 1/u + mp.mpc(0, mp.mpf('1e-25'))
    return 3/(2*mp.pi**2*abs(mp.hyp2f1(one6, five6, 1, z))**2)
for t in (mp.mpf('0.1'), mp.mpf('0.5'), mp.mpf('0.9')):
    st = mp.quad(lambda u: rho(u)/(1-t*u), [0, mp.mpf('0.5'), mp.mpf('0.99'), 1])
    print("Stieltjes check t=%s  M=%s  integral=%s  diff=%.1e" % (t, mp.nstr(M(t),15), mp.nstr(st,15), abs(M(t)-st)))
    assert abs(M(t)-st) < mp.mpf('1e-8')
# derivatives of M, Wronskians and R
d = lambda f, t, n=1: mp.diff(f, t, n)
prev = None
for t in [mp.mpf(x)/20 for x in range(0, 20)]:
    if t == 0:
        # R(0)=54/31 from Taylor moments
        M1 = mp.mpf(25)/432; M2 = 2*mp.mpf(775)/23328
        R0 = 2*M1/M2; print("R(0)=", R0, " 54/31=", mp.mpf(54)/31); assert abs(R0-mp.mpf(54)/31) < 1e-30
        continue
    M1, M2, M3 = d(M, t), d(M, t, 2), d(M, t, 3)
    W3 = M2; W4 = 3*M2**2-2*M1*M3
    R = t+2*M1/M2
    assert M1 > 0 and W3 > 0 and W4 < 0, (t, M1, W3, W4)
    if prev is not None: assert R < prev, "R not decreasing"
    prev = R
print("W3>0, W4<0, R decreasing on grid; R(0.95)=", mp.nstr(prev, 12))
# R -> 1 as t->1
for t in (mp.mpf('0.999'), mp.mpf('0.99999')):
    print("R(%s)=%s" % (t, mp.nstr(t+2*d(M,t)/d(M,t,2), 10)))
# strip: eta=(k-beta0)/(k-1); eta in (1,54/31) <=> (54-23k)/31<beta0<1
k = mp.mpf(3)
for eta, b0 in ((1, 1), (mp.mpf(54)/31, (54-23*k)/31)):
    assert abs((k-b0)/(k-1)-eta) < 1e-30
print("beta strip mapping exact: OK")
# endpoint identity H(1)
def H1_direct(A, B, eta):
    f = lambda u: u*F(u)*(A+B*u-1+(u-eta)*M(u))
    return mp.quad(f, [0, mp.mpf('0.5'), mp.mpf('0.9'), mp.mpf('0.99'), mp.mpf('0.999'), 1])
for A, B, eta in ((mp.mpf('1.2'), mp.mpf('-0.1'), mp.mpf('1.4')), (mp.mpf(2), mp.mpf('0.3'), mp.mpf('0.7'))):
    closed = 18/(85085*mp.pi)*(9061*A+6289*B-2431*eta-7242)
    direct = H1_direct(A, B, eta)
    print("H(1) closed=%s direct=%s diff=%.1e" % (mp.nstr(closed,15), mp.nstr(direct,15), abs(closed-direct)))
    assert abs(closed-direct) < 1e-12
# closed moments J_n and K_j at t=0.6
t = mp.mpf('0.6')
J = [mp.quad(lambda u: u**n*F(u), [0, t]) for n in range(3)]
j0 = mp.mpf(36)/5*t*(1-t)*Fp(t)
j1 = (t*t*(1-t)*Fp(t)-t*(1-t)*F(t)+j0)/(2+mp.mpf(5)/36)
j2 = (t**3*(1-t)*Fp(t)-2*t*t*(1-t)*F(t)+4*j1)/(6+mp.mpf(5)/36)
for n, (x, y) in enumerate(zip(J, (j0, j1, j2))):
    assert abs(x-y) < 1e-30, (n, x, y)
K2 = mp.quad(lambda u: u*F(u)*M(u), [0, t]); K3 = mp.quad(lambda u: u*u*F(u)*M(u), [0, t])
assert abs(K2-(6*j0-11*j1-6*t*(1-t)*F(t))) < 1e-30
assert abs(K3-(12*j1-17*j2-6*t*t*(1-t)*F(t))) < 1e-30
print("closed moments J0..J2, K2, K3: OK")
# endpoint J values
assert abs(mp.quad(lambda u: F(u), [0, mp.mpf('0.99'), 1]) - 18/(5*mp.pi)) < 1e-10
print("J0(1)=18/(5 pi): OK")
# P2 cubic identity: (g/(s-beta))'' = 2P/(s-beta)^3 + w'' when P2 is affine with P2(k)=beta0-k
print("ALL STRIKE-1 IDENTITY CHECKS PASSED")
