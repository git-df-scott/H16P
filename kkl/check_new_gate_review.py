"""Bounded exact-only reviewer checks. No orbit integration or cycle count."""
import os
import resource
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
resource.setrlimit(resource.RLIMIT_CPU, (10, 10))
os.nice(10)

import time
import sympy as S

started = time.process_time()
c, x, y, u, s, m = S.symbols('c x y u s m')
Q = S.Rational
J = 305 + 634*c - 11*c**2 - 1000*c**3
d, e = 16-10*c, 1+2*c
m0 = 210/(11*c-5)
sigma = 5*e/21
den = 1000*c*c+1231*c+145
n11 = 1000*c*c+769*c+355
num = S.Matrix([[n11,-10*e*(11*c-5)],[-2100*e,-n11]])
M = num/den
M_formula = S.Matrix([[m0*sigma**2-1, -2*sigma],
                      [-2*m0*sigma, 1-m0*sigma**2]])/(1+m0*sigma**2)
assert all(S.cancel(v)==0 for v in M-M_formula)
assert all(S.expand(v)==0 for v in num*num-den**2*S.eye(2))
assert S.expand(num.det()+den**2)==0
z = S.Matrix([x,y])
# Clear the common denominator before expanding the quadratic identity.
linear = S.Matrix([[0,11*c-5],[-210,0]])
assert all(S.expand(v)==0 for v in linear*num+num*linear)
quadratic = S.Matrix([x*x+x*y,-10*x*x+Q(11,5)*x*y+c*y*y])
nz = num*z
residual = quadratic.subs({x:nz[0],y:nz[1]},simultaneous=True)+den*num*quadratic
for entry in residual:
    for coefficient in S.Poly(S.expand(entry), x, y).coeffs():
        assert S.rem(coefficient, J, c) == 0
assert J.subs(c,Q(241,250)) == Q(25281,2500)
assert J.subs(c,Q(39,40)) == -Q(11333,800)

F = 10*s+Q(11,5)*s*s/(s-1)-c*s**3/(s-1)**2
sH = 21/d
mH = S.factor(F.subs(s,sH))
assert S.cancel(mH-21*(1000*c*c+1021*c+481)/(25*d*e*e)) == 0
KH = S.factor(mH*(Q(11,5)*c-1)-42)
assert S.cancel(KH+441*J/(125*d*e*e)) == 0

W = m+(2*m+10)*x+(m+Q(111,5))*x*x+(Q(61,5)-c)*x**3
B = x*(21+d*x)/(5*(1+x))
C = x*W/(1+x)
A = (c+1)/(1+x)
Fs = S.diff(F,s)
assert S.cancel(S.diff(C,x).subs({m:F,x:-s})-s*(s-1)*Fs) == 0
assert S.cancel(S.diff(C,x,2).subs({m:F,x:-s})+2*(2*s-1)*Fs+s*(s-1)*S.diff(Fs,s)) == 0
assert S.expand(5*e**2*(1-c)*(61-5*c)-c*(c+1)*d*d-J) == 0
# The original numerator reduces to B'(x_H)*s_H*2J/(25e^2).
BpH = 21*d/(25*e)
assert S.cancel(BpH*sH*2*J/(25*e*e)-882*J/(625*e**3)) == 0

N = S.expand((d*(1+x)+(c+1)*(21+d*x))*W-(1+x)*(21+d*x)*S.diff(W,x))
point = {c:Q(1001,1000), m:Q(196,5)}
Nu = S.Poly(S.expand(N.subs(point).subs(x,u-1)),u)
expected = (3006504501,592888296,-920742694,-122238304,6708201)
assert [Nu.nth(j) for j in range(5)] == [Q(v,10**8) for v in expected]
assert Nu.eval(1) == Q(32039,1250)
assert Nu.eval(4) == -Q(15459777419,10**8)
assert Nu.eval(64) == Q(76770488465461,10**8)
assert S.cancel(F.subs({c:point[c],s:Q(88,25)})-point[m]) == -Q(102412,1771875)
assert S.cancel(F.subs({c:point[c],s:Q(353,100)})-point[m]) == Q(32990493,581900000)
assert Q(88,25) > S.cancel(sH.subs(c,point[c]))
Kpoint = S.cancel(point[m]*(Q(11,5)*point[c]-1)-42)
assert Kpoint == Q(32039,6250)

print('PASS exact original-coordinate reversing involution modulo J')
print('PASS remote trace threshold and original-coordinate Hopf numerator')
print('PASS rational c>1 checkpoint isolation values and multiplier quartic')
print('No ODEs, numerical discovery, return validation, or cycle-count claim')
print('CPU seconds:', round(time.process_time()-started, 6))
