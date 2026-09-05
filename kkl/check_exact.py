"""Small exact algebra replay, with one thread and a ten-second CPU fuse."""
import os
import resource
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
resource.setrlimit(resource.RLIMIT_CPU,(10,10))
import sympy as s
from geometry import gates

x,y,v,c,m=s.symbols('x y v c m')
u=1+x
P=y+x*x+x*y
Q=-10*x*x+s.Rational(11,5)*x*y+c*y*y-m*x
W=m+(2*m+10)*x+(m+s.Rational(111,5))*x*x+(s.Rational(61,5)-c)*x**3
B=x*(21+(16-10*c)*x)/(5*u)
ddx=s.expand(s.diff(P,x)*P+s.diff(P,y)*Q)
assert s.factor(ddx.subs(y,(v-x*x)/u)-((c+1)*v*v/u+B*v-x*W/u))==0
assert s.factor(s.diff(P,x)+s.diff(Q,y)-B-(1+2*c)*P/u)==0
d=16-10*c
N=s.expand((d*u+(c+1)*(21+d*x))*W-u*(21+d*x)*s.diff(W,x))
K=m*(s.Rational(11,5)*c-1)-42
assert s.factor(N.subs(x,0)-5*K)==0
star=s.Poly(N.subs({c:s.Rational(7,10),m:80}),x)
assert star.as_expr()==6+s.Rational(753,5)*x+s.Rational(7821,25)*x*x+s.Rational(12291,100)*x**3-s.Rational(621,20)*x**4
assert star.eval(-1)==s.Rational(357,25)
assert star.eval(-s.Rational(1,2))==-s.Rational(13431,1600)
on_hopf=s.Poly(s.factor(N.subs(m,210/(11*c-5))),x)
delta=s.factor(on_hopf.nth(2)-d*on_hopf.nth(1)/7)
J=305+634*c-11*c*c-1000*c**3
assert s.factor(delta-4*J/(5*(11*c-5)))==0
assert J.subs(c,s.Rational(33,40))==s.Rational(103619,400)
assert J.subs(c,s.Rational(9,10))==s.Rational(13769,100)
gates('7/10','-80')
print('PASS exact Lienard reduction, divergence identity, multiplier quartic and local quintic algebra')
print('PASS exact rational remote stable-focus gate at starting field')
print('No cycle count or interval-return certificate asserted')
