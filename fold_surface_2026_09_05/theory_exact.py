#!/usr/bin/env python3
"""Exact identities for the fold-component theory lane. Zero ODE calls."""
import json
from pathlib import Path
import sympy as s

c, K, x, y, u, v, X, Y, b, d = s.symbols('c K x y u v X Y b d')
J = 305+634*c-11*c*c-1000*c**3
e = 11*c-5
m = 210/e
sigma = 5*(2*c+1)/21
D = 1+m*sigma**2
A = (10*c*sigma*m+5*sigma**3*m*m-sigma**2*m-5*sigma*m+100*sigma+11)/(5*D**2)
B = sigma*(5*c-50*sigma**2-6*sigma-5)/(5*D**2)
C = (5*c*sigma**3*m*m+11*sigma**2*m+5*sigma*m-50*sigma+5)/(5*D**2)

# Rational interval isolation and sign checks for the chosen algebraic root.
center_interval=(s.Rational(96862063355349,10**14),s.Rational(96862063355350,10**14))
assert J.subs(c,center_interval[0])>0>J.subs(c,center_interval[1])
def imul(a,b):
    values=[aa*bb for aa in a for bb in b]
    return min(values),max(values)
def polynomial_interval(expr):
    value=(s.S.Zero,s.S.Zero)
    for coefficient in s.Poly(expr,c).all_coeffs():
        value=imul(value,center_interval)
        value=(value[0]+coefficient,value[1]+coefficient)
    return value
def rational_interval(expr):
    nn,dd=s.fraction(s.cancel(expr))
    ni,di=polynomial_interval(nn),polynomial_interval(dd)
    assert not di[0]<=0<=di[1]
    return imul(ni,(1/di[1],1/di[0]))
bi,di=rational_interval(B*m/A),rational_interval(C/A)
assert bi[0]>-s.Rational(3,10) and bi[1]<-s.Rational(29,100)
assert di[0]>s.Rational(101,100) and di[1]<s.Rational(51,50)
assert rational_interval(A)[0]>0

def at_center_zero(expr):
    numerator = s.fraction(s.cancel(expr))[0]
    assert s.rem(numerator, J, c) == 0

xx=(v-sigma*u)/D
yy=(u+m*sigma*v)/D
P=yy+xx**2+xx*yy
Q=-m*xx-10*xx**2+s.Rational(11,5)*xx*yy+c*yy**2
for polynomial in (Q-m*sigma*P-(-m*v+A*u*v), P+sigma*Q-(u+B*u*u+C*v*v)):
    for coef in s.Poly(s.cancel(polynomial),u,v).coeffs():
        at_center_zero(coef)
at_center_zero(A+2*C-s.Rational(21,5))

# Exact first integral in the fixed normalized center coordinates.
ss=1-X
a0=(1+b)/d
a1=2*(1+2*b)/(2*d-1)
a2=b/(d-1)
qform=Y**2+a0-a1*ss+a2*ss**2
P0=-ss*Y
Q0=X+b*X**2+d*Y**2
# H=(qform)*ss**(-2d)-qform(0,0). Clear the positive power.
assert s.factor((s.diff(qform,X)+2*d*qform/ss)*P0+s.diff(qform,Y)*Q0)==0

# Integrating-factor divergence and integration-by-parts identity.
k, ac, am, dc, dm = s.symbols('k ac am dc dm')
nu=2*d+1
P1=ac*dc*(X+k*Y)**2+am*dm*(k*X-Y)
Q1=k*P1
weighted_div=s.expand(ss*(s.diff(P1,X)+s.diff(Q1,Y))+nu*P1)
even=s.expand((weighted_div+weighted_div.subs(Y,-Y))/2)
target=ac*dc*(2*(1+k*k)*ss*X+nu*(X*X+k*k*Y*Y))+am*dm*k*nu*X
assert s.expand(even-target)==0
test=s.factor(s.diff(Y/ss,X)*P0+s.diff(Y/ss,Y)*Q0)
assert s.factor(test-(X+b*X*X+(d-1)*Y*Y)/ss)==0

# Exact multiplier-polynomial exclusion at c=1.
N1=(5*K-s.Rational(132,5))*u**3-s.Rational(48,5)*u*u+6*u+30
K0=s.Rational(6292,1125)
assert s.factor(N1.subs(K,K0)-s.Rational(2,225)*(11*u+15)*(4*u-15)**2)==0
assert s.diff(N1,K)==5*u**3

# General N polynomial in u=1+x and explicit sufficient no-cycle gate.
mm=5*(K+42)/e
W=(s.Rational(61,5)-c)*u**3+(mm-s.Rational(72,5)+3*c)*u*u+(s.Rational(11,5)-3*c)*u+c
dd=16-10*c
N=s.Poly(s.cancel((dd*u+(c+1)*(21+dd*(u-1)))*W-u*(21+dd*(u-1))*s.diff(W,u)),u)
n=[s.factor(N.coeff_monomial(u**j)) for j in range(5)]
T=c*(c+1)*(40*c-43)**2/(20*(2*c+1))
K2=-(8800*c**5-58120*c**4+92089*c**3+39218*c*c-53831*c-31720)/(500*(c-1)*(2*c+1)**2)
K3=-(2200*c**4-19095*c**3+40912*c*c-20175*c-3050)/(50*c*(5*c-8))
assert s.factor(n[2].subs(K,K2)-T)==0
assert s.factor(n[3].subs(K,K3))==0
assert s.factor(n[1]**2/(4*n[0])-T)==0

# A global two-parameter residual linearization in Lienard coordinates.
kappa=c/(2*c+1); gamma=(c+1)/(2*c+1)
Fpoly=-dd*u*u/(5*(1-c))+(21-2*dd)*u/(5*c)+(dd-21)/(5*(c+1))
ares=-J/(25*(c-1)*(2*c+1)**2)
bres=(5*K+(26-11*c)*ares)/e
assert s.factor(W+kappa*gamma*(dd*u+21-dd)*Fpoly/5-u*u*(ares*u+bres))==0

# c=8/5 endpoint: exact infinity exponents and separatrix jets.
z0,z,aa,bb,mm=s.symbols('z0 z aa bb mm')
pz=-10+s.Rational(6,5)*z+s.Rational(3,5)*z*z
az=5*(mm+z0*z0)/(11*(1+z0))
bz=(5*az*z0-8*az*az)/(16*(1+z0))
jet=z0+az*v+bz*v*v
Z=pz.subs(z,jet)-v*(mm+jet*jet)
V=-v*(1+jet)-jet*v*v
remainder=s.Poly(s.cancel(Z-s.diff(jet,v)*V),v)
for degree in (0,1,2):
    num=s.fraction(s.cancel(remainder.coeff_monomial(v**degree)))[0]
    assert s.rem(num,3*z0*z0+6*z0-50,z0)==0
disc=s.sqrt(40*c-s.Rational(964,25))
rho=(disc+2*c-s.Rational(16,5))/(disc-2*c+s.Rational(16,5))
assert s.simplify(rho.subs(c,s.Rational(8,5))-1)==0
assert s.simplify(s.diff(rho,c).subs(c,s.Rational(8,5))-10/s.sqrt(159))==0
tt=2*c-s.Rational(16,5)
nu_forward=2*(c-1)*disc/(disc-tt)
nu_backward=2*(c-1)*disc/(disc+tt)
nu_difference=10*disc*(s.Rational(8,5)-c)/(61-5*c)
assert s.simplify(nu_backward-nu_forward-nu_difference)==0
assert s.simplify(nu_forward/nu_backward-rho)==0
assert s.simplify(nu_forward.subs(c,s.Rational(8,5))-s.Rational(6,5))==0
assert s.simplify(nu_backward.subs(c,s.Rational(8,5))-s.Rational(6,5))==0
assert s.simplify(-s.diff(nu_difference,c).subs(c,s.Rational(8,5))-12/s.sqrt(159))==0

# No invariant conic with the required pair of infinity directions.
lc,nc,tc,acon,bcon,gcon=s.symbols('lc nc tc acon bcon gcon')
PP=y+x*x+x*y
QQ=-10*x*x+s.Rational(11,5)*x*y+s.Rational(8,5)*y*y-mm*x
conic=y*y+2*x*y-s.Rational(50,3)*x*x+acon*x+bcon*y+gcon
inv=s.Poly(s.expand(s.diff(conic,x)*PP+s.diff(conic,y)*QQ-(lc*x+nc*y+tc)*conic),x,y)
assert s.solve([coef for mon,coef in inv.terms() if sum(mon)==3],(lc,nc),dict=True)==[{lc:s.Rational(16,5),nc:s.Rational(16,5)}]
coeff=s.Poly(inv.as_expr().subs({lc:s.Rational(16,5),nc:s.Rational(16,5)}),x,y)
e1=coeff.coeff_monomial(x*x);e2=coeff.coeff_monomial(x*y);e3=coeff.coeff_monomial(y*y)
forced=s.solve([e1-e2,e3],(bcon,tc),dict=True)[0]
assert forced=={bcon:s.Rational(20,11),tc:-s.Rational(10,11)}
assert coeff.coeff_monomial(1).subs(forced)==s.Rational(10,11)*gcon
forced[gcon]=0;forced[acon]=-s.Rational(200,121);forced[mm]=-s.Rational(100,121)
assert coeff.coeff_monomial(x).subs(forced)==0
assert coeff.coeff_monomial(y).subs(forced)==0
assert e1.subs(forced)==-s.Rational(10180,363)

result={'status':'EXACT_IDENTITIES_PASS','orbit_evaluations':0,
        'checks':['reversible coordinates modulo J','A+2C=21/5',
                  'Darboux first integral','Melnikov even divergence',
                  'moment integration-by-parts identity','c=1 no-cycle factorization',
                  '1<c<8/5 sufficient no-cycle thresholds',
                  'global Lienard residual becomes (u-1)u^(-c)(a*u+b)',
                  'c=8/5 infinity separatrix jet through v^2',
                  'infinity product derivative 10/sqrt(159)',
                  'half-passage exponents 6/5 and difference slope 12/sqrt(159)',
                  'no invariant conic carrying both c=8/5 saddles'],
        'c1_no_cycle_K_threshold':str(K0),
        'algebraic_center_isolation':[str(z) for z in center_interval],
        'proved_normalized_bounds':{'b':'(-3/10,-29/100)','d':'(101/100,51/50)','A':'positive'},
        'K2':str(K2),'K3':str(K3),
        'scope':'No global at-most-two theorem or complete fold-component exclusion.'}
out=Path(__file__).with_name('theory_exact.json')
out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
