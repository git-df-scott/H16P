"""Exact a-connection for the remaining function; no sign theorem."""
import sympy as S,resource,time,json
from pathlib import Path
resource.setrlimit(resource.RLIMIT_CPU,(60,60))
a,t,F,G,J,Jp=S.symbols('a t F G J Jp'); c=S.Rational(5,36)
j0=S.Rational(36,5)*t*(1-t)*G
K=(t*t*(1-t)*G-t*(1-t)*F+j0)/(2+c)
b=1-a*t;d=1-t;src=-K/(1152*t*t*d)
A=(1-a)/(2*b*d);B=-c*a/(b*d);srcn=src/(b*d)
der=lambda z:S.diff(z,t)+S.diff(z,F)*G+S.diff(z,G)*((2*t-1)*G+c*F)/(t*d)+S.diff(z,J)*Jp+S.diff(z,Jp)*(A*Jp+B*J+srcn)
h=-d/(a*(1-a))
deg=5;zs=S.symbols('z:'+str(2*deg+1))
U=(sum(zs[i]*t**i for i in range(deg))*F+sum(zs[deg+i]*t**i for i in range(deg))*G)/(t*b)
ans=h*Jp+zs[-1]*J+U
assert S.factor(h*S.diff(A,t)+S.diff(h,t)*A+S.diff(h,t,2)-S.diff(A,a))==0
assert S.factor(h*S.diff(B,t)+2*S.diff(h,t)*B-S.diff(B,a))==0
L=S.cancel(b*d*der(der(U))-(1-a)*der(U)/2+c*a*U)
expected=S.cancel(b*d*(S.diff(srcn,a)-h*der(srcn)-(2*S.diff(h,t)+zs[-1])*srcn))
res=S.Poly(S.together(L-expected).as_numer_denom()[0],t,F,G)
# Match the actual center data of J=Y_K0.
ell=S.Rational(9,3080);y1=-S.Rational(3,2)*(1+a)*ell
y2=S.factor(((1-a)*y1/2-c*a*ell-S.Rational(1,2304))/2)
series={F:1+c*t+S.Rational(385,5184)*t*t,G:c+S.Rational(385,2592)*t,J:ell+y1*t+y2*t*t,Jp:y1+2*y2*t}
cent=S.Poly(S.series(S.cancel(t*ans.subs(series)),t,0,3).removeO(),t)
centm=cent.coeff_monomial(1)
cent0=cent.coeff_monomial(t)
cent1=cent.coeff_monomial(t*t)
sol=S.linsolve(res.coeffs()+[centm,cent0,cent1+S.Rational(3,2)*ell],zs)
assert sol!=S.EmptySet
assert len(sol)==1
if sol!=S.EmptySet:
 result=S.factor(ans.subs(dict(zip(zs,next(iter(sol))))))
 cc=S.factor(S.diff(result,J))
 uu=S.factor(result-cc*J-h*Jp)
 assert S.factor(S.diff(result,Jp)-h)==0
 substitutions=dict(zip(zs,next(iter(sol))))
 assert all(S.factor(e.subs(substitutions))==0 for e in res.coeffs())
 u0=uu.subs({t:0,F:1,G:c})
 u1=(S.diff(uu,t)+S.diff(uu,F)*G+S.diff(uu,G)*S.Rational(385,2592)).subs({t:0,F:1,G:c})
 assert S.factor(cc*ell+h.subs(t,0)*y1+u0)==0
 assert S.factor(cc*y1+S.diff(h,t)*y1+h.subs(t,0)*2*y2+u1+S.Rational(3,2)*ell)==0
 print('Exact a-connection, compatibility, and center data: PASS',flush=True)
 Path(__file__).with_suffix('.json').write_text(json.dumps({'status':'EXACT_IDENTITY_ONLY','connection':str(result),'c':str(cc),'h':str(h),'U':str(uu)},indent=2)+'\n')
