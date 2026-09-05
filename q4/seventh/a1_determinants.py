"""Exact limiting reconstruction and boundary determinant factorization."""
import sympy as S,json,resource,time
from pathlib import Path
resource.setrlimit(resource.RLIMIT_CPU,(60,60))
t=S.symbols('t');F,G,M=S.symbols('F G M'); c=S.Rational(5,36)
der=lambda z:S.diff(z,t)+S.diff(z,F)*G+S.diff(z,G)*((2*t-1)*G+c*F)/(t*(1-t))
j0=S.Rational(36,5)*t*(1-t)*G
j1=(t*t*(1-t)*G-t*(1-t)*F+j0)/(2+c)
j2=(t**3*(1-t)*G-2*t*t*(1-t)*F+4*j1)/(6+c)
K=list(map(S.factor,[j1,j2,6*j0-11*j1-6*t*(1-t)*F,12*j1-17*j2-6*t*t*(1-t)*F]))
Y=[9*(410*F*t+385*F+11448*G*t*t-11653*G*t)/1185800,
   9*(24436*F*t+11088*F+420239*G*t*t-426528*G*t)/52412360,
   3*(60*F*t+385*F+6408*G*t*t-6438*G*t)/215600,
   3*(7213*F*t+5544*F+176586*G*t*t-177228*G*t)/6166160]
Y=list(map(S.factor,Y))
center=[S.Rational(9,3080),S.Rational(162,85085),S.Rational(3,560),S.Rational(27,10010)]
for j in range(4):
 assert S.factor((1-t)**2*der(der(Y[j]))+c*Y[j]+K[j]/(1152*t*t*(1-t)))==0
 assert Y[j].subs({t:0,F:1,G:c})==center[j]
 yp=S.factor(der(Y[j]));lim=S.limit(yp.subs(G,c+S.Rational(385,2592)*t).subs(F,1+c*t),t,0)
 assert lim==-3*center[j]+(S.Rational(1,192) if j==2 else 0),(j,lim)
start=time.process_time();rows={};dets={}
for name,row in [('confluent',[S.factor(der(x)) for x in K]),('endpoint',list(map(S.Integer,[9061,6289,2431,1819])))]:
 mat=S.Matrix([K,row,Y,list(map(lambda z:S.factor(der(z)),Y))])
 det=S.factor(mat.det(method='domain-ge'))
 mm=S.factor(det.subs(G,(1-M)*F/(6*(1-t))))
 print(name,'F,G:',det,flush=True);print(name,'M:',mm,flush=True)
 rows[name]={'determinant_F_G':str(det),'determinant_M':str(mm)}
 dets[name]=mm
Q=S.factor(-dets['confluent']*2620618000/(81*F**4*t**3))
z,v=S.symbols('z v',positive=True)
den=1296*(v+1)**4*(z+1)**3*(z+6)**4
subQ=S.factor(Q.subs({t:z/(1+z),M:(S.Rational(1,6)+v*(1+z)/(6+z))/(1+v)}))
P=S.Poly(S.cancel(subQ*den/(125*z)),v,z)
assert all(co>0 for co in P.coeffs())
assert len(P.terms())==29
Qe=t*(M**2+4*M)-3*M-2
assert S.factor(dets['endpoint']-S.Rational(27,1600)*F**3*t**2*(M-1)*(t-1)*Qe)==0
assert S.factor(Qe-((M+2)*(M-1)-(1-t)*(M**2+4*M)))==0
u=S.symbols('u')
assert S.factor(S.diff((1-u)/(1-t*u),u,2)+2*t*(1-t)/(1-t*u)**3)==0
print('Confluent positivity polynomial:',P.as_expr(),flush=True)
Path(__file__).with_suffix('.json').write_text(json.dumps(dict(
 status='EXACT_A1_FACE_SIGNS_RELATIVE_TO_STIELTJES_REPRESENTATION',
 scope='a=1 only; no conclusion for all finite a follows by continuity',
 Y=list(map(str,Y)),determinants=rows,confluent_Q=str(Q),
 positivity_polynomial=str(P.as_expr()),positive_coefficient_count=len(P.terms()),
 cpu_seconds=time.process_time()-start),indent=2)+'\n')
