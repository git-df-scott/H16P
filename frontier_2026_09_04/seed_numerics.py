"""Independent numerical seed audit; floating-point evidence, not proof."""
import json, sys
from pathlib import Path
import numpy as np
import sympy as s
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
x,y=s.symbols('x y'); R=s.Rational
families={
 'gt_remote':(-s.Rational(1,10**200)*x-y-10*x*x+(5-s.Rational(1,10**13))*x*y+y*y,x+x*x+(-25-8*s.Rational(1,10**52)+9*s.Rational(1,10**13))*x*y,0,[(.04,.05)]),
 'shi_visual':(-R(2,10**8)*x-y-10*x*x+R(49,10)*x*y+y*y,x+x*x-R(30135,1250)*x*y,0,[(.002,.005),(.005,.012),(.026,.034),(.043,.10)]),
 'chen_visual':(-R(1,50000)*x-y-3*x*x+R(99,100)*x*y+y*y,x+R(2,9)*x*x-3*x*y,0,[(.02,.1),(.1,.14),(.176,.210),(.36,.49)]),
 'yu_zeng':(y*(1-R(30,7)*x)+x/R(20000),-x+x*x-R(671,210)*y*y-R(500001,10**10)*y+R(49182857,96810000000)*x*y,1,[(-.06,-.002),(-.25,-.06),(-.5,-.25),(2,400)]),
 'kkl':(y+x*x+x*y,-10*x*x+R(11,5)*x*y+R(7,10)*y*y-R(363889,5000)*x+R(3,2000)*y,1,[(.428,.785),(1.438,2.637),(8.858,16.238),(-3945,-2893)])}

def run(name,tol=1e-12):
 P,Q,axis,brackets=families[name]
 f=s.lambdify((x,y),[P,Q],'numpy'); div=s.lambdify((x,y),s.diff(P,x)+s.diff(Q,y),'numpy')
 def field(t,z): return [*f(z[0],z[1]),div(z[0],z[1])]
 def initial(r): return np.array([0.,r,0.] if axis==0 else [r,0.,0.])
 def ret(r):
  z=initial(r); direction=np.sign(field(0,z)[axis])
  def opposite(t,z): return z[axis]
  opposite.direction=-direction; opposite.terminal=True
  a=solve_ivp(field,[0,100],z,method='DOP853',events=opposite,rtol=tol,atol=tol*.01,max_step=.03)
  if not len(a.t_events[0]): raise RuntimeError('no opposite crossing')
  def desired(t,z): return z[axis]
  desired.direction=direction; desired.terminal=True
  b=solve_ivp(field,[a.t[-1],100],a.y[:,-1],method='DOP853',events=desired,rtol=tol,atol=tol*.01,max_step=.03)
  if not len(b.t_events[0]): raise RuntimeError('no return')
  return b.y[1-axis,-1]-r,b.t[-1],np.exp(b.y[2,-1]), np.max(np.abs(np.concatenate([a.y[:2],b.y[:2]],axis=1)),axis=1).tolist()
 data={'name':name,'evidence':'NONRIGOROUS','tolerance':tol,'P':str(P),'Q':str(Q),'cycles':[],'scan':[]}
 # A small scan resolves section positions, without counting discontinuous sign changes as roots.
 samples= sorted(set([v for br in brackets for v in br]+list(np.geomspace(.0005,.99,45)))) if axis==0 else sorted(set(v for br in brackets for v in br))
 for r in samples:
  try: val=ret(r); data['scan'].append([r,*val[:3]])
  except Exception as e: data['scan'].append([r,str(e)])
 for left,right in brackets:
  try:
   fl,fr=ret(left)[0],ret(right)[0]
   root=brentq(lambda r:ret(r)[0],left,right,xtol=5e-13,rtol=3e-14)
   res,T,mult,bounds=ret(root)
   data['cycles'].append({'bracket':[left,right],'endpoint_displacements':[fl,fr],'section_coordinate':root,'residual':res,'period':T,'multiplier':mult,'max_abs_xy':bounds})
  except Exception as e: data['cycles'].append({'bracket':[left,right],'error':str(e)})
 # Exact eliminant and all numerical real equilibria.
 elim=s.factor(s.resultant(P,Q,y)); data['equilibrium_eliminant_x']=str(elim)
 sols=[]
 for xx in s.polys.polytools.intervals(elim,eps=R(1,10**25)):
  lo,hi=xx[0]; xv=float((lo+hi)/2)
  for yv in s.nroots(P.subs(x,xv),maxsteps=200):
   if abs(s.im(yv))<1e-10 and abs(float(Q.subs({x:xv,y:s.re(yv)})))<1e-5:
    yy=float(s.re(yv)); J=np.array(s.Matrix([P,Q]).jacobian([x,y]).subs({x:xv,y:yy}),dtype=float)
    sols.append({'x':xv,'y':yy,'trace':float(np.trace(J)),'determinant':float(np.linalg.det(J)),'eigenvalues':[str(z) for z in np.linalg.eigvals(J)]})
 data['equilibria']=sols
 v=s.symbols('v'); P2=s.Poly(P,x,y).homogeneous_component(2) if False else sum(c*x**i*y**j for (i,j),c in s.Poly(P,x,y).terms() if i+j==2)
 Q2=sum(c*x**i*y**j for (i,j),c in s.Poly(Q,x,y).terms() if i+j==2)
 angular=s.factor(Q2.subs({x:1,y:v})-v*P2.subs({x:1,y:v})); data['infinity_angular_polynomial']=str(angular)
 data['infinity']=[]
 if P2.subs({x:0,y:1})==0:
  u=s.symbols('u'); vertical=s.expand(P2.subs({x:u,y:1})-u*Q2.subs({x:u,y:1})); data['infinity'].append({'direction':'vertical','chart_eigenvalues':[float(s.diff(vertical,u).subs(u,0)),float(-Q2.subs({x:0,y:1}))]})
 for vv in s.nroots(angular,maxsteps=200):
  if abs(s.im(vv))<1e-12:
   vr=float(s.re(vv)); data['infinity'].append({'slope':vr,'chart_eigenvalues':[float(s.diff(angular,v).subs(v,vr)),float(-P2.subs({x:1,y:vr}))]})
 (Path(__file__).resolve().parent/'data'/(name+'.json')).write_text(json.dumps(data,indent=2)+'\n'); print(json.dumps(data),flush=True)
if __name__=='__main__':run(sys.argv[1],float(sys.argv[2]) if len(sys.argv)>2 else 1e-12)
