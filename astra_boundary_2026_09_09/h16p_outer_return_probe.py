import numpy as np,json,time,sys
from scipy.integrate import solve_ivp
from pathlib import Path
OUT=Path(__file__).resolve().parent
def half(s,A,direction,delta=.02,side=1,rtol=3e-12,max_step=.05,method='DOP853',m=8/9,c=159/80):
 al=A*delta**2;e0=-delta**3/3;e1=-c*delta**3;ga=m*delta**3/3
 def rhs(t,z):
  R,U=z
  if R>100 or abs(U)>1e30:raise ValueError('trial-stage range failure; no orbit-exit conclusion')
  r=np.exp(R)
  return direction*np.array([side*(2*U-e0*r*r-ga*U*U),side/3+2*r/3-side*5*r*r/12+side*al*U*U+e1*r*U-side*e0*U*r*r-side*ga*U**3])
 z0=[-s,0.];dep=solve_ivp(rhs,[0,1e-6],z0,rtol=rtol,atol=rtol/100,method=method)
 def event(t,z):return z[1]
 event.terminal=True;event.direction=-direction*side
 sol=solve_ivp(rhs,[1e-6,2000],dep.y[:,-1],events=event,rtol=rtol,atol=rtol/100,method=method,max_step=max_step)
 if not sol.success or not len(sol.t_events[0]):raise ValueError('unresolved: '+sol.message)
 return dict(R=float(sol.y_events[0][0][0]),t=float(sol.t_events[0][0]),steps=len(sol.t),max_abs_U=float(np.max(np.abs(sol.y[1]))))
def run():
 start=time.monotonic();rows=[]
 for A in [-2.,-1.,-.65,-.60,-.585,-.5824,-.58,-.55,-.5,0.,1.]:
  for s in [4.04,5.,8.,12.,20.,40.,80.,160.,320.,640.,1280.,2560.,5120.]:
   row=dict(A=A,s=s)
   try:
    f=half(s,A,1);b=half(s,A,-1);row.update(D=f['R']-b['R'],forward=f,backward=b)
   except (ValueError,RuntimeError) as e:row['failure']=str(e)
   rows.append(row)
   if time.monotonic()-start>45:break
  print(json.dumps(dict(A=A,values=[(r['s'],r.get('D',r.get('failure'))) for r in rows if r['A']==A])),flush=True)
  (OUT/'h16p_outer_probe_replay.json').write_text(json.dumps(dict(scope='NUM only; exp(R) can underflow for R<-745; no interval or absence claim',delta=.02,rtol=3e-12,wall_seconds=time.monotonic()-start,rows=rows),indent=2))
  if time.monotonic()-start>45:break
if __name__=='__main__':run()
