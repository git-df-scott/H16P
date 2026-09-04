"""Bounded return-map audit on beta=0; no parameter optimization or proof."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
def ret(x0,c,alpha,beta=0,b=2.2):
 def f(t,z):
  x,y,w=z
  return [y+x*x+x*y,-10*x*x+b*x*y+c*y*y+alpha*x+beta*y,(2+b)*x+(1+2*c)*y+beta]
 direction=np.sign(f(0,[x0,0,0])[1])
 def e1(t,z):return z[1]
 e1.direction=-direction;e1.terminal=True
 a=solve_ivp(f,[0,50],[x0,0,0],method='DOP853',events=e1,rtol=3e-12,atol=3e-13,max_step=.02)
 if not len(a.t_events[0]):raise RuntimeError('no return')
 def e2(t,z):return z[1]
 e2.direction=direction;e2.terminal=True
 d=solve_ivp(f,[a.t[-1],50],a.y[:,-1],method='DOP853',events=e2,rtol=3e-12,atol=3e-13,max_step=.02)
 if not len(d.t_events[0]):raise RuntimeError('no return')
 return float(d.y[0,-1]-x0),float(np.exp(d.y[2,-1]))
if __name__=='__main__':
 rows=[]
 for c in [.7,1.,1.3]:
  for alpha in [-120,-80,-72.7778,-50]:
   row={'c':c,'alpha':alpha,'K':-alpha*(2.2*c-1)-42,'samples':[],'roots':[]}
   samples=[.05,.2,.5,1.,2.,5.,10.,20.,50.,100.]
   for r in samples:
    try:row['samples'].append([r,*ret(r,c,alpha)])
    except Exception as e:row['samples'].append([r,str(e)])
   for a,b in zip(row['samples'],row['samples'][1:]):
    if len(a)==3 and len(b)==3 and a[1]*b[1]<0:
     try:
      r=brentq(lambda r:ret(r,c,alpha)[0],a[0],b[0],xtol=1e-10); row['roots'].append([r,*ret(r,c,alpha)])
     except:pass
   rows.append(row);print(json.dumps(row),flush=True)
 (Path(__file__).resolve().parent/'data'/'hopf_probe.json').write_text(json.dumps({'status':'NONRIGOROUS bounded 12-point probe; nonreturns and unsampled intervals unresolved','rows':rows},indent=2)+'\n')
