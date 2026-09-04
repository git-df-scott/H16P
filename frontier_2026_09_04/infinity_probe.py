"""Nonrigorous two separatrix matching probe on the KKL Hopf section beta=0."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
def branch(c,alpha,sgn,Y=1e4,b=2.2):
 # Formal x(y)=-1-1/((c+1)y)-(b+2)/((c+1)(2c+1)y²)+O(y^-3).
 y0=sgn*Y; x0=-1-1/((c+1)*y0)-(b+2)/((c+1)*(2*c+1)*y0*y0)
 def rhs(y,z):
  x=z[0]; q=-10*x*x+b*x*y+c*y*y+alpha*x
  return [(y+x*x+x*y)/q]
 sol=solve_ivp(rhs,[y0,0],[x0],method='DOP853',rtol=2e-12,atol=2e-13)
 if not sol.success: raise RuntimeError('branch turns or fails before y=0')
 return float(sol.y[0,-1])
def splitting(c,a,Y=1e4):return branch(c,a,1,Y)-branch(c,a,-1,Y)
if __name__=='__main__':
 rows=[]
 for c in [.5,.7,.85,.95]:
  for a in [-200,-100,-77.7778,-72.7778,-40,-15]:
   try:r={'c':c,'alpha':a,'S':splitting(c,a),'upper':branch(c,a,1),'lower':branch(c,a,-1)}
   except Exception as e:r={'c':c,'alpha':a,'error':str(e)}
   rows.append(r);print(json.dumps(r),flush=True)
 (Path(__file__).resolve().parent/'data'/'infinity_probe.json').write_text(json.dumps({'status':'NONRIGOROUS; asymptotic truncation not bounded','rows':rows},indent=2)+'\n')
