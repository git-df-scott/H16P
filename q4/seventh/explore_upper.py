"""Bounded diagnostics for a first-Duhamel upper bound, not a certificate."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import sys,json,resource,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'sixth'))
from boundary_diagnostic import anchored_pair,center,mp,np,solve_ivp,hyp2f1

def run(a,r,pair):
 co=np.array(pair,dtype=float); yc=np.array([float(center(c)) for c in pair]);D0=-yc[0]/192
 # H_E,H_V,Y_E,Y_V,W_E,W_V,e,z,W_e,W_z,D,Dlinear
 state=np.r_[np.zeros(2),yc,-1.5*(1+a)*yc+co[:,2]/192,1,0,-1.5*(1+a),1,D0,D0]
 def rhs(x,S):
  d=np.exp(-x);t=-np.expm1(-x);b=1-a*t
  H,Y,W=S[:2],S[2:4],S[4:6];hom,homw=S[6:8],S[8:10]
  F=hyp2f1(1/6,5/6,1,t);Fp=5/36*hyp2f1(7/6,11/6,2,t);M=1-6*d*Fp/F
  Hd=d*t*F*(co[:,0]+co[:,1]*t+(co[:,2]+co[:,3]*t)*M)
  hr=(co[:,0]+co[:,2]/6)/2 if x==0 else H/t**2
  Wd=-W+(1-a)*W/(2*b)-5*a*d*Y/(36*b)-hr/(1152*b)
  homwd=-homw+(1-a)*homw/(2*b)-5*a*d*hom/(36*b)
  fac=1/(1152*b**1.5*np.sqrt(d))
  Dd=(hr[1]*Y[0]-hr[0]*Y[1])*fac
  Ld=((hr[1]*yc[0]-hr[0]*yc[1])*hom[0]-hr[0]*hom[1]/192)*fac
  return np.r_[Hd,W,Wd,homw,homwd,Dd,Ld]
 sol=solve_ivp(rhs,(0,-np.log1p(-float(r))),state,rtol=3e-12,atol=2e-15,method='DOP853',max_step=.08)
 assert sol.success
 return dict(D0=D0,D=float(sol.y[-2,-1]),upper=float(sol.y[-1,-1]))

def main():
 resource.setrlimit(resource.RLIMIT_CPU,(45,45));os.nice(10);mp.mp.dps=55
 start=time.process_time();rows=[]
 for rt in ('.825','.99','.99999'):
  r=mp.mpf(rt)
  for face in ('confluent','endpoint'):
   pair=anchored_pair(r,face)
   for a in (.875,.99,1.):
    row=dict(r=rt,face=face,a=a,**run(a,r,pair));rows.append(row);print(row,flush=True)
 Path(__file__).with_suffix('.json').write_text(json.dumps(dict(status='NUMERICAL_ONLY',rows=rows,cpu_seconds=time.process_time()-start),indent=2)+'\n')
if __name__=='__main__':main()
