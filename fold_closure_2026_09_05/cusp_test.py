"""Augmented degeneracy test from known finite fold, not a generic scan.
K=m(Bc-1)-10(B+2) is held positive; B is the xy coefficient in Q.
F=G=G_z=0 would be a triple-return-root candidate, requiring verification.
"""
import numpy as np,json
from pathlib import Path
from closure_budget import call,H
K=1/512
# Horizontal amplitude transported from the center limit plus known tangent.
x=np.array([np.log(6.7594),.9688884793906646,2.2]);records=[]
def ev(x,why):
 return call(dict(r=str(np.exp(x[0])),c=str(x[1]),K=str(K),B=str(x[2]),tol='2e-18'),why)
def val(a):return np.array([float(a[k]) for k in ('F','G','G_z')])/K
for it in range(4):
 a=ev(x,'fixed-focus augmented cusp residual')
 if a['status']!='NUMERICAL_TWO_HALF_PASSAGES':records.append(dict(status='UNRESOLVED',x=x.tolist(),result=a));break
 v=val(a)
 if max(abs(v))<1e-7:
  records.append(dict(status='TRIPLE_ROOT_NUMERICAL_SIGNAL',x=x.tolist(),result=a));break
 J=np.empty((3,3));steps=[1e-3,2e-6,1e-3];pert=[]
 for j,h in enumerate(steps):
  y=x.copy();y[j]+=h;b=ev(y,'finite-difference augmented cusp Jacobian');pert.append(b)
  if b['status']!='NUMERICAL_TWO_HALF_PASSAGES':break
  J[:,j]=(val(b)-v)/h
 if len(pert)!=3 or any(b['status']!='NUMERICAL_TWO_HALF_PASSAGES' for b in pert):records.append(dict(status='UNRESOLVED_JACOBIAN',result=a,perturbations=pert));break
 # First two analytic columns improve the rank test.
 J[0,:2]=[float(a['F_z'])/K,float(a['F_c'])/K];J[1,:2]=[float(a['G_z'])/K,float(a['G_c'])/K]
 dx=np.linalg.solve(J,-v);scale=min(1,.5/max(abs(dx[0]),1e-30),.05/max(abs(dx[1]),1e-30),.75/max(abs(dx[2]),1e-30));step=dx*scale
 records.append(dict(status='AUGMENTED_STEP',x=x.tolist(),residual=v.tolist(),jacobian=J.tolist(),newton_step=dx.tolist(),damping=scale,result=a,perturbations=pert));x+=step
(H/'cusp_attempt.json').write_text(json.dumps(dict(K=K,records=records,next_predictor=x.tolist(),component_excluded=False),indent=2))
