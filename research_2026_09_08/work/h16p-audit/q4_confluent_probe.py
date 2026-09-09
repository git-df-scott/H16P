"""Confluent original-Q4 basis probe via linear PF reconstruction.
Tests determinant zeros that could permit a fourth-order compact zero.
Floating ODE only; not a global Wronskian sign theorem.
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import hyp2f1
P=Path(__file__).resolve().parent

def probe(a,rtol):
    y0=3*np.array([1326.,864.,-2431.,-102.])/1361360
    v0=-1.5*(1+a)*y0-np.array([0.,0.,1.,0.])/192
    initial=np.array([np.zeros(4),y0,v0,np.zeros(4)]).ravel() # H,Y,V,X
    def rhs(t,flat):
        H,Y,V,X=flat.reshape(4,4)
        F=hyp2f1(1/6,5/6,1,t);Fp=5/36*hyp2f1(7/6,11/6,2,t)
        M=1-6*(1-t)*Fp/F
        q=np.array([1.,t,-M,t*M-1])
        h_over_t2=np.array([1.,0.,-1/6,-1.])/2 if t==0 else H/(t*t)
        vp=(-h_over_t2/(1152*(1-t))+(1-a)*V/2-5*a*Y/36)/((1-a*t)*(1-t))
        return np.array([t*F*q,V,vp,Y/(1-a*t)**1.5]).ravel()
    points=[.001,.003,.01,.03,.1,.2,.4,.6,.8,.9,.95,.98,.99,.995,.999,.9999,.99999,.999999]
    sol=solve_ivp(rhs,[0,points[-1]],initial,method='DOP853',rtol=rtol,atol=rtol*1e-5,max_step=.01,dense_output=True)
    if not sol.success:return {'status':'failed','message':sol.message}
    rows=[]
    for t in points:
        H,Y,V,X=sol.sol(t).reshape(4,4)
        matrix=np.array([X,Y,V,H]);norms=np.linalg.norm(matrix,axis=1);scaled=matrix/norms[:,None]
        det=float(np.linalg.det(scaled));sv=np.linalg.svd(scaled,compute_uv=False)
        rows.append({'t':t,'normalized_state_determinant':det,'smallest_singular_value':float(sv[-1]),'states_H_Y_V_X':[H.tolist(),Y.tolist(),V.tolist(),X.tolist()]})
    return {'status':'passed','rtol':rtol,'nfev':sol.nfev,'rows':rows}

def main():
    out={'scope':'Nonvalidated PF Wronskian proxy det[X,Y,V,H]; W_X=f^3*c*det where f=(1-at)^(-3/2), c=-1/[1152t²(1-at)(1-t)²]. No original-basis global sign proof.','records':[]}
    for r in [.03125,.125,.5,1.,2.,8.,32.,128.]:
        a=r*r/(1+r*r);runs=[probe(a,tol) for tol in [2e-11,3e-13]]
        record={'r':r,'a':a,'runs':runs};out['records'].append(record)
        (P/'q4_confluent_probe.json').write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps({'r':r,'statuses':[run['status'] for run in runs],'determinants':[[row['normalized_state_determinant'] for row in run.get('rows',[])] for run in runs]}),flush=True)
if __name__=='__main__':main()
