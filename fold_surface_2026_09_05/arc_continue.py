"""Pseudo-arclength continuation in (log r,c,K), including vertical K tangents.

Newton rank is checked at every corrected point. Regular-branch orientation
comes from the nullspace of the two exact return-equation derivatives.
No failed correction is promoted to a mathematical endpoint.
"""
import json,sys,math
import numpy as np
from budget import HERE,call
from continue_fold import paired_profile,topology

def eval_at(v,purpose):
    return call(dict(r=math.exp(v[0]),c=format(v[1],'.17g'),K=format(v[2],'.17g'),tol=2e-12),purpose,engine='log_variational.py')
def equations(a):
    return np.array([a['L'],a['L_z']]),np.array([[a['L_z'],a['L_c'],a['L_K']],[a['L_zz'],a['L_zc'],a['L_zK']]])
def tangent(J,old):
    _,s,vt=np.linalg.svd(J);t=vt[-1]
    if t@old<0:t=-t
    return t,s
def correct(pred,tangent0):
    v=pred.copy();history=[]
    for _ in range(9):
        if not (v[2]>0 and 5/11<v[1]<12.2 and -12<v[0]<95):return dict(status='CHART_GUARD',history=history)
        a=eval_at(v,'arclength_fold_corrector');history.append(a)
        if a['status']!='NUMERICAL_ONLY':return dict(status='UNRESOLVED_RETURN',history=history)
        f,J=equations(a);plane=float(tangent0@(v-pred))
        if max(abs(f))<2e-8 and abs(plane)<1e-8 and abs(a['first_derivative_discrepancy'])<1e-5:
            t,s=tangent(J,tangent0)
            return dict(status='ACCEPTED_NUMERICAL_FOLD',r=float(math.exp(v[0])),log_r=float(v[0]),
                        c=float(v[1]),K=float(v[2]),arc_tangent=t.tolist(),rank_singular_values=s.tolist(),
                        return_data=a,history=history)
        A=np.vstack([J,tangent0]);rhs=-np.r_[f,plane]
        try:step=np.linalg.solve(A,rhs)
        except np.linalg.LinAlgError:return dict(status='AUGMENTED_RANK_FAILURE',history=history)
        scale=min(1.,.7/max(abs(step[0]),1e-99),.03/max(abs(step[1]),1e-99),.6/max(abs(step[2]),1e-99))
        v+=scale*step
    return dict(status='CORRECTOR_ITERATION_LIMIT',history=history)

def run():
    point=json.loads((HERE/'events_increasing.json').read_text())[-1]['fold']
    v=np.array([point['log_r'],point['c'],point['K']]);
    t=np.r_[point['tangent'],1.];t/=np.linalg.norm(t)
    events=[];ds=.4
    for i in range(35):
        pred=v+ds*t;fold=correct(pred,t);event=dict(step=i,arc_step=ds,predictor=pred.tolist(),fold=fold)
        if fold['status']=='ACCEPTED_NUMERICAL_FOLD':
            v=np.array([fold['log_r'],fold['c'],fold['K']]);t=np.array(fold['arc_tangent'])
            event['topology']=topology(fold['c'],fold['K'])
            event['pair']=paired_profile(fold,full=True)
            if len(event['pair']['root_brackets'])>=3:
                events.append(event);(HERE/'K1_TRIGGER.json').write_text(json.dumps(event,indent=2));break
            ds=min(.7,ds*1.15)
        else:
            ds*=.5
        events.append(event);(HERE/'events_arclength.json').write_text(json.dumps(events,indent=2))
        if ds<.003:break
    (HERE/'events_arclength.json').write_text(json.dumps(events,indent=2))

if __name__=='__main__':run()
