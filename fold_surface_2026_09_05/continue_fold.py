"""Predictor/corrector continuation of D=D_r=0, never a simple-cycle trace.

Fixed-K charts use analytic second flow variations. Near chart singularity
the corrector returns a chart-change event rather than claiming an endpoint.
Every accepted point receives a paired-side return profile; coverage beyond
the computed log-radius interval remains explicitly unresolved.
"""
import json,math,sys
from fractions import Fraction as F
import numpy as np
from budget import call,HERE,used

def evaluate(z,c,K,purpose,tol=2e-12):
    return call(dict(r=math.exp(z),c=format(c,'.17g'),K=format(K,'.17g'),beta='0',tol=tol),purpose,engine='log_variational.py')

def correct(K,z,c,purpose,maxiter=9):
    history=[]
    for _ in range(maxiter):
        a=evaluate(z,c,K,purpose)
        history.append(a)
        if a['status']!='NUMERICAL_ONLY':return dict(status='UNRESOLVED_RETURN',history=history)
        f=np.array([a['L'],a['L_z']]);J=np.array([[a['L_z'],a['L_c']],[a['L_zz'],a['L_zc']]])
        scale=max(abs(K),1e-12)
        if max(abs(f))/scale<2e-7 and abs(a['first_derivative_discrepancy'])<1e-6:
            tangent=np.linalg.solve(J,-np.array([a['L_K'],a['L_zK']]))
            return dict(status='ACCEPTED_NUMERICAL_FOLD',K=K,c=c,log_r=z,r=math.exp(z),
                        return_data=a,tangent=tangent.tolist(),jacobian=J.tolist(),history=history)
        try:step=np.linalg.solve(J,-f)
        except np.linalg.LinAlgError:return dict(status='FIXED_K_CHART_SINGULAR',history=history)
        factor=min(1.,.7/max(abs(step[0]),1e-99),.03/max(abs(step[1]),1e-99))
        z+=factor*step[0];c+=factor*step[1]
        if abs(z)>25 or c<=5/11 or c>12.2:return dict(status='NUMERICAL_DOMAIN_GUARD',history=history)
    return dict(status='CORRECTOR_ITERATION_LIMIT',history=history)

def topology(c,K):
    c,K=F(str(c)),F(str(K));alpha=-5*(K+42)/(11*c-5)
    A,B,C,D=c-F(61,5),alpha-F(111,5),2*alpha-10,alpha
    discr=18*A*B*C*D-4*B**3*D+B*B*C*C-4*A*C**3-27*A*A*D*D
    J=305+634*c-11*c*c-1000*c**3
    KH=-441*J/(125*(16-10*c)*(1+2*c)**2) if c!=F(8,5) else None
    return dict(c=str(c),K=str(K),alpha=str(alpha),cubic_discriminant=str(discr),
                one_remote_equilibrium_certified_by_discriminant=discr<0,
                J=str(J),K_H=str(KH) if KH else None,
                remote_stable_gate_in_inherited_box=(bool(K>KH) if KH is not None and F(1,2)<=c<=F(3,2) else None))

def root_refine(c,K,left,right,fa,fb,purpose):
    # Safeguarded log-coordinate secant; preserve numerical endpoint signs.
    records=[]
    for _ in range(9):
        z=(left*fb-right*fa)/(fb-fa)
        z=max(left+.1*(right-left),min(right-.1*(right-left),z))
        a=call(dict(r=math.exp(z),c=format(c,'.17g'),K=format(K,'.17g'),beta='0',log_radius_cap=100),purpose,engine='compact')
        records.append(a)
        if a['status']!='NUMERICAL_ONLY':break
        f=a['log_displacement']
        if abs(f)<1e-9 and abs(a['log_displacement_derivative'])>1e-6:
            return dict(status='NUMERICAL_SIMPLE_ROOT',root=a,log_bracket=[left,right],signs=[fa,fb],history=records)
        if fa*f<=0:right,fb=z,f
        else:left,fa=z,f
    return dict(status='NUMERICAL_SIGN_BRACKET',log_bracket=[left,right],signs=[fa,fb],history=records)

def paired_profile(point,full=True):
    a=point['return_data'];K=point['K'];c=point['c'];z0=point['log_r']
    # Move minimum downward; parameter step scales with nonzero K near center.
    shift=-math.copysign(min(abs(K)*.002,.01*abs(a['L_zz']/a['L_c'])),a['L_c'])
    cp=c+shift
    offsets=[-5.,-2.,-1.,-.35,0.,.35,1.,2.,4.,7.,11.,16.] if full else [-1.,0.,1.,5.,12.]
    profile=[]
    for dz in offsets:
        z=z0+dz
        if z>95:continue
        out=call(dict(r=math.exp(z),c=format(cp,'.17g'),K=format(K,'.17g'),beta='0',log_radius_cap=100),
                 'fold_pair_profile',engine='compact')
        profile.append(dict(log_r=z,result=out))
    brackets=[]
    for aa,bb in zip(profile,profile[1:]):
        x,y=aa['result'],bb['result']
        if x['status']!='NUMERICAL_ONLY' or y['status']!='NUMERICAL_ONLY':continue
        fa,fb=x['log_displacement'],y['log_displacement']
        if fa*fb<0:
            brackets.append(root_refine(cp,K,aa['log_r'],bb['log_r'],fa,fb,'pair_root_isolation'))
    stationary_brackets=[]
    for aa,bb in zip(profile,profile[1:]):
        x,y=aa['result'],bb['result']
        if x['status']=='NUMERICAL_ONLY' and y['status']=='NUMERICAL_ONLY':
            if x['log_displacement_derivative']*y['log_displacement_derivative']<0:
                stationary_brackets.append([aa['log_r'],bb['log_r']])
    return dict(c=cp,K=K,shift=shift,profile=profile,root_brackets=brackets,
                stationary_brackets=stationary_brackets,
                full_root_coverage_proved=False,unresolved='unsampled intervals and return-domain endpoints')

def run(direction,targets):
    seed=json.loads((HERE/'derivative_controls.json').read_text())['seed']
    z=math.log(seed['r']);c=float(F(seed['c']));K=float(F(seed['K']))
    point=correct(K,z,c,'seed_chart_corrector',3)
    events=[]
    for target in targets:
        if point['status']!='ACCEPTED_NUMERICAL_FOLD':break
        delta=target-point['K'];pred=np.array([point['log_r'],point['c']])+delta*np.array(point['tangent'])
        candidate=correct(target,*pred,'fold_'+direction)
        event=dict(direction=direction,target_K=target,fold=candidate)
        if candidate['status']=='ACCEPTED_NUMERICAL_FOLD':
            event['topology']=topology(candidate['c'],target)
            event['pair']=paired_profile(candidate,full=True)
            point=candidate
        events.append(event)
        (HERE/('events_'+direction+'.json')).write_text(json.dumps(events,indent=2))
        if len(event.get('pair',{}).get('root_brackets',[]))>=3:
            (HERE/'K1_TRIGGER.json').write_text(json.dumps(event,indent=2))
            print('K1 numerical sign trigger: stop broad exploration',flush=True);break
        if candidate['status']!='ACCEPTED_NUMERICAL_FOLD':break

if __name__=='__main__':
    direction=sys.argv[1]
    targets=([1/1024,1/2048,1/4096,1/8192] if direction=='decreasing' else
             [1/256,1/128,1/64,1/32,1/16,.1,.15,.2,.3,.45,.65,.9,1.2,1.6,2.,3.,4.,6.])
    run(direction,targets)
