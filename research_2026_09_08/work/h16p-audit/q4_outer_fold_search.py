"""Bounded numerical solution attempt for an outer double return-map zero.
Three interior equalities are fitted at every field evaluation. Center/large-radius
smallness is deflated by epsilon and a nonzero reference displacement.
"""
import json,time
import numpy as np
from scipy.optimize import root,least_squares
from q4_finite_continuation import P,difference,COUNTS
from q4_radial_variation import radial,COUNTS as RADIAL_COUNTS

def main(use_domain=False,outname="q4_outer_fold_search.json",negative=False):
    seeds={a['r']:a for a in json.loads((P/('q4_finite_negative.json' if negative else 'q4_finite_extended.json')).read_text())['records'] if a['epsilon']==(-.2 if negative else 4)}
    reference_eps=-.1 if negative else .1
    refs={a['r']:a for a in json.loads((P/('q4_finite_negative.json' if negative else 'q4_finite_continuation.json')).read_text())['records'] if a['epsilon']==reference_eps}
    out=dict(scope=__doc__,records=[]);start=time.perf_counter()
    for r,seed in seeds.items():
        anchors=seed['anchors'];outer=anchors[-1]
        for initial_factor in [1.3,2.]:
            fit_cache={seed['epsilon']:(np.array(seed['normalized_controls']),seed['scaled_fit_residual'])};fits=[];calls=[];domains={}
            def fitted(eps):
                if eps in fit_cache:return fit_cache[eps]
                warm=min(fit_cache,key=lambda e:abs(np.log(e/eps)));controls=fit_cache[warm][0];cached=None;history=[]
                def fun(c):
                    nonlocal cached
                    if cached is not None and np.array_equal(c,cached[0]):return cached[1:]
                    rows=[difference(r,eps,c,R,True) for R in anchors];history.append(dict(controls=c.tolist(),rows=rows))
                    if any('difference' not in q for q in rows):raise ArithmeticError('inner anchor passage unresolved')
                    res=np.array([q['difference'] for q in rows])/eps;jac=np.array([q['gradient'] for q in rows])/eps
                    cached=(c.copy(),res,jac);return res,jac
                try:
                    fit=root(lambda c:fun(c)[0],controls,jac=lambda c:fun(c)[1],options={'xtol':1e-9,'maxfev':60});res=fun(fit.x)[0]
                    if max(abs(res))>2e-10:raise ArithmeticError('inner fitting residual too large')
                    if np.linalg.norm(fit.x)<1e-5:raise ArithmeticError('near-center control gate')
                    tau,u,w=eps*fit.x
                    if negative and np.linalg.norm([tau,w,r+u])<1e-4:raise ArithmeticError('known reversible center distance gate')
                    fit_cache[eps]=(fit.x,res.tolist());return fit_cache[eps]
                finally:fits.append(dict(epsilon=eps,calls=history))
            def objective(z):
                eps=float((-1 if negative else 1)*np.exp(z[0]));c,res=fitted(eps)
                proposed_gap=float(outer*np.exp(z[1]));R=outer+proposed_gap
                if use_domain:
                    if eps not in domains:
                        from q4_finite_equilibria import equilibria
                        from q4_saddle_separatrices import field,shoot
                        eq=equilibria(r,eps,c);saddles=[q for q in eq['equilibria'] if q['determinant']<0];rs=None;shots=[]
                        if len(saddles)==1:
                            saddle=np.array(saddles[0]['point']);f,j=field(dict(r=r,epsilon=eps,normalized_controls=c));values,vectors=np.linalg.eig(j(saddle));assert max(abs(values.imag))<1e-10
                            idx=np.argmin(values.real);vec=vectors[:,idx].real;vec*=1 if vec[0]>=0 else -1
                            shots=[shoot(f,saddle,vec,-1,b,1e-9,3e-13) for b in [-1,1]]
                            passed=[q['radius'] for q in shots if q['status']=='negative_section']
                            if len(passed)!=1:raise ArithmeticError('saddle boundary unresolved')
                            rs=passed[0]
                            if rs<=outer*1.01:raise ArithmeticError('boundary too close to outer anchor')
                        domains[eps]=dict(stable_radius=rs,shots=shots)
                    rs=domains[eps]['stable_radius']
                    if rs is not None:R=outer+proposed_gap/(1+proposed_gap/(rs-outer))
                candidate=radial(r,eps,c,R);reference=radial(r,reference_eps,refs[r]['normalized_controls'],R)
                item=dict(epsilon=eps,radius=R,controls=c.tolist(),anchor_residual=res,candidate=candidate,reference=reference);calls.append(item)
                if 'radial_derivative' not in candidate or 'difference' not in reference:raise ArithmeticError('outer/reference passage unresolved')
                scale=eps*reference['difference']/reference_eps
                if negative:
                    tau,u,w=eps*c;rt,ru,rw=reference_eps*np.array(refs[r]['normalized_controls'])
                    distance=float(np.linalg.norm([tau,w,r+u]));reference_distance=float(np.linalg.norm([rt,rw,r+ru]))
                    item['reversible_center_distance']=distance
                    scale*=distance/reference_distance
                if abs(scale)<1e-11:raise ArithmeticError('reference scale too small')
                value=np.array([candidate['difference'],R*candidate['radial_derivative']])/scale
                item['scaled_fold_residual']=value.tolist();return value
            result=dict(r=r,initial_factor=initial_factor,fit_calls=fits,outer_calls=calls,domain_records=domains,use_domain=use_domain,negative=negative)
            try:
                sol=least_squares(objective,np.log([.2 if negative else 4.,initial_factor-1]),bounds=(np.log([.01 if negative else .1,.1]),np.log([.8 if negative else 4. if r==.5 and use_domain else 8. if use_domain else 12.,9.])),diff_step=1e-3,max_nfev=25,xtol=1e-8,ftol=1e-8,gtol=1e-8)
                final=objective(sol.x);result.update(status='terminated',solver_success=bool(sol.success),message=sol.message,solution_coordinates=sol.x.tolist(),scaled_fold_residual=final.tolist(),candidate_gate=bool(max(abs(final))<1e-6))
            except Exception as e:result.update(status='failed_gate',reason=repr(e))
            out['records'].append(result);out['counts']=COUNTS.copy();out['radial_counts']=RADIAL_COUNTS.copy();out['wall_seconds']=time.perf_counter()-start
            print(json.dumps({k:v for k,v in result.items() if k not in ['fit_calls','outer_calls','domain_records']}|{'outer_evaluations':len(calls)}),flush=True)
            (P/outname).write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
