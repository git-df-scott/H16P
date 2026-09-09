"""Continue three finite-cycle anchors by a common shift in base energy.
No validated zeros, manifold enclosures, or global cycle-count assertion.
"""
import json,time
import numpy as np
from scipy.optimize import root
from q4_finite_continuation import P,anchor,difference,COUNTS
from q4_finite_equilibria import equilibria
from q4_saddle_separatrices import field,shoot

def main(amounts=(.05,.1,.15,.2),outname="q4_anchor_placement.json",warm=None):
    seeds={q['r']:q for q in json.loads((P/'q4_finite_extended.json').read_text())['records'] if q['epsilon']==4}
    out={'scope':__doc__,'records':[],'failures':[]};start=time.perf_counter()
    for r,seed in seeds.items():
        eps=seed['epsilon']
        for direction in [-1,1]:
            controls=np.array(warm[(r,direction)] if warm is not None else seed['normalized_controls'])
            for amount in amounts:
                fractions=np.array([.25,.5,.75])+direction*amount
                anchors=[anchor(r,t) for t in fractions];calls=[];cache=None
                def evaluate(c):
                    nonlocal cache
                    if cache is not None and np.array_equal(c,cache[0]):return cache[1:]
                    rows=[difference(r,eps,c,R,True) for R in anchors]
                    calls.append(dict(controls=c.tolist(),rows=rows))
                    if any('difference' not in q for q in rows):raise ArithmeticError('anchor return unresolved')
                    residual=np.array([q['difference'] for q in rows])/eps
                    jac=np.array([q['gradient'] for q in rows])/eps
                    cache=(c.copy(),residual,jac);return residual,jac
                try:
                    fit=root(lambda c:evaluate(c)[0],controls,jac=lambda c:evaluate(c)[1],options={'xtol':2e-8,'maxfev':70})
                    residual=evaluate(fit.x)[0]
                    if max(abs(residual))>2e-9:raise ArithmeticError('fit residual exceeds threshold')
                    controls=fit.x
                    item=dict(r=r,epsilon=eps,shift=direction*amount,energy_fractions=fractions.tolist(),anchors=anchors,normalized_controls=controls.tolist(),solver_success=bool(fit.success),solver_message=fit.message,scaled_fit_residual=residual.tolist(),fit_calls=calls)
                    eq=equilibria(r,eps,controls);item['geometry']=eq
                    saddles=[q for q in eq['equilibria'] if q['determinant']<0]
                    rs=None
                    if len(saddles)==1:
                        s=np.array(saddles[0]['point']);f,j=field(item);values,vectors=np.linalg.eig(j(s));assert max(abs(values.imag))<1e-10
                        idx=np.argmin(values.real);v=vectors[:,idx].real;v*=1 if v[0]>=0 else -1
                        shots=[shoot(f,s,v,-1,b,1e-9,3e-13) for b in [-1,1]]
                        item['stable_shots']=shots
                        radii=[q['radius'] for q in shots if q['status']=='negative_section']
                        if len(radii)==1:rs=radii[0]
                    top=rs*(1-1e-7) if rs is not None else 1e5
                    radii=sorted(set(np.geomspace(max(1e-3,min(anchors)*.15),top,43).tolist()+[q*f for q in anchors for f in [.9,1.1] if q*f<top]+([rs-(rs-anchors[-1])*s for s in [1e-2,1e-4,1e-6]] if rs is not None and rs>anchors[-1] else [])))
                    rows=[difference(r,eps,controls,R) for R in radii]
                    brackets=[[a['radius'],b['radius']] for a,b in zip(rows,rows[1:]) if 'difference' in a and 'difference' in b and a['difference']*b['difference']<0 and min(abs(a['difference']),abs(b['difference']))>1e-11]
                    item.update(profile=rows,brackets=brackets,unresolved=sum('difference' not in a for a in rows),stable_radius=rs)
                    out['records'].append(item)
                    print(json.dumps(dict(r=r,shift=direction*amount,controls=controls.tolist(),brackets=brackets,unresolved=item['unresolved'],stable_radius=rs)),flush=True)
                except Exception as e:
                    out['failures'].append(dict(r=r,shift=direction*amount,reason=repr(e),calls=calls));print(json.dumps(dict(r=r,shift=direction*amount,failure=repr(e))),flush=True)
                    break
                finally:
                    out['counts']=COUNTS.copy();out['wall_seconds']=time.perf_counter()-start
                    (P/outname).write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
