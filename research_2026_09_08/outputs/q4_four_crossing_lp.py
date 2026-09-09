"""Finite-sample four-crossing feasibility over the full Q4 control space.
LP and frozen integrals are floating point; no continuum exclusion is implied.
"""
from pathlib import Path
from itertools import combinations
import json,hashlib,time
import numpy as np
from scipy.optimize import linprog
P=Path(__file__).resolve().parent

def solve_patterns(M,ts,npoints):
    row_scale=np.linalg.norm(M,axis=1);N=M/row_scale[:,None]
    _,sing,Vh=np.linalg.svd(N,full_matrices=False)
    transform=Vh.T/sing
    H=N@transform;H/=np.linalg.norm(H,axis=1)[:,None]
    sign=(-1.)**np.arange(npoints)
    best=None;failures=[];positive=[];count=0
    for inds in combinations(range(len(ts)),npoints):
        # Sign reversal is equivalent to replacing all controls by their negative.
        B=sign[:,None]*H[list(inds)]
        result=linprog([0,0,0,0,-1],A_ub=np.column_stack([-B,np.ones(npoints)]),b_ub=np.zeros(npoints),bounds=[(-1,1)]*4+[(0,None)],method='highs',options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9})
        count+=1
        if not result.success:
            failures.append({'indices':inds,'status':result.status,'message':result.message});continue
        gamma=float(result.x[-1]);d=transform@result.x[:4]
        if np.max(abs(d))>0:d/=np.max(abs(d))
        actual=M@d;record={'indices':inds,'fractions':[ts[i] for i in inds],'conditioned_margin':gamma,'direction':d.tolist(),'physical_values':actual.tolist(),'selected_signed_values':(sign*actual[list(inds)]).tolist(),'lp_residual_max':float(np.max(result.ineqlin.residual*-1))}
        if best is None or gamma>best['conditioned_margin']:best=record
        if gamma>1e-8:positive.append(record)
    return {'tested_patterns':count,'best':best,'positive_above_1e_8':positive,'failures':failures,'normalized_matrix_singular_values':sing.tolist()}

def main():
    source=P/'q4_modulus_boundary_constraints.json';ref=P/'q4_modulus_constraint_refine.json'
    records=json.loads(source.read_text())['records'];refined={r['modulus']:r for r in json.loads(ref.read_text())['records']}
    out={'scope':'Nonvalidated LP over discrete energies; tests for four sign changes, not a root bound or a five-cycle field','source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,ref]},'records':[]}
    for record in records:
        r=record['modulus'];run=refined.get(r,record)['runs'][-1];rows=run['rows'];assert all('M' in x for x in rows)
        M=np.array([x['M'] for x in rows]);ts=[x['fraction'] for x in rows]
        four=solve_patterns(M,ts,5)
        # Independent positive control: three alternating crossings are possible.
        inds=[2,5,9,12];three=solve_patterns(M[inds],np.array(ts)[inds].tolist(),4)
        entry={'modulus':r,'quadrature_order':run['order'],'four_crossings':four,'three_crossing_control':three};out['records'].append(entry)
        (P/'q4_four_crossing_lp.json').write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps({'modulus':r,'patterns':four['tested_patterns'],'best_margin':four['best']['conditioned_margin'],'positive':len(four['positive_above_1e_8']),'failures':len(four['failures']),'control_margin':three['best']['conditioned_margin']}),flush=True)
if __name__=='__main__':main()
