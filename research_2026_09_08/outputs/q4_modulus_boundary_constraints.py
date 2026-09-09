"""Bounded modulus scan of the doubly constrained first variation."""
from pathlib import Path
import json
import numpy as np
import sympy as S
from q4_modulus_probe import compact,boundary
ROOT=Path(__file__).resolve().parent

def exact_exponent():
    r=S.symbols('r',positive=True);z,u,v,w=S.symbols('z u v w')
    p=-(2+r*r)+(2*r+2*u+w)*z+z*z;q=r+u+(-1-3*r*r+v)*z-(r+u)*z*z;R=S.expand(q-z*p);base={u:0,v:0,w:0}
    values=[]
    for control in [u,v,w]:
        terms=[]
        for root in [-r-S.sqrt(1+r*r),-r+S.sqrt(1+r*r)]:
            assert S.simplify(R.subs(base).subs(z,root))==0
            a=S.simplify(S.diff(R,z).subs(base).subs(z,root));b=S.simplify(p.subs(base).subs(z,root));dz=-S.diff(R,control).subs(base).subs(z,root)/a
            terms.append(S.simplify((S.diff(R,z,control).subs(base).subs(z,root)+S.diff(R,z,2).subs(base).subs(z,root)*dz)/a-(S.diff(p,control).subs(base).subs(z,root)+S.diff(p,z).subs(base).subs(z,root)*dz)/b))
        values.append(S.simplify(terms[0]-terms[1]))
    expected=[20*r*r/(1+r*r)**S.Rational(3,2),5*r/(1+r*r)**S.Rational(3,2),(6*r*r+1)/(1+r*r)**S.Rational(3,2)]
    assert all(S.simplify(a-b)==0 for a,b in zip(values,expected))
    return [str(q) for q in expected]

def scan(r,order):
    fractions=[.02,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,.98,.99]
    b=boundary(r,order);l=np.array([0.,20*r*r,5*r,6*r*r+1]);A=np.array([b/np.linalg.norm(b),l/np.linalg.norm(l)])
    _,sv,vh=np.linalg.svd(A,full_matrices=True);basis=vh[2:].T
    rows=[]
    for t in fractions:
        try:
            m,diag=compact(r,t,order);rows.append({'fraction':t,'M':m.tolist(),**diag})
        except Exception as e:rows.append({'fraction':t,'failure':repr(e)})
    valid=[row for row in rows if 'M' in row];M=np.array([row['M'] for row in valid]);projected=M@basis
    angles=sorted(set(float(np.arctan2(-q[0],q[1])%np.pi) for q in projected));sectors=[]
    for i,a in enumerate(angles):
        end=angles[(i+1)%len(angles)]+(np.pi if i==len(angles)-1 else 0)
        theta=(a+end)/2;direction=basis@np.array([np.cos(theta),np.sin(theta)]);direction/=max(abs(direction));values=M@direction
        # No bracket may bridge a failed energy evaluation.
        brackets=[[valid[j]['fraction'],valid[j+1]['fraction']] for j in range(len(valid)-1) if fractions.index(valid[j+1]['fraction'])==fractions.index(valid[j]['fraction'])+1 and values[j]*values[j+1]<0 and min(abs(values[j]),abs(values[j+1]))>1e-11]
        sectors.append({'direction':direction.tolist(),'brackets':brackets,'values':values.tolist(),'constraint_residuals':(A@direction).tolist()})
    return {'order':order,'boundary':b.tolist(),'rows':rows,'nullspace':basis.tolist(),'constraint_singular_values':sv.tolist(),'sectors':sectors,'max_sampled_crossings':max(len(q['brackets']) for q in sectors)}

def main():
    out={'evidence':'Exact local exponent formula; nonvalidated compact and boundary quadrature; no global root bound','exponent_derivatives_u_v_w':exact_exponent(),'records':[]}
    for r in [.125,.25,.5,1.,2.,4.,8.]:
        runs=[scan(r,order) for order in [256,512]]
        changes=[max(abs(np.array(a['M'])-b['M'])) for a,b in zip(runs[0]['rows'],runs[1]['rows']) if 'M' in a and 'M' in b]
        record={'modulus':r,'runs':runs,'maximum_absolute_quadrature_change':float(max(changes)),'boundary_change':float(max(abs(np.array(runs[0]['boundary'])-runs[1]['boundary'])))}
        out['records'].append(record);(ROOT/'q4_modulus_boundary_constraints.json').write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps({'modulus':r,'max_crossings':[q['max_sampled_crossings'] for q in runs],'quadrature_change':record['maximum_absolute_quadrature_change'],'failures':sum('failure' in row for run in runs for row in run['rows'])}),flush=True)
if __name__=='__main__':main()
