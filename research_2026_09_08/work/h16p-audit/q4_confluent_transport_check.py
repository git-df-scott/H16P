"""Independent area quadrature check of all four PF reconstruction columns."""
from pathlib import Path
import sys,json
import numpy as np
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'H16P'/'q4'))
from q4_integrals import alpha_beta_from_mu,basis_mp
k=2.;a=.5;d=k-1
mapping=np.column_stack([alpha_beta_from_mu(k,row)[1:] for row in np.eye(4)])
# Columns correspond to A,B,eta,c in the homogeneous four-control extension.
T=np.array([[-1,-2*k/d,0,0],[0,1/d,0,0],[0,0,-d,k],[0,0,0,1.]])
mu=np.linalg.solve(mapping,T)
source=json.loads((P/'q4_confluent_probe.json').read_text());run=next(r for r in source['records'] if r['r']==1)['runs'][-1]
rows=[]
for t in [.2,.6,.9]:
    state=next(row for row in run['rows'] if row['t']==t);X=np.array(state['states_H_Y_V_X'][3])
    pf=-a*np.pi/np.sqrt(d)*np.sqrt(1-a*t)*X/2
    area=np.array([float(x) for x in basis_mp(k,k-d*t,dps=45)])@mu
    change=area-pf
    assert max(abs(change))<2e-11
    rows.append({'t':t,'pf':pf.tolist(),'independent_area':area.tolist(),'difference':change.tolist()})
    print(json.dumps({'t':t,'maximum_absolute_difference':float(max(abs(change)))}),flush=True)
out={'scope':'Nonvalidated inherited independent area evaluator; all four universal columns checked; fixed220-step bisection and endpoint clipping remain numerical limitations','mu_matrix':mu.tolist(),'rows':rows}
(P/'q4_confluent_transport_check.json').write_text(json.dumps(out,indent=2)+'\n')
