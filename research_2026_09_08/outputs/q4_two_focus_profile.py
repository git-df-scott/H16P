"""Finite Q4 profiles on a between-foci section, with unchanged coefficients."""
from pathlib import Path
import json
import numpy as np
from q4_finite_equilibria import equilibria
from reversible_finite_hopf import profile,count
P=Path(__file__).resolve().parent

def model(row):
    r,eps=row['r'],row['epsilon'];tau,u,w=eps*np.array(row['normalized_controls']);v=-eps
    eq=equilibria(r,eps,row['normalized_controls'])
    if eq['unhandled_x0_y1_degeneracy']:raise ArithmeticError('near degenerate elimination geometry')
    foci=[q for q in eq['equilibria'][1:] if q['determinant']>0 and q['discriminant']<0]
    if len(foci)!=1 or len(eq['equilibria'])!=2:raise ArithmeticError('not two isolated numerical foci')
    origin=np.zeros(2);remote=np.array(foci[0]['point']);distance=float(np.linalg.norm(origin-remote));e=(origin-remote)/distance;n=np.array([e[1],-e[0]])
    def rhs(t,z):
        x,y=z
        return np.array([tau*x-y-(2+r*r)*x*x+(2*r+2*u+w)*x*y+y*y,x+tau*y+(r+u)*x*x+(-1-3*r*r+v)*x*y-(r+u)*y*y])
    return row['rational_parameters'],eq['equilibria'],origin,remote,e,n,distance,rhs

def main():
    source=json.loads((P/'q4_finite_negative.json').read_text());out={'scope':'NUM same-field two-focus profiles; section joins equilibria, complete paired half-passages; no rigorous cycle count','records':[],'gates':[]}
    for row in source['records']:
        if len(row['brackets'])!=3:continue
        try:m=model(row)
        except ArithmeticError as error:
            out['gates'].append({'r':row['r'],'epsilon':row['epsilon'],'reason':str(error)});continue
        origin=profile(m,np.geomspace(.01,1e5,49),1,3e-13);remote=profile(m,np.geomspace(.01,1e7,33),-1,3e-13)
        result={'r':row['r'],'epsilon':row['epsilon'],'parameters':m[0],'equilibria':m[1],'origin':origin,'remote':remote,'origin_brackets':count(origin),'remote_brackets':count(remote),'unresolved':sum('difference' not in r for r in origin+remote)}
        out['records'].append(result);(P/'q4_two_focus_profile.json').write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps({k:result[k] for k in ['r','epsilon','origin_brackets','remote_brackets','unresolved']}),flush=True)
    (P/'q4_two_focus_profile.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
