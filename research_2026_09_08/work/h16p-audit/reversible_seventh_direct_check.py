"""Independent original Cartesian return check of two seventh-order candidates."""
from pathlib import Path
from fractions import Fraction as F
import json
import numpy as np
from scipy.integrate import solve_ivp
from reversible_finite_hopf import field
P=Path(__file__).resolve().parent
source=json.loads((P/'reversible_third_focus_probe.json').read_text());records=[]
for record in source['records']:
    expected=record['local_order7'][-1]['return_coefficients_c2_to_order'][-1]
    if abs(expected)<.001:continue
    p=record['field'];model=field(F(p['k']),F(p['m']),F(p['a']),F(p['b']));center=model[2];rhs=model[-1]
    x,y=center;a,b,e2=[float(F(p[k])) for k in ['a','b','epsilon2']];A=2*x;B=1-b+2*b*y+e2*x;omega=np.sqrt(float(F(p['determinant'])))
    rows=[]
    for radius in [.12,.08,.05,.04]:
        initial=center+np.array([radius,-A*radius/B])
        def section(t,z):return -(A*(z[0]-x)+B*(z[1]-y))/omega
        section.terminal=True;section.direction=1
        kick=solve_ivp(rhs,[0,1e-5],initial,method='DOP853',rtol=3e-14,atol=1e-15)
        sol=solve_ivp(rhs,[1e-5,3*np.pi/omega],kick.y[:,-1],method='DOP853',rtol=3e-14,atol=1e-15,events=section,max_step=.02/omega)
        assert sol.success and len(sol.t_events[0])==1
        final=sol.y_events[0][0];returned=final[0]-x;assert returned>0
        rows.append({'radius':radius,'displacement':float(returned-radius),'displacement_over_radius7':float((returned-radius)/radius**7),'time':float(sol.t_events[0][0]),'section_residual':float(section(0,final))})
    records.append({'field':p,'predicted_seventh_return_coefficient':expected,'direct_returns':rows});print(json.dumps({'a':p['a'],'predicted':expected,'rows':rows}),flush=True)
(P/'reversible_seventh_direct_check.json').write_text(json.dumps({'scope':'NUM original Cartesian Poincare returns; finite-radius ratios approach formal seventh coefficient; no interval validation','records':records},indent=2)+'\n')
