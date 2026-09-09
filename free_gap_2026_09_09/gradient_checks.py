"""Check the implicit, centre-distance-normalized objective derivative."""
from pathlib import Path
import sys,json
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'moving_cycle_2026_09_09'))
import engine as e
from center_geometry import center_coefficients
src=json.loads((HERE/'strong_divided_corrected_result.json').read_text())
s=src['states'][0];amp=src['fixed_amplitude'];c=np.array(s['c']);w=np.array(s['w']);W=np.array(s['W'])

def objective(c,w):
    w=w.copy()
    for _ in range(6):
        q=np.r_[c[:2],-amp/e.TAU,w[:2]*amp/e.TAU]
        rr=[e.pair(t,q,rtol=8e-14,regularized=True) for t in w[2]+np.array([0,c[2],c[2]+c[3]])]
        assert all(r['status']=='passed' for r in rr)
        G=np.array([r['derivatives'] for r in rr]);D=np.array([r['D'] for r in rr])
        if np.max(abs(D/G[:,0]))<1e-7:break
        w+=np.linalg.solve(np.c_[G[:,4:6]*amp/e.TAU,G[:,0]],-D)
    else:raise RuntimeError('control corrector failed')
    p=e.pair(w[2]+c[2]+c[3]+c[4],q,rtol=8e-14,regularized=True)
    dist=np.linalg.norm((e.coefficients(q)[2:]-center_coefficients(c[0],c[1]))/amp)
    return p['D']/(amp*dist)

checks=[]
for j in [0,2]:
    expected=s['gradient'][j]/amp
    for h in [1e-4,3e-5]:
        d=np.zeros(5);d[j]=h
        plus=objective(c+d,w+W@d);minus=objective(c-d,w-W@d)
        fd=(plus-minus)/(2*h);error=fd-expected
        assert abs(error)<max(1e-7,abs(expected)*.003)
        checks.append({'control_index':j,'step':h,'implicit_derivative':expected,'finite_difference':fd,'error':error})
out={'evidence':'NUM derivative implementation check','checks':checks,'counts':e.COUNTS}
(HERE/'gradient_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
