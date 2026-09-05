import json
from scipy.optimize import brentq
from budget import call,HERE
rows=[]
def fun(K,v0=1e-4,tol=2e-12):
 b=[call(dict(K=str(K),branch=s,v0=v0,tol=tol),'infinity separatrix connection refinement',engine='graphic_shoot.py') for s in (-1,1)]
 if any(q['status']!='NUMERICAL_SEPARATRIX_PASSAGE' for q in b):raise RuntimeError('unresolved separatrix')
 d=b[0]['endpoint_x']-b[1]['endpoint_x'];rows.append(dict(K=K,v0=v0,tol=tol,splitting=d,branches=b));return d
root=brentq(fun,7,7.5,xtol=2e-11)
for v in [1e-3,1e-5,1e-6]:fun(root,v,3e-13)
(HERE/'graphic_connection.json').write_text(json.dumps(dict(root_K=root,c='8/5',rows=rows,certified=False),indent=2))
