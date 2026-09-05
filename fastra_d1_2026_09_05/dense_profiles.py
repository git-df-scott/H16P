"""Deterministic origin profiles supplement the unchanged Fable counter.

Two tolerance levels, positive horizontal ray, rational scaling for large m,
smaller starting radii, and a longer scan. Disagreement is marked unresolved.
These are sampled numerical profiles, not interval root counts.
"""
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
import json, math, sys
from pathlib import Path
from fractions import Fraction as Q
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from run_d1 import HERE,rm,finite

def work(row):
    v=[Q(x) for x in row['coefficient_vector']];m=-v[7]
    s=Q(10)**max(0,int(math.floor(math.log10(float(m))/2)))
    c=np.array([0,0,1,float(1/s),1,0,0,float(-m/s**2),float(v[8]/s),float(-10/s**2),float(Q(11,5)/s),float(v[11])])
    grid=np.arange(-25.,46.0001,.125);runs=[]
    for tol in [1e-12,1e-13]:
        points=[];boundary=None
        for i in range(0,len(grid),16):
            us=grid[i:i+16]
            u1,_,stat=rm.returns_log(c[None],np.zeros((1,2)),us[None],th0=0,rtol=tol,umax=60,Smax=10000,maxsteps=1000000)
            fail=np.flatnonzero(stat[0]!=0);n=int(fail[0]) if len(fail) else len(us)
            points.extend((float(u),float(d)) for u,d in zip(us[:n],(u1[0]-us)[:n]))
            if n<len(us):boundary=dict(log_r=float(us[n]),status=int(stat[0,n]));break
        runs.append(dict(rtol=tol,points=points,boundary=boundary))
    common=min(len(runs[0]['points']),len(runs[1]['points']));samples=[]
    for (u,d0),(_,d1) in zip(runs[0]['points'][:common],runs[1]['points'][:common]):
        # Empirical resolution check, explicitly not a rigorous error bound.
        e=max(5e-12,10*abs(d0-d1));sig=1 if d1>e else -1 if d1 < -e else 0
        samples.append(dict(log_r=u,D=d1,tolerance_difference=d1-d0,resolved_sign=sig))
    brackets=[]
    previous=None
    for p in samples:
        if not p['resolved_sign']:continue
        if previous and p['resolved_sign']!=previous['resolved_sign']:
            brackets.append(dict(left=previous,right=p,stability='S' if previous['resolved_sign']>0 else 'U'))
        previous=p
    outside=brackets[-1]['right'] if brackets else None
    edge=samples[-1] if samples else None
    mismatch=bool(outside and edge and outside['resolved_sign']*edge['resolved_sign']==-1)
    return finite(dict(label=row['label'],kind=row['kind'],coefficient_vector=row['coefficient_vector'],coordinate_scale=str(s),section='y=0,x>0; x unchanged by scaling',scan_log_radius=[-25,46],runs=runs,samples=samples,brackets=brackets,stability=''.join(b['stability'] for b in brackets),outside_outer=outside,edge=edge,edge_kind='integration_failure' if runs[-1]['boundary'] else 'scan_cap',sign_mismatch=mismatch,exhaustive=False))

def main():
    rows=[json.loads(s) for s in (HERE/'fields.jsonl').read_text().splitlines()]
    done=set()
    path=HERE/'dense.jsonl'
    if path.exists():done={(r['label'],r['kind']) for r in map(json.loads,path.read_text().splitlines())}
    rows=[r for r in rows if r['kind'] in ['pair','hopf'] and (r['label'],r['kind']) not in done]
    with ProcessPoolExecutor(max_workers=3) as pool:
        for result in pool.map(work,rows):
            with path.open('a') as f:f.write(json.dumps(result,allow_nan=False)+'\n')
            print(result['label'],result['kind'],result['stability'],result['edge_kind'],flush=True)

if __name__=='__main__':main()
