"""Bounded floating exploration of the open two-anchor determinant. No proof."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ[k]='1'
import sys, json, resource, time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import mpmath as mp
from q4_threshold_path import primitive_basis_closed
from q4_reconstruction import reconstruct, green_coordinates, center_data
resource.setrlimit(resource.RLIMIT_CPU,(60,60)); os.nice(10)
mp.mp.dps=55

def coefficients(r,s):
    rows=[primitive_basis_closed(t) for t in (r,s,1)]
    B=tuple(mp.lu_solve(mp.matrix([[x[0],x[1],-x[2]] for x in rows]),mp.matrix([x[0]-x[3] for x in rows])))
    bv=tuple(mp.lu_solve(mp.matrix([[x[0],x[1]] for x in rows[:2]]),mp.matrix([-x[2] for x in rows[:2]])))
    v=mp.mpf(9)/3080*(bv[0]+mp.mpf(144)/221*bv[1]+mp.mpf(11)/6)
    yb=center_data(mp.mpf('.5'),*B)[0]; lc=-yb/v
    C=(B[0]+lc*bv[0],B[1]+lc*bv[1],B[2]-lc)
    return B,C,yb,v,lc

def main():
    records=[]; start=time.process_time()
    for r in map(mp.mpf,('.25','.6','.9','.99','.9999')):
      for gap in map(mp.mpf,('.01','.5','.99')):
        s=r+(1-r)*gap
        B,C,yb,v,lc=coefficients(r,s)
        for a in (.5,.875,.95,.999):
          bsol=reconstruct(a,*map(float,B),t_end=float(r)); csol=reconstruct(a,*map(float,C),t_end=float(r))
          zb,pb=map(float,green_coordinates(a,bsol,float(r))); zc,pc=map(float,green_coordinates(a,csol,float(r)))
          det=pb*zc-zb*pc
          row=dict(r=str(r),s=str(s),a=a,P_B=pb,Z_B=zb,P_C=pc,Z_C=zc,K=det,first_gate=pb>0 and det>0)
          records.append(row)
          if row['first_gate']: print('PASS FIRST GATE',json.dumps(row),flush=True)
    record=dict(status='NUMERICAL_ONLY',cpu_seconds=time.process_time()-start,records=records)
    Path(__file__).with_name('determinant_exploration.json').write_text(json.dumps(record,indent=2)+'\n')
    print('SUMMARY',len(records),'points',sum(x['first_gate'] for x in records),'pass first gate; seconds',record['cpu_seconds'])
    print('largest K',max(records,key=lambda x:x['K']))
if __name__=='__main__':main()
