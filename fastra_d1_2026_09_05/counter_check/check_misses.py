"""Same-section audit of the FIVE fields with archived successful full returns.

No continuation, sweeps or counter patches. Original Fable evaluate outputs
are retained. Extended-grid samples are labeled as not visited by the counter.
"""
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
import sys,json,math,subprocess,time
from pathlib import Path
from fractions import Fraction as Q
import numpy as np
import mpmath as mp
HERE=Path(__file__).resolve().parent;D1=HERE.parent
sys.path.insert(0,str(D1))
from run_d1 import dec,st,rm,finite
mp.mp.dps=65
def save(name,row):
    with (HERE/name).open('a') as f:f.write(json.dumps(row,allow_nan=False)+'\n')
def quad(row,r,theta=0,span=0,tol=None,rounded=False):
    p=row['parameters'];m,c,b=dec(p['m']),dec(p['c']),dec(p['beta'])
    if rounded:m,c,b=[dec(str(Q.from_float(float(v)))) for v in [m,c,b]]
    tol=tol or ('2e-27' if row['label']=='positive_extension_1' else '2e-25')
    req=dict(r=st(r),c=st(c),m=st(m),beta=st(b),tol=tol,theta0=st(dec(theta)),span=st(dec(span)))
    t=time.time()
    try:
        p=subprocess.run(['/tmp/d1_full_ray'],input=' '.join(req.values())+'\n',text=True,capture_output=True,timeout=55)
        out=json.loads(p.stdout) if p.returncode==0 else dict(status='UNRESOLVED',error='exit '+str(p.returncode))
    except Exception as e:out=dict(status='UNRESOLVED',error=str(e))
    save('calls.jsonl',dict(label=row['label'],kind=row['kind'],request=req,rounded=rounded,coefficient_vector=row['coefficient_vector'],result=out,seconds=time.time()-t))
    return out
def native(row,us,theta):
    cf=np.array([float(Q(x)) for x in row['coefficient_vector']]);u=np.array([us],float)
    a,S,status=rm.returns_log(cf[None],np.zeros((1,2)),u,th0=theta,rtol=1e-12,umax=45,Smax=2000,maxsteps=300000)
    return finite([dict(u=float(x),r=float(np.exp(x)),u_return=a[0,i],D=a[0,i]-x,S=S[0,i],status=status[0,i]) for i,x in enumerate(us)])
def refine(row,lo,hi,theta,fa,fb):
    a,b=mp.log(lo),mp.log(hi);last=None
    for i in range(8):
        z=(a*fb-b*fa)/(fb-fa)
        if not a+.02*(b-a)<z<b-.02*(b-a):z=(a+b)/2
        q=quad(row,mp.exp(z),theta);last=q
        if q['status']!='NUMERICAL_ONLY':break
        f=dec(q['L'])
        if f*fa>0:a,fa=z,f
        else:b,fb=z,f
        # Positive endpoint's full map is poorly conditioned even in binary128.
        target=mp.mpf('2e-5') if row['label']=='positive_extension_1' else mp.mpf('1e-22')
        if abs(f)<target:break
    return dict(log_bracket=[st(a),st(b)],signs=[st(fa),st(fb)],approx_r=last.get('r') if last else None,approx_D=last.get('L') if last else None,return_data=last)
def check(row):
    nest=next(n for n in row['fable']['nests'] if abs(n['pt'][0])+abs(n['pt'][1])<1e-7)
    theta=nest['theta'];r0=dec(row['parameters']['fold_r']);mapped=[]
    for off in ['-.8','0','.8']:
        p=quad(row,r0*mp.exp(dec(off)),span=-theta)
        if p['status']!='NUMERICAL_ONLY':raise RuntimeError('section transfer failed')
        r=dec(p['return_coordinate']);q=quad(row,r,theta)
        if q['status']!='NUMERICAL_ONLY':raise RuntimeError('mapped full return failed')
        mapped.append(dict(horizontal_log_offset=off,transfer=p,r=st(r),full=q))
    roots=[]
    for i in range(2):
        aa,bb=mapped[i:i+2];fa,fb=dec(aa['full']['L']),dec(bb['full']['L'])
        assert fa*fb<0
        refined=refine(row,dec(aa['r']),dec(bb['r']),theta,fa,fb)
        ustar=math.log(float(refined['approx_r']));start=float(np.log(.001));j=math.floor((ustar-start)/.25)
        us=[start+j*.25,start+(j+1)*.25];grid=[]
        ns=native(row,us,theta)
        for u,nd in zip(us,ns):
            rr=mp.exp(dec(repr(u)));h=quad(row,rr,theta)
            grid.append(dict(u_float_exact_ratio=str(Q.from_float(u)),u_decimal=repr(u),native=nd,quad=h,inside_default_grid=(u<40),visited_in_original_run=any(abs(u-v)<1e-10 for call in row['fable']['calls'] if abs(call['focus'][0])+abs(call['focus'][1])<1e-7 for v in call['u'])))
        roots.append(dict(index=i+1,stability='S' if fa>0 else 'U',section_theta=repr(theta),initial_bracket=[aa['r'],bb['r']],refined=refined,grid=grid))
        print(row['label'],row['kind'],'root',i+1,'r',refined['approx_r'],'grid',[(g['u_decimal'],g['native']['status'],g['native']['D'],g['quad'].get('L')) for g in grid],flush=True)
    record=dict(label=row['label'],kind=row['kind'],parameters=row['parameters'],coefficient_vector=row['coefficient_vector'],original_fable=nest,mapped_brackets=mapped,roots=roots)
    save('missed_roots.jsonl',record)
    return record
def main():
    rows=[json.loads(s) for s in (D1/'fields.jsonl').read_text().splitlines()]
    selected=[r for r in rows if (r['label'],r['kind']) in [('center_0.0000000001','pair'),('center_0.0000000001','hopf'),('positive_extension_1','pair'),('negative_extension_1','pair'),('negative_extension_1','hopf')]]
    path=HERE/'missed_roots.jsonl';done={(r['label'],r['kind']) for r in map(json.loads,path.read_text().splitlines())} if path.exists() else set()
    for row in selected:
        if (row['label'],row['kind']) not in done:check(row)
if __name__=='__main__':main()
