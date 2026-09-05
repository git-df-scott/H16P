"""D1 deterministic fold continuation and full compactified profiles.

No random sweep or descent. Fable evaluate is called without modification.
Its sampled count is a lower bound, not an exhaustive zero isolation theorem.
All fields are recorded as exact rational vectors before float conversion.
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import sys, json, math, time, subprocess
from pathlib import Path
from fractions import Fraction as Q
import mpmath as mp
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
OLD=ROOT/'fold_surface_2026_09_05'
mp.mp.dps=65
sys.path.insert(0,str(ROOT/'audit/fable_engine'))
saved=sys.argv[:];sys.argv=['d1','kklstar','unused','0','0']
import sweep_log as sl
import retmap as rm
sys.argv=saved
def dec(v):
    q=Q(str(v));return mp.mpf(q.numerator)/q.denominator
def st(v):return mp.nstr(v,50)
def vector(c,m,b=0):return list(map(str,[0,0,1,1,1,0,0,-Q(str(m)),Q(str(b)),-10,Q(11,5),Q(str(c))]))
def append(name,row):
    with (HERE/name).open('a') as f:f.write(json.dumps(row,allow_nan=False)+'\n')
def hp(r,c,m,b=0,tol='2e-25',full=False,purpose=''):
    req=dict(r=st(dec(r)),c=st(dec(c)),m=st(dec(m)),beta=st(dec(b)),tol=tol)
    try:
        p=subprocess.run(['/tmp/d1_full' if full else '/tmp/d1_half'],input=' '.join(req.values())+'\n',text=True,capture_output=True,timeout=float(os.environ.get('D1_HP_TIMEOUT','25')))
        a=json.loads(p.stdout) if p.returncode==0 else dict(status='UNRESOLVED',error='process exit '+str(p.returncode))
    except Exception as e:a=dict(status='UNRESOLVED',error=str(e))
    append('hp_calls.jsonl',dict(purpose=purpose,request=req,coefficient_vector=vector(req['c'],req['m'],req['beta']),full=full,result=a))
    return a
def correct(r,c,t,chart='K'):
    z=mp.log(dec(r));c=dec(c);t=dec(t);hist=[]
    for it in range(12):
        m=5*(t+42)/(11*c-5) if chart=='K' else t
        a=hp(mp.exp(z),c,m,purpose='fold corrector');hist.append(a)
        if a['status']!='NUMERICAL_TWO_HALF_PASSAGES':return None,hist
        f=mp.matrix([a['F'],a['G']]);fz=dec(a['F_z']);gz=dec(a['G_z'])
        mc=-11*m/(11*c-5) if chart=='K' else 0
        fc=dec(a['F_c'])+mc*dec(a['F_m']);gc=dec(a['G_c'])+mc*dec(a['G_m'])
        if abs(f[0])<mp.mpf('2e-22') and abs(f[1])<mp.mpf('2e-19'):
            a['coefficient_vector']=vector(a['c'],a['m']);a['chart']=chart;return a,hist
        d=mp.lu_solve(mp.matrix([[fz,fc],[gz,gc]]),f)
        factor=min(mp.mpf(1),mp.mpf('.6')/max(abs(d[0]),mp.mpf('1e-80')),mp.mpf('.08')/max(abs(d[1]),mp.mpf('1e-80')))
        z-=factor*d[0];c-=factor*d[1]
    return None,hist
def correct_fixed_r(r,c,m):
    c=dec(c);m=dec(m);hist=[]
    for it in range(10):
        a=hp(r,c,m,tol='2e-29',purpose='fixed-radius fold continuation');hist.append(a)
        if a['status']!='NUMERICAL_TWO_HALF_PASSAGES':return None,hist
        f=mp.matrix([a['F'],a['G']])
        curvature=mp.exp(dec(a['backward_log_sensitivity']))*abs(dec(a['G_z']))
        if abs(f[0])<min(mp.mpf('2e-25'),curvature*mp.mpf('1e-5')) and abs(f[1])<mp.mpf('2e-19'):
            a['coefficient_vector']=vector(a['c'],a['m']);a['chart']='r';return a,hist
        d=mp.lu_solve(mp.matrix([[a['F_c'],a['F_m']],[a['G_c'],a['G_m']]]),f);c-=d[0];m-=d[1]
    return None,hist

def finite(a):
    if isinstance(a,np.ndarray):return finite(a.tolist())
    if isinstance(a,(list,tuple)):return [finite(x) for x in a]
    if isinstance(a,dict):return {k:finite(v) for k,v in a.items()}
    if isinstance(a,(float,np.floating)):return float(a) if np.isfinite(a) else None
    if isinstance(a,np.integer):return int(a)
    return a
def fable_count(v,normalize=False):
    cf=[Q(x) for x in v];scale=Q(1)
    if normalize:
        scale=Q(10)**int(math.floor(math.log10(float(-cf[7]))/2))
        # x=X, y=sY, tau=s*t; exact rational conjugacy and positive time change.
        cf=[cf[0]/scale,cf[1]/scale,cf[2],cf[3]/scale,cf[4],cf[5]*scale,
            cf[6]/scale**2,cf[7]/scale**2,cf[8]/scale,cf[9]/scale**2,cf[10]/scale,cf[11]]
    c=np.array([float(x) for x in cf]);calls=[];original=rm.returns_log
    def capture(coef,foc,u0,**kw):
        result=original(coef,foc,u0,**kw)
        calls.append(finite(dict(focus=np.asarray(foc)[0],theta=kw.get('th0',0),u=np.asarray(u0)[0],D=(result[0]-u0)[0],status=result[2][0])))
        return result
    rm.returns_log=capture
    try:result=sl.evaluate(c)
    finally:rm.returns_log=original
    nests=[]
    for n in result:
        cs=[q for q in calls if np.linalg.norm(np.array(q['focus'])-np.array(n['pt']))<1e-8*(1+np.linalg.norm(n['pt']))]
        points={u:(d,s) for q in cs for u,d,s in zip(q['u'],q['D'],q['status'])}
        good=sorted((u,d) for u,(d,s) in points.items() if s==0 and d is not None)
        fail=sorted(u for u,(d,s) in points.items() if s!=0)
        edge=math.log(n['redge']) if n['redge'] else None
        # Use exactly the edge returned by the counter; ignore valid samples beyond its first failure.
        selected=[(u,d) for u,d in good if edge is not None and u<=edge+1e-10]
        outside=None
        if n['roots']:
            outer=math.log(n['roots'][-1]);outside=next((dict(log_r=u,D=d) for u,d in selected if u>outer),None)
        n.update(edge_kind='integration_failure_bracket' if fail else 'scan_cap',edge_D=selected[-1][1] if selected else None,outside_outer=outside,theta=cs[0]['theta'] if cs else None,first_failure_log_r=fail[0] if fail else None)
        nests.append(n)
    return finite(dict(coefficient_vector=list(map(str,cf)),coordinate_scale=str(scale),nests=nests,calls=calls,exhaustive=False))

def profile(a,label):
    c=dec(a['c']);m=dec(a['m']);r=dec(a['r']);K=dec(a['K'])
    fold=dict(label=label,kind='fold',coefficient_vector=vector(st(c),st(m)),parameters={k:a[k] for k in ['c','K','m','r']})
    fold['fable']=fable_count(fold['coefficient_vector'])
    if m>1000:fold['scaled_fable']=fable_count(fold['coefficient_vector'],True)
    append('fields.jsonl',fold)
    curvature=mp.exp(dec(a['backward_log_sensitivity']))*dec(a['G_z'])
    # Fixed actual m transverse displacement. Local expected roots at log offsets +/-0.4.
    cp=c-mp.mpf('.08')*curvature/dec(a['F_c'])
    pair=hp(r,cp,m,purpose='pair side center test')
    if pair['status']!='NUMERICAL_TWO_HALF_PASSAGES':return
    btrial=-mp.sign(K)*min(mp.mpf('1e-6'),abs(K)*mp.mpf('1e-5'))
    trial=hp(r,cp,m,btrial,purpose='beta response calibration')
    if trial['status']!='NUMERICAL_TWO_HALF_PASSAGES':return
    Fb=(dec(trial['F'])-dec(pair['F']))/btrial
    bcap=-mp.sign(K)*min(abs(btrial),abs(dec(pair['F'])/Fb)/10)
    for kind,b in [('pair',mp.mpf(0)),('hopf_small',bcap/10),('hopf',bcap)]:
        row=dict(label=label,kind=kind,coefficient_vector=vector(st(cp),st(m),st(b)),parameters=dict(c=st(cp),m=st(m),K=st(m*(11*cp-5)/5-42),beta=st(b),fold_r=st(r)),origin_equilibrium_stability=('S' if b<0 else 'U') if b else ('U' if K>0 else 'S'))
        row['fable']=fable_count(row['coefficient_vector'])
        if m>1000:row['scaled_fable']=fable_count(row['coefficient_vector'],True)
        # Independent brackets at +/- .8 and at fold radius; count still comes from full scan above.
        checks=[]
        for off in [-.8,0,.8]:
            h=hp(r*mp.exp(off),cp,m,b,purpose=kind+' pair bracket')
            checks.append(dict(log_offset=off,result=h))
        row['half_pair_checks']=checks
        append('fields.jsonl',row)
        ns=row.get('scaled_fable',row['fable'])['nests'];origin=next((n for n in ns if np.linalg.norm(n['pt'])<1e-7),{})
        print(label,kind,'c',st(cp)[:18],'m',st(m)[:15],'beta',st(b)[:13],'origin',origin.get('stab'),'edge',origin.get('edge_D'),flush=True)
        if len(origin.get('roots',[]))>=4:
            append('triggers.jsonl',row);raise RuntimeError('FOUR_ORIGIN_TRIGGER_REQUIRES_REVIEW')

def seeds():
    out=[]
    # Correct all selected inherited checkpoints, and insert new targets using nearest seeds.
    inc=json.load(open(OLD/'events_increasing.json'))
    for i in [0,2,4,6,8,10,12,14,15,16,17]:
        a=inc[i]['fold'];out.append(('positive_K_'+str(a['K']),a['r']*.973,a['c'],a['K'],'K'))
    center=json.load(open(OLD/'center_binary128.json'))['rows'][-1]['result']
    for k in ['0.001953125','0.0001','0.00001','0.000001','0.00000001','0.0000000001']:
        out.append(('center_'+k,center['r'],st(dec('0.96862063355349428616412539953798547325')+dec(k)*dec('.137109611')),k,'K'))
    for file,indices in [('events_half.json',[0,3,7]),('events_quad.json',[0,3,6,9]),('events_negative.json',[0,3,6,9,12,15,17,19]),('events_m.json',[0,2,4,5]),('events_logm.json',[0,3,6,10])]:
        es=[e for e in json.load(open(OLD/file)) if e.get('status')=='ACCEPTED']
        for i in indices:
            a=es[i]['fold'];m=a.get('m') or st(5*(dec(a['K'])+42)/(11*dec(a['c'])-5))
            out.append((file[:-5]+'_'+str(i),a['r'],a['c'],m,'m'))
    return out

def main():
    done=set()
    if (HERE/'continuation.jsonl').exists():done={x['label'] for x in map(json.loads,(HERE/'continuation.jsonl').read_text().splitlines()) if x.get('completed')}
    for label,r,c,t,chart in seeds():
        if label in done:continue
        a,h=correct(r,c,t,chart)
        if a is None:append('continuation.jsonl',dict(label=label,status='UNRESOLVED_CORRECTOR',history=h));continue
        append('accepted.jsonl',dict(label=label,fold=a,history=h))
        profile(a,label)
        append('continuation.jsonl',dict(label=label,completed=True,fold=a))
    # Extensions beyond both inherited large-radius endpoints.
    for file,kind in [('events_quad.json','positive_extension'),('events_logm.json','negative_extension')]:
        a=[e['fold'] for e in json.load(open(OLD/file)) if e.get('status')=='ACCEPTED'][-1]
        for i in range(2):
            label=kind+'_'+str(i)
            m=a.get('m') or st(5*(dec(a['K'])+42)/(11*dec(a['c'])-5))
            if kind=='positive_extension':b,h=correct_fixed_r(dec(a['r'])*mp.exp(1),a['c'],m)
            else:b,h=correct(dec(a['r'])*mp.exp(1.3),a['c'],dec(m)*mp.exp(1.5),'m')
            if b is None:append('continuation.jsonl',dict(label=label,status='UNRESOLVED_CORRECTOR',history=h));break
            a=b;append('accepted.jsonl',dict(label=label,fold=a,history=h));profile(a,label);append('continuation.jsonl',dict(label=label,completed=True,fold=a))

if __name__=='__main__':main()
