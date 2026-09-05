"""D1 deterministic fold continuation and full compactified profiling helpers.
No sweep or descent is run. Exact rational vectors are the authoritative fields.
"""
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import sys,json,time,subprocess,hashlib
from pathlib import Path
from fractions import Fraction as Q
import numpy as np
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1]; HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'audit/fable_engine'))
# The inherited module parses CLI arguments at import; N=0, __name__ != __main__.
argv=sys.argv[:];sys.argv=['sweep_log.py','kklstar','/dev/null','0','1']
import sweep_log as sw
sys.argv=argv
rm=sw.rm;mp.mp.dps=50

def st(v): return mp.nstr(v,40)
def rat(v): return str(Q(str(v)))
def vector(c,K=None,m=None,beta='0'):
 c=Q(str(c)); M=Q(str(m)) if m is not None else 5*(Q(str(K))+42)/(11*c-5)
 return [str(v) for v in [0,0,1,1,1,0,0,-M,Q(str(beta)),-10,Q(11,5),c]]
def floats(v): return np.array([float(Q(x)) for x in v])
def append(name,obj):
 def clean(x):
  if isinstance(x,(float,np.floating)) and not np.isfinite(x):return None
  if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
  if isinstance(x,list):return [clean(v) for v in x]
  return x
 with open(HERE/name,'a') as f:f.write(json.dumps(clean(obj),allow_nan=False)+'\n');f.flush()
def half(r,c,p,chart='K',tol='2e-25'):
 src=ROOT/'fold_surface_2026_09_05'/('half_quad.cpp' if chart=='K' else 'half_m_quad.cpp')
 exe=HERE/('.'+src.stem)
 if not exe.exists():subprocess.run(['g++','-O2','-std=c++17','-fext-numeric-literals',str(src),'-o',str(exe),'-lquadmath'],check=True)
 req=dict(r=st(r),c=st(c),**{chart:st(p)},tol=tol)
 t=time.time()
 try:
  run=subprocess.run([str(exe)],input=f'{st(r)} {st(c)} {st(p)} {tol}\n',text=True,capture_output=True,timeout=12)
  res=json.loads(run.stdout) if run.returncode==0 else dict(status='UNRESOLVED',returncode=run.returncode)
 except Exception as e:res=dict(status='UNRESOLVED',error=str(e))
 append('shooting.jsonl',dict(request=req,result=res,wall_seconds=time.time()-t))
 return res

def correct(z,c,p,chart='K',fixed='p',tol='2e-25'):
 history=[]
 for _ in range(10):
  a=half(mp.exp(z),c,p,chart,tol);history.append(a)
  if a['status']!='NUMERICAL_TWO_HALF_PASSAGES':return None,history
  F=mp.matrix([a['F'],a['G']])
  keys=['z','c'] if fixed=='p' else ['c',chart]
  J=mp.matrix([[a['F_'+k] for k in keys],[a['G_'+k] for k in keys]])
  if abs(F[0])<mp.mpf('2e-21') and abs(F[1])<mp.mpf('2e-17'):
   return a,history
  try: d=mp.lu_solve(J,F)
  except Exception:return None,history
  if fixed=='p':
   if abs(d[0])>mp.mpf('.8'):d*=mp.mpf('.8')/abs(d[0])
   z-=d[0];c-=d[1]
  else:c-=d[0];p-=d[1]
 return None,history

def returns(cf,us,th=0.,rtol=1e-13,umax=60,focus=(0.,0.),maxsteps=500000):
 us=np.asarray(us,float)
 v,s,status=rm.returns_log(cf[None],np.array([focus]),us[None],th0=th,rtol=rtol,umax=umax,Smax=10000.,maxsteps=maxsteps)
 return v[0]-us,status[0],s[0]

def profile(cf,z=None,du=.1,rtol=1e-13,umin=-11.5,umax=40.,focus=(0.,0.),th=0.):
 # Every point is scanned over the whole configured domain; extra fold points only supplement it.
 us=np.arange(umin,umax+du/2,du)
 if z is not None:us=np.unique(np.r_[us,np.arange(z-.6,z+.601,.03)])
 ds=[];uu=[];edge=None
 for j in range(0,len(us),16):
  ch=us[j:j+16];d,status,S=returns(cf,ch,th,rtol,umax+12,focus)
  status=np.where((status==0)&~np.isfinite(d),5,status)
  n=len(ch) if (status==0).all() else int(np.flatnonzero(status)[0])
  uu.extend(ch[:n]);ds.extend(d[:n])
  if n<len(ch):
   edge=dict(kind='integration_failure',status=int(status[n]),u_failed=float(ch[n]),u_valid=float(uu[-1]) if uu else None)
   if uu:
    lo=uu[-1];hi=ch[n]
    for _ in range(12):
     mid=(lo+hi)/2;dd,ss,_=returns(cf,[mid],th,rtol,umax+12,focus)
     if ss[0]==0:lo=mid;de=float(dd[0])
     else:hi=mid
    if lo>uu[-1]:uu.append(lo);ds.append(de)
    edge.update(u_valid=float(lo),u_failed=float(hi))
   break
 if edge is None:edge=dict(kind='scan_cap',u_valid=float(uu[-1]))
 roots=[];uncertain=[];noise=max(5e-10,rtol*500)
 for i in range(len(ds)-1):
  if ds[i]*ds[i+1]>=0:continue
  entry=dict(u_lo=float(uu[i]),u_hi=float(uu[i+1]),D_lo=float(ds[i]),D_hi=float(ds[i+1]),stability='S' if ds[i]>0 else 'U')
  if min(abs(ds[i]),abs(ds[i+1]))<noise:uncertain.append(entry);continue
  lo,hi=uu[i:i+2];dl,dh=ds[i:i+2]
  for _ in range(22):
   mid=(lo+hi)/2;dd,ss,_=returns(cf,[mid],th,min(rtol,2e-14),umax+12,focus)
   if ss[0]:entry['refinement_failure']=int(ss[0]);break
   if dd[0]*dl>0:lo=mid;dl=dd[0]
   else:hi=mid;dh=dd[0]
  entry.update(r=float(np.exp((lo+hi)/2)),u_bracket=[float(lo),float(hi)],D_bracket=[float(dl),float(dh)])
  roots.append(entry)
 edge.update(r_valid=float(np.exp(uu[-1])) if uu else None,D=float(ds[-1]) if ds else None)
 if ds:edge['sign']=int(np.sign(ds[-1])) if abs(ds[-1])>noise else 0
 outside=None
 if roots:
  uo=roots[-1]['u_hi']+.05
  if uo<uu[-1]:
   dd,ss,_=returns(cf,[uo],th,rtol,umax+12,focus)
   outside=dict(u=float(uo),r=float(np.exp(uo)),status=int(ss[0]),D=float(dd[0]) if ss[0]==0 else None)
   if ss[0]==0:outside['sign']=int(np.sign(dd[0])) if abs(dd[0])>noise else 0
 return dict(section_angle=th,focus=list(focus),rtol=rtol,noise=noise,roots=roots,stability=''.join(r['stability'] for r in roots),uncertain_sign_changes=uncertain,edge=edge,outside_outer=outside,grid=[dict(u=float(u),D=float(d)) for u,d in zip(uu,ds)])

def field_record(label,vec,z=None,umax=40.,supplement=True):
 cf=floats(vec);z=float(z) if z is not None else None;t=time.time();baseline=sw.evaluate(cf)
 row=dict(label=label,coefficient_order='P:1,x,y,x^2,xy,y^2; Q:same',rational_vector=vec,baseline=baseline)
 if supplement:row['origin_profile']=profile(cf,z,umax=umax)
 if supplement and (z is not None and (z>=15 or -cf[7]>=10000) or len(row['origin_profile']['roots'])>=4):
  from matching import matching_profile
  print(label,'binary128 full-grid recheck',flush=True)
  row['matching_profile']=matching_profile(vec,z or 2.,umax=umax,tol='2e-29' if (z or 0)<36 else '2e-30',noise=1e-25 if (z or 0)<36 else 1e-27)
  row['origin_profile']['interpretation']='unvalidated double-precision count; use binary128 matching_profile'

 row['wall_seconds']=time.time()-t
 append('fields.jsonl',row)
 if supplement and len(row.get('matching_profile',row['origin_profile'])['roots'])>=4:
  (HERE/'FOUR_ORIGIN_TRIGGER.json').write_text(json.dumps(row,indent=2));raise RuntimeError('STOP: four origin roots require hostile reproduction')
 print(label, 'baseline',[(len(n['roots']),''.join(n['stab'])) for n in baseline], 'horizontal',row.get('matching_profile',row.get('origin_profile',{})).get('stability'), 'seconds',round(row['wall_seconds'],2),flush=True)
 return row
