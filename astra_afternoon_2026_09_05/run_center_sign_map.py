"""Deterministic centerward sheet only. 24 K values, exact rational unfoldings."""
from pathlib import Path
import json,subprocess,time,csv
from fractions import Fraction as Q
import mpmath as mp
from full_return128 import full_returns
mp.mp.dps=60
HERE=Path(__file__).resolve().parent
st=lambda x:mp.nstr(x,43)
def half(r,c,k,tol='2e-27'):
 exe=HERE/'.inherited_half_quad'
 if not exe.exists():subprocess.run(['g++','-O2','-std=c++17','-fext-numeric-literals',str(HERE/'inherited_half_quad.cpp'),'-o',str(exe),'-lquadmath'],check=True)
 p=subprocess.run([str(exe)],input=f'{st(r)} {st(c)} {st(k)} {tol}\n',text=True,capture_output=True,check=True,timeout=30)
 a=json.loads(p.stdout)
 if a['status']!='NUMERICAL_TWO_HALF_PASSAGES':raise RuntimeError(a)
 return a

def fold(k,seed):
 z=mp.log(mp.mpf(seed['r']));c=mp.mpf(seed['c']);old=mp.mpf(seed['K'])
 c+=mp.mpf(seed['F_K'])/-mp.mpf(seed['F_c'])*(k-old)
 history=[]
 for i in range(8):
  a=half(mp.exp(z),c,k);history.append(a)
  J=mp.matrix([[a['F_z'],a['F_c']],[a['G_z'],a['G_c']]])
  d=mp.lu_solve(J,mp.matrix([a['F'],a['G']]))
  if abs(d[0])<mp.mpf('2e-12') and abs(d[1])<mp.mpf('2e-23'):return a,history
  z-=d[0];c-=d[1]
 raise RuntimeError('fold correction failed')

def vector(c,k):
 c=Q(st(c));k=Q(st(k));M=5*(k+42)/(11*c-5)
 return list(map(str,[0,0,1,1,1,0,0,-M,0,-10,Q(11,5),c]))
def ev(v,us,tol='2e-25'):
 return full_returns(v,'0',[st(x) for x in us],tolerance=tol,y_scale='6',max_evaluations=250000)['points']
def value(p):return mp.mpf(p['log_displacement']) if p['status']=='OK_NUMERICAL' else None
def sign(p):
 v=value(p);return int(mp.sign(v)) if v is not None and abs(v)>mp.mpf('1e-22') else None

def main():
 events=[json.loads(l) for l in (HERE/'inherited_center_events.jsonl').read_text().splitlines()]
 ks=[mp.mpf(x) for x in ['1e-10','2e-10','5e-10','1e-9','2e-9','5e-9','1e-8','2e-8','5e-8','1e-7','2e-7','5e-7','1e-6','2e-6','5e-6','1e-5','2e-5','5e-5','0.0001','0.0002','0.0004','0.0008','0.0012','0.001953125']]
 records=[]
 for idx,k in enumerate(ks):
  start=time.time();seed=min((e['fold'] for e in events),key=lambda a:abs(mp.log(k/mp.mpf(a['K']))))
  a,h=fold(k,seed);z=mp.log(mp.mpf(a['r']));c=mp.mpf(a['c'])
  # Same inherited local fold unfolding, curvature-adapted; keeps K exact.
  cp=c-mp.mpf('.02')*mp.exp(mp.mpf(a['backward_log_sensitivity']))*mp.mpf(a['G_z'])/mp.mpf(a['F_c'])
  v=vector(cp,k)
  local_us=[z+mp.mpf(j)/20 for j in range(-10,11)]
  local=ev(v,local_us)
  brackets=[(local_us[j],local_us[j+1]) for j in range(len(local)-1) if value(local[j]) is not None and value(local[j+1]) is not None and value(local[j])*value(local[j+1])<0]
  roots=[]
  for lo,hi in brackets:
   pl=ev(v,[lo])[0]
   for j in range(27):
    mid=(lo+hi)/2;pm=ev(v,[mid])[0]
    if value(pm) is None:raise RuntimeError('root bracket lost chart')
    if value(pl)*value(pm)<=0:hi=mid
    else:lo=mid;pl=pm
   roots.append({'log_bracket':[st(lo),st(hi)],'radius':st(mp.exp((lo+hi)/2)),'endpoint_returns':ev(v,[lo,hi],'2e-28')})
  if len(roots)!=2:raise RuntimeError(('pair not recovered',st(k),len(roots)))
  outer_hi=mp.mpf(roots[-1]['log_bracket'][1]);outside_u=outer_hi+mp.mpf('.02')
  outside=ev(v,[outside_u])[0];outside_tight=ev(v,[outside_u],'2e-28')[0]
  # Full-return domain diagnostics up to e^40. Failure is NOT certified nonreturn.
  tail_us=[outside_u]+[mp.mpf(x) for x in [3,4,6,8,12,16,20,24,28,32,36,40] if mp.mpf(x)>outside_u]
  tail=ev(v,tail_us);first_fail=next((j for j,p in enumerate(tail) if value(p) is None),None)
  boundary=[]
  if first_fail is None:
   edge_u=tail_us[-1];edge=tail[-1];edge_kind='configured_grid_end_not_mathematical_boundary'
  else:
   lo=tail_us[first_fail-1];hi=tail_us[first_fail];edge=tail[first_fail-1]
   for j in range(16):
    mid=(lo+hi)/2;p=ev(v,[mid])[0];boundary.append(p)
    if value(p) is None:hi=mid
    else:lo=mid;edge=p
   edge_u=lo;edge_kind='last_success_before_unresolved_angular_or_evaluation_failure'
  edge_tight=ev(v,[edge_u],'2e-28')[0]
  row={'index':idx,'K':st(k),'fold':a,'fold_correction_history':h,'pair_vector':v,'coefficient_order':'P:1,x,y,x^2,xy,y^2; Q:same','local_profile':local,'roots':roots,'outside':outside,'outside_tight':outside_tight,'edge':edge,'edge_tight':edge_tight,'edge_kind':edge_kind,'tail':tail,'boundary_refinement':boundary,'outside_sign':sign(outside_tight),'edge_sign':sign(edge_tight),'wall_seconds':time.time()-start}
  row['comparison']='agree' if row['outside_sign'] is not None and row['edge_sign']==row['outside_sign'] else 'unresolved' if row['outside_sign'] is None or row['edge_sign'] is None else 'DISAGREE'
  records.append(row)
  with (HERE/'center_sign_map.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
  print(idx,st(k),[r['radius'][:12] for r in roots],row['outside_sign'],row['edge_sign'],edge_kind,'seconds',round(time.time()-start,1),flush=True)
 cols=['index','K','c_pair','inner_radius','outer_radius','outside_D','edge_log_radius','edge_D','outside_sign','edge_sign','comparison','edge_kind','rational_vector']
 with (HERE/'center_sign_map.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=cols);w.writeheader()
  for r in records:w.writerow(dict(index=r['index'],K=r['K'],c_pair=r['pair_vector'][-1],inner_radius=r['roots'][0]['radius'],outer_radius=r['roots'][1]['radius'],outside_D=r['outside_tight']['log_displacement'],edge_log_radius=r['edge']['log_radius'],edge_D=r['edge_tight']['log_displacement'],outside_sign=r['outside_sign'],edge_sign=r['edge_sign'],comparison=r['comparison'],edge_kind=r['edge_kind'],rational_vector=json.dumps(r['pair_vector'])))
if __name__=='__main__':main()
