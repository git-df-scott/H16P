"""Finish requested boundary diagnostics on the SAME 24 fields; no new fields."""
import json,time
from pathlib import Path
import mpmath as mp
from full_return128 import full_returns
mp.mp.dps=60
H=Path(__file__).parent
st=lambda x:mp.nstr(x,42)
def ev(vec,u,tol='2e-24',budget=250000):
 return full_returns(vec,'0',[st(u)],y_scale='6',tolerance=tol,max_evaluations=budget)['points'][0]
rows=[json.loads(l) for l in (H/'center_sign_map.jsonl').read_text().splitlines()]
for r in rows:
 start=time.time();v=r['pair_vector'];lo=mp.mpf(40);history=[];last=r['edge']
 for u in [60,100,160,250,400,600,900,1300]:
  p=ev(v,mp.mpf(u));history.append(p)
  if p['status']!='OK_NUMERICAL':hi=mp.mpf(u);failure=p;break
  lo=mp.mpf(u);last=p
 else:hi=None
 if hi is not None:
  for _ in range(9):
   mid=(lo+hi)/2;p=ev(v,mid);history.append(p)
   if p['status']=='OK_NUMERICAL':lo=mid;last=p
   else:hi=mid;failure=p
 # Back off one log unit to avoid selecting controller failure as a boundary sign.
 sample=lo-1;p=ev(v,sample);tight=ev(v,sample,'2e-27',1000000)
 sign=int(mp.sign(mp.mpf(tight['log_displacement']))) if tight['status']=='OK_NUMERICAL' else None
 q={'index':r['index'],'K':r['K'],'vector':v,'last_success_log_radius':st(lo),'first_failure_log_radius':st(hi) if hi else None,'first_failure_status':failure['status'] if hi else None,'sample_log_radius':st(sample),'sample':p,'sample_tight':tight,'edge_sign':sign,'outside_sign':r['outside_sign'],'comparison':'agree' if sign==r['outside_sign'] else 'unresolved' if sign is None else 'DISAGREE','interpretation':'numerically resolved returning endpoint; failure beyond it is not a proof of nonreturn','history':history,'wall_seconds':time.time()-start}
 with (H/'domain_edge_map.jsonl').open('a') as f:f.write(json.dumps(q)+'\n')
 print(r['index'],r['K'],'edge bracket',st(lo),st(hi) if hi else None,q['first_failure_status'],'sign',sign,'seconds',round(time.time()-start,1),flush=True)
