"""Nonrigorous sampled continuation of the exact GT rational point."""
import json, math
from pathlib import Path
from gt_taylor import compute
from seed_numerics import families,run,s,x,y
rows=[]
for v in ['1e-8','1e-4','1e-2']:
 f=float(v);v2=abs(5*(-f)*(5-f)*(95-36*f)/12);v3=35625/8
 pairs=[(.6*1e4*f**6,.8*1e4*f**6),(.7*math.sqrt(f**4/v2),1.3*math.sqrt(f**4/v2)),(.5*math.sqrt(v2/v3),1.05*math.sqrt(v2/v3))]
 for pair in pairs:
  for r in pair:
   d=compute(format(r,'.17g'),128,900,v);d['s']=v;rows.append(d)
p=Path(__file__).resolve().parent/'data'/'gt_continuation.json'
p.write_text(json.dumps({'status':'NONRIGOROUS sampled path; no interval certificate','family':'delta=-s, epsilon=-s^4, lambda=-10^8*s^16','samples':rows},indent=2)+'\n')
u=s.Rational(1,100)
families['gt_relaxed']=(-10**8*u**16*x-y-10*x*x+(5-u)*x*y+y*y,x+x*x+(-25-8*u**4+9*u)*x*y,0,[(.04,.06)])
run('gt_relaxed')
