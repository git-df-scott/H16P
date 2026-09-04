"""Nonrigorous 1200-bit replay of original Shi lambda=-10^-250 point."""
from pathlib import Path
import json
from gt_taylor import compute
from seed_numerics import families,run,s,x,y
rows=[]
for order in [128,144]:
 for r in ['7.0e-100','7.2e-100','2.2e-21','2.3e-21','6.5e-8','6.8e-8']:
  d=compute(r,order,1200,lambda_value='-1e-250');rows.append(d)
p=Path(__file__).resolve().parent/'data'/'shi_original_taylor.json'
p.write_text(json.dumps({'status':'NONRIGOROUS MPFR; no interval remainder','lambda':'-1e-250','samples':rows},indent=2)+'\n')
P,Q,_,_=families['gt_remote']
P=P+(s.Rational(1,10**200)-s.Rational(1,10**250))*x
families['shi_original_remote']=(P,Q,0,[(.04,.05)])
run('shi_original_remote')
