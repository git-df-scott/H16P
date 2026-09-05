#!/usr/bin/env python3
"""Exact algebra and saved negative-sheet checks; no ODE calls."""
import json,hashlib
from pathlib import Path
from fractions import Fraction as Q
import sympy as s
HERE=Path(__file__).resolve().parent
c,u,Y=s.symbols('c u Y');x=u-1
a=2/(1-2*c);b=1/c;d0=16-10*c
assert s.factor((-2*c*(Y*Y+a*u+b)/u+a)*u*Y+2*Y*(1-u+c*Y*Y))==0
assert s.factor(s.Rational(21,5)*x-(2*c+1)*x*x/u-x*(d0*x+21)/(5*u))==0
pR=(3-4*c)/(2-4*c);pL=(2*c+1)/(2*c)
assert s.factor(pR-pL-(3*c-1)/(2*c*(1-2*c)))==0
assert s.factor((pR-pL)*(1-2*c)-(3*c-1)/(2*c))==0
assert s.diff((3*c-1)/(2*c),c).subs(c,s.Rational(1,3))==s.Rational(9,2)
CR=s.Rational(4,5)*s.Rational(38,3)*s.Rational(1,36)*s.Rational(4,5)
CL=4*s.Rational(5,3)*s.Rational(1,9)*s.Rational(2,5)
assert CR==s.Rational(152,675) and CL==s.Rational(8,27) and CL/CR==s.Rational(25,19)
mH=21*(1000*c*c+1021*c+481)/(50*(2*c+1)**2*(8-5*c))
mh,cv=s.symbols('m cv')
W=(s.Rational(61,5)-c)*u**3+(mh-s.Rational(72,5)+3*c)*u*u+(s.Rational(11,5)-3*c)*u+c
vH=5*(1+2*c)/(16-10*c)
assert s.factor(W.subs({u:-vH,mh:mH}))==0
assert mH.subs(c,s.Rational(5,11))==s.Rational(1738,75)
out={'status':'EXACT_IDENTITIES_AND_SAVED_DATA_ONLY','orbit_evaluations':0,
     'right_Melnikov_constant_c_one_third':str(CR),'left_Melnikov_constant_c_one_third':str(CL),
     'conditional_delta_log_radius_limit':'(2/9)*log(25/19)',
     'conditional_delta_log_radius_limit_decimal':str(s.N(s.Rational(2,9)*s.log(s.Rational(25,19)),20)),
     'm_H':str(mH),'files':{}}
def q(v):return Q(str(v))
def hm(c):return 21*(1000*c*c+1021*c+481)/(50*(2*c+1)**2*(8-5*c))
for name in ['events_negative.json','events_m.json','events_logm.json']:
    raw=(HERE/name).read_bytes();events=json.loads(raw);rows=[];failures=[]
    for i,e in enumerate(events):
        if e.get('status')!='ACCEPTED':failures.append({'index':i,'status':e.get('status')});continue
        f=e['fold'];cc=q(f['c']);mm=q(f['m']) if 'm'in f else -q(f['alpha'])
        assert Q(1,3)<=cc<=1 and mm>=37
        KK=mm*(11*cc-5)/5-42
        pair=e['pair_profile'];brackets=pair['root_sign_brackets']
        assert len(brackets)==2
        labels=[]
        for bb in brackets:
            assert q(bb['left']['F'])*q(bb['right']['F'])<0
            assert q(bb['left']['r'])<q(bb['right']['r'])
            labels.append(bb.get('stability','unassigned'))
        sample=next(t['result'] for t in pair['samples'] if 'alpha' in t['result'])
        pc,pm=q(sample['c']),-q(sample['alpha'])
        assert mm>hm(cc) and pm>hm(pc)
        assert q(f['G_z'])<0 and KK<0
        rows.append({'r':f['r'],'c':f['c'],'m':str(mm),'K_from_actual_m':str(KK),
                     'remote_Hopf_gap_m_minus_mH':float(mm-hm(cc)),
                     'curvature_negative':True,'root_brackets':len(brackets),'stability_labels':labels})
    out['files'][name]={'sha256':hashlib.sha256(raw).hexdigest(),'accepted':len(rows),
                       'failures':failures,'rows':rows}
path=HERE/'large_m_full_reproduction.json'
if path.exists():
    raw=path.read_bytes(); replay=json.loads(raw); rows=[]; fields=set()
    for rr in replay['rows']:
        left,right=[z['result'] for z in rr['checks'][:2]]
        assert left['status']==right['status']=='NUMERICAL_ONLY'
        assert q(left['L'])*q(right['L'])<0 and q(left['r'])<q(right['r'])
        assert q(left['L_z'])*q(right['L_z'])>0
        for z in [left,right]:fields.add((q(z['c']),q(z['m'])))
        rows.append({'left_radius':left['r'],'right_radius':right['r'],
                     'left_L':left['L'],'right_L':right['L'],
                     'endpoint_derivative_sign':'positive' if q(left['L_z'])>0 else 'negative'})
    assert len(fields)==1
    cc,mm=next(iter(fields));assert mm>hm(cc)
    out['full_return_reproduction']={'sha256':hashlib.sha256(raw).hexdigest(),
          'field':{'c':str(cc),'m':str(mm)},'two_brackets_verified':len(rows)==2,
          'rows':rows,'interval_certified':False,'third_cycle_trigger':False}
(HERE/'theory_negative_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='files'},indent=2))
for name,z in out['files'].items():print(name,'accepted',z['accepted'],'failed',len(z['failures']))
