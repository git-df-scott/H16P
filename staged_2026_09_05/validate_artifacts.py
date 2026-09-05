"""Static evidence consistency gates; runs zero ODE evaluations."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math
HERE=Path(__file__).resolve().parent
def read(name):return json.loads((HERE/name).read_text())
def lines(name):return [json.loads(s) for s in (HERE/name).read_text().splitlines()]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

rows=lines('kkl_returns.jsonl');shi=lines('shi_returns.jsonl')
assert len(rows)==400 and len(shi)==150
assert [x['stage_evaluation'] for x in rows]==list(range(1,401))
assert sum(x['result']['status']=='UNRESOLVED' for x in rows)==6
assert sum(x['status']=='failed' for x in shi)==11
sources={sha(HERE/p) for p in ('compact_return.py','compact_return_v1.py','cartesian_check.py')}
assert all(x['evaluator_sha256'] in sources for x in rows)
for row in rows:
    req=row['request'];c,a,b=(F(str(req.get(k,0))) for k in ('c','alpha','beta'))
    assert b==0
    K=-a*(F(11,5)*c-1)-42
    assert F(1,512)-F(1,10**10)<=K<=F(6,5)+F(1,10**10)

controls=read('controls.json')
assert abs(controls[0]['old']['R_r']-.809691136)<2e-7
assert abs(controls[3]['old']['R_r']-1.902808492)<2e-7
assert max(abs(x['error']) for x in read('derivative_validation.json')[:2])<2e-5
final=read('final_verification.json');pair=read('fold_pair_coefficients.json')
assert [1 if final[i]['result']['log_displacement']>0 else -1 for i in (1,2,3)]==[1,-1,1]
for i in (1,2,3):
    a=final[i]['result'];assert abs(a['log_displacement'])>1e-7
    assert abs(a['c']-float(F(pair['field']['c'])))<1e-15
assert abs(final[0]['result']['log_displacement'])<1e-9
assert abs(final[0]['result']['multiplier']-1)<1e-8
for i in (6,7):assert abs(final[i]['result']['log_displacement'])<1e-8
assert len(read('profiles.json'))==162
# Exact remote-focus rejection at the rational pair field.
c=F(pair['field']['c']);m=-F(pair['field']['alpha']);K=F(1,512)
assert m*(F(11,5)*c-1)-42==K
J=305+634*c-11*c*c-1000*c**3
KH=-441*J/(125*(16-10*c)*(1+2*c)**2)
assert J<0 and K<KH
lo,hi=F(3324674,10**6),F(3324675,10**6)
def restoring(s):return 10*s+F(11,5)*s*s/(s-1)-c*s**3/(s-1)**2
def trace(s):return s*((1+2*c)*s/(s-1)-F(21,5))
assert restoring(lo)<m<restoring(hi)
assert 0<trace(hi)<trace(lo)<F(6,1000)
Fp_lower=10+F(11,5)*lo*(lo-2)/(hi-1)**2-c*hi*hi*(hi-3)/(lo-1)**3
assert lo*(lo-1)*Fp_lower>88
prior=len((HERE.parent/'kkl/data/returns.jsonl').read_text().splitlines())
assert prior==206
out={'status':'STATIC_NUMERICAL_EVIDENCE_CHECKS_PASS_NOT_INTERVAL_CERTIFICATE',
     'prior_evaluations':prior,'new_kkl_evaluations':len(rows),'new_shi_evaluations':len(shi),
     'total_used':prior+len(rows)+len(shi),'remaining_of_4096':4096-prior-len(rows)-len(shi),
     'kkl_historical_source_hashes_verified':True,'shi_historical_source_hashes_available':False,
     'exact_rational_pair_remote_focus':'UNSTABLE: 0<trace<0.006, determinant>88',
     'new_kkl_completed':394,'new_kkl_unresolved':6,'new_shi_completed':139,'new_shi_failed':11,
     'kkl_recorded_success_cpu_seconds':sum(x['result'].get('cpu_seconds',0) for x in rows),
     'new_ode_evaluations_from_this_check':0}
(HERE/'validation.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
