"""Audit exact field identity, counting scope, and independent sign evidence."""
import json,hashlib,subprocess
from pathlib import Path
from fractions import Fraction as Q
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
def rows(name):return [json.loads(x) for x in (HERE/name).read_text().splitlines()]
fields=rows('fields.jsonl');accepted=rows('accepted.jsonl');dense=rows('dense.jsonl')
assert len({a['label'] for a in accepted})==len(accepted)==44
assert len(fields)==176 and len(dense)==88
assert len({(r['label'],r['kind']) for r in fields})==176
stages={'fold','pair','hopf_small','hopf'}
for a in accepted:
    assert {r['kind'] for r in fields if r['label']==a['label']}==stages
pair_checks=0
for row in fields:
    v=list(map(Q,row['coefficient_vector']))
    assert len(v)==12 and v[:7]==list(map(Q,[0,0,1,1,1,0,0]))
    assert v[9:11]==[Q(-10),Q(11,5)] and v[7]<0
    assert row['fable']['coefficient_vector']==row['coefficient_vector']
    if row['kind']=='fold':continue
    p=row['parameters']
    assert v[11]==Q(p['c']) and v[7]==-Q(p['m']) and v[8]==Q(p['beta'])
    expected_K=-v[7]*(11*v[11]-5)/5-42
    assert abs(expected_K-Q(p['K']))<Q(1,10**25)*max(1,abs(expected_K))
    hs=[x['result'] for x in row['half_pair_checks']]
    assert len(hs)==3 and all(x['status']=='NUMERICAL_TWO_HALF_PASSAGES' for x in hs)
    signs=[Q(x['F']) for x in hs]
    assert signs[0]*signs[1]<0 and signs[1]*signs[2]<0
    pair_checks+=1
example=json.loads((HERE/'rational_3_plus_1.json').read_text())
assert len(example['replays'])==4
for r in example['replays']:
    signs=[]
    for q in r['checks']:
        assert q['cartesian_status']==0
        assert q['log_D']*q['cartesian_D_over_r']>0
        signs.append(q['log_D'])
    assert signs[0]*signs[1]<0
    assert ('S' if signs[0]>0 else 'U')==r['stability']
for file in ['retmap.py','retmap.c','retmap_log.c','sweep_log.py']:
    path='audit/fable_engine/'+file
    original=subprocess.check_output(['git','show','fable/current:'+path],cwd=ROOT)
    assert original==(ROOT/path).read_bytes(),path
source_paths=list(HERE.glob('*.py'))+list(HERE.glob('*.cpp'))+list(HERE.glob('*.sh'))+[ROOT/'audit/fable_engine'/f for f in ['retmap.py','retmap.c','retmap_log.c','sweep_log.py']]
manifest={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(source_paths)}
(HERE/'source_manifest.json').write_text(json.dumps(manifest,indent=2))
result=dict(status='PASS',accepted=44,field_records=176,dense_profiles=88,pair_side_and_beta_fields_with_independent_two_half_sign_brackets=pair_checks,rational_example_cross_engine_brackets=4,fable_sources_unchanged=True,interval_certification=False,exhaustive_root_isolation=False)
(HERE/'audit.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result))
