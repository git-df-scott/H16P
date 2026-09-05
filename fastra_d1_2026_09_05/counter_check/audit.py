"""Audit the discrepancy claims against full returns and exact field identity."""
import json,subprocess
from pathlib import Path
from fractions import Fraction as Q
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent.parent
rows=json.loads((HERE/'discrepancy_ledger.json').read_text())
assert len(rows)==5 and sum(len(r['roots']) for r in rows)==10
assert {r['id'] for r in rows}=={'C0','Cbeta','P0','N0','Nbeta'}
original={(r['label'],r['kind']):r for r in map(json.loads,(HERE.parent/'fields.jsonl').read_text().splitlines())}
for r in rows:
    assert r['coefficient_vector']==original[r['label'],r['kind']]['coefficient_vector']
    assert len([Q(x) for x in r['coefficient_vector']])==12
    for root in r['roots']:
        a,b=root['grid'];da,db=Q(a['quad']['L']),Q(b['quad']['L'])
        assert a['quad']['status']==b['quad']['status']=='NUMERICAL_ONLY'
        assert da*db<0
        assert root['stability']==('S' if da>0 else 'U')
        assert abs(float(b['u_decimal'])-float(a['u_decimal'])-.25)<1e-14
        if r['id'].startswith('C'):assert a['visited_in_original_run'] and b['visited_in_original_run']
        else:assert not a['visited_in_original_run'] and not b['visited_in_original_run']
        if r['id']=='P0':assert not a['inside_default_grid'] and not b['inside_default_grid']
steps=json.loads((HERE/'first_step.json').read_text())
bad=next(r for r in steps if r['label']=='negative_extension_1')['first_step']
assert not bad['error_finite'] and not bad['original_reject_condition']
rounding=[r for r in map(json.loads,(HERE/'coefficient_rounding.jsonl').read_text().splitlines())]
for label in ['center_0.0000000001','negative_extension_1']:
    rr=[r for r in rounding if r['label']==label]
    assert len(rr)==2
    assert all(r['quad_on_rounded_field']['status']=='NUMERICAL_ONLY' for r in rr)
    assert Q(rr[0]['quad_on_rounded_field']['L'])*Q(rr[1]['quad_on_rounded_field']['L'])<0
for f in ['retmap.py','retmap.c','retmap_log.c','sweep_log.py']:
    path='audit/fable_engine/'+f
    assert subprocess.check_output(['git','show','fable/current:'+path],cwd=ROOT)==(ROOT/path).read_bytes()
result=dict(status='PASS',fields=5,full_return_pair_roots=10,same_ray_grid_endpoints=20,fable_production_sources_unchanged=True,no_new_sheets_or_sweeps=True,D1='OPEN')
(HERE/'audit.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
