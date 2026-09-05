"""Separate append-only fold-component ledger under inherited total4096."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,subprocess,sys,time
HERE=Path(__file__).resolve().parent;LEDGER=HERE/'returns.jsonl'
PRIOR=756
def used():return len(LEDGER.read_text().splitlines()) if LEDGER.exists() else 0
def call(req,purpose,engine='variational_return.py'):
    n=used()
    if PRIOR+n>=4096:raise RuntimeError('inherited 4096 total evaluation ceiling reached')
    source=HERE/engine
    if engine=='compact':
        source=HERE.parent/'staged_2026_09_05/compact_return.py'
        from fractions import Fraction as F
        c,K=F(str(req['c'])),F(str(req['K']))
        req=req|{'alpha':str(-5*(K+42)/(11*c-5))}
    start=time.perf_counter()
    try:
        p=subprocess.run([sys.executable,str(source)],input=json.dumps(req),text=True,capture_output=True,timeout=20)
        result=json.loads(p.stdout) if p.returncode==0 else dict(status='UNRESOLVED',error='process exit '+str(p.returncode),stderr=p.stderr[-500:])
    except Exception as e:result=dict(status='UNRESOLVED',error=str(e))
    row=dict(evaluation=n+1,campaign_evaluation=PRIOR+n+1,purpose=purpose,request=req,result=result,
             source=str(source.relative_to(HERE.parent)),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
             recorded_at=datetime.now(timezone.utc).isoformat(),wall_seconds=time.perf_counter()-start)
    if engine in ('angular_ld.py','half_ld.py','half_quad.py','angular_quad.py','half_m.py','angular_m_quad.py','half_m_quad.py'):
        dependency=source.with_suffix('.cpp')
        row['dependency_sha256']={dependency.name:hashlib.sha256(dependency.read_bytes()).hexdigest()}
    with LEDGER.open('a') as f:f.write(json.dumps(row,allow_nan=False)+'\n')
    result=result|{'evaluation':n+1}
    print(json.dumps(dict(eval=n+1,purpose=purpose,status=result['status'],r=result.get('r'),c=result.get('c'),K=result.get('K'),L=result.get('L',result.get('log_displacement')))),flush=True)
    return result
