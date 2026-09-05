"""Append-only charged supervisor. Serial use only, including failures.

The prior 206 calls remain frozen. New KKL calls have a 400-call ceiling;
all staged ledgers plus the legacy ledger count against 4096.
"""
from pathlib import Path
import hashlib,json,subprocess,sys,time
from datetime import datetime,timezone
HERE=Path(__file__).resolve().parent
LEDGER=HERE/'kkl_returns.jsonl'
def count(path):return len(path.read_text().splitlines()) if path.exists() else 0
def call(req,purpose,engine='compact_return.py'):
    used=count(LEDGER)
    total=count(HERE.parent/'kkl/data/returns.jsonl')+sum(count(p) for p in HERE.glob('*returns.jsonl'))
    if used>=400 or total>=4096:raise RuntimeError('evaluation budget exhausted')
    source=HERE/engine if engine!='legacy' else HERE.parent/'kkl/return_map.py'
    start=time.perf_counter()
    try:
        result=subprocess.run([sys.executable,str(source)],input=json.dumps(req),
                              text=True,capture_output=True,timeout=20)
        out=json.loads(result.stdout) if result.returncode==0 else dict(status='UNRESOLVED',error='process exit '+str(result.returncode),stderr=result.stderr[-1000:])
    except Exception as exc:out=dict(status='UNRESOLVED',error=str(exc))
    row=dict(stage_evaluation=used+1,prior_completed_calls=total,purpose=purpose,
             request=req,result=out,engine=str(source.relative_to(HERE.parent)),
             evaluator_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
             recorded_at_utc=datetime.now(timezone.utc).isoformat(),wall_seconds=time.perf_counter()-start)
    with LEDGER.open('a') as f:f.write(json.dumps(row,allow_nan=False)+'\n')
    print(json.dumps(dict(evaluation=used+1,purpose=purpose,status=out['status'],
                         D=out.get('log_displacement',out.get('D')),mu=out.get('multiplier',out.get('R_r')),error=out.get('error'))),flush=True)
    return out
