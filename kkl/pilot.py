"""Serial bounded numerical pilot; every evaluator request is recorded.

No interval certification or full parameter-space coverage is claimed.
"""
import argparse
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
LEDGER=HERE/'data'/'returns.jsonl'


def call(request,purpose):
    LEDGER.parent.mkdir(exist_ok=True)
    count=sum(1 for _ in LEDGER.open()) if LEDGER.exists() else 0
    if count>=4096: raise RuntimeError('full-strike evaluation budget exhausted')
    start=time.perf_counter()
    try:
        child=subprocess.run([sys.executable,str(HERE/'return_map.py')],input=json.dumps(request),
                             text=True,capture_output=True,timeout=20)
        result=json.loads(child.stdout) if child.returncode==0 else {
            'status':'UNRESOLVED','error':'evaluator exit '+str(child.returncode),'stderr':child.stderr[-1000:]}
    except Exception as exc:
        result={'status':'UNRESOLVED','error':str(exc)}
    record={'evaluation':count+1,'purpose':purpose,'request':request,'result':result,
            'total_wall_seconds':time.perf_counter()-start,
            'recorded_at_utc':datetime.now(timezone.utc).isoformat(),
            'evaluator_sha256':hashlib.sha256((HERE/'return_map.py').read_bytes()).hexdigest()}
    with LEDGER.open('a') as f: f.write(json.dumps(record,allow_nan=False)+'\n')
    print(json.dumps({key:record[key] for key in ('evaluation','purpose')} |
                     {key:result.get(key) for key in ('status','r','D','R_r','R_c','R_alpha','error')}),flush=True)
    return result


def used_steps():
    total=0
    for path in (HERE/'data').glob('*.jsonl'):
        if path.name in ('returns.jsonl','continuation_events.jsonl'): continue
        total+=sum(1 for _ in path.open())
    return total


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--request',required=True)
    parser.add_argument('--purpose',required=True)
    args=parser.parse_args()
    call(json.loads(args.request),args.purpose)
