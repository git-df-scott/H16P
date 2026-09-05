"""Frozen inherited accounting:756+3297+24=4077,19 calls remain."""
import json,time,hashlib,subprocess,sys
from pathlib import Path
H=Path(__file__).resolve().parent;P=H/'returns.jsonl';PRIOR=4077

def call(q,purpose,engine=None):
 n=len(P.read_text().splitlines()) if P.exists() else 0
 if n>=19:raise RuntimeError('shared4096 evaluation cap reached')
 source=H/'cusp_return.py' if engine is None else H.parent/'fold_surface_2026_09_05'/engine;t=time.perf_counter()
 try:
  p=subprocess.run([sys.executable,str(source)],input=json.dumps(q),text=True,capture_output=True,timeout=20)
  a=json.loads(p.stdout) if p.returncode==0 else dict(status='UNRESOLVED',error=p.stderr[-300:])
 except Exception as e:a=dict(status='UNRESOLVED',error=str(e))
 row=dict(evaluation=n+1,campaign_evaluation=PRIOR+n+1,purpose=purpose,request=q,result=a,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),cpp_sha256=hashlib.sha256(source.with_suffix('.cpp').read_bytes()).hexdigest(),wall_seconds=time.perf_counter()-t)
 with P.open('a') as f:f.write(json.dumps(row)+'\n')
 print(json.dumps(dict(evaluation=n+1,purpose=purpose,result=a)),flush=True)
 return a
