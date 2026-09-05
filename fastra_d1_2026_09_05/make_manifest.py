from pathlib import Path
import json,hashlib,subprocess
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
files={}
for p in sorted(HERE.rglob('*')):
 if not p.is_file() or p.name.startswith('.') or '__pycache__' in p.parts or p.name=='MANIFEST.json':continue
 files[str(p.relative_to(ROOT))]=dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
for name in ['FASTRA_D1_REPORT_2026_09_05.md','audit/fable_engine/sweep_log.py','audit/fable_engine/retmap.py','audit/fable_engine/retmap.c','audit/fable_engine/retmap_log.c','fold_surface_2026_09_05/half_quad.cpp','fold_surface_2026_09_05/half_m_quad.cpp']:
 p=ROOT/name;files[name]=dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
manifest=dict(task='FASTRA D1, 2026-09-05',fable_input_commit='4cece20',astra_fold_input_commit='7db8597cb7d9bb34e119e85bec3f229270eaf1aa',integration_parent=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),numerical_only=True,exact_rational_vectors=True,files=files)
(HERE/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Manifest:',len(files),'files')
