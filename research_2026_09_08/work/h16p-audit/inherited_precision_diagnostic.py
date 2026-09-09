"""Reproduce inherited import-order precision loss without changing the repo."""
import hashlib
import json
import sys
from pathlib import Path
import mpmath as mp

sys.dont_write_bytecode = True
root = Path(__file__).resolve().parent
sys.path.insert(0, str(root.parent / 'H16P' / 'audit'))
mp.mp.dps = 15
import claude_green_tools as G
mp.mp.dps = 60
rows = []
for txt in ['0.2', '0.7', '0.99']:
    t = mp.mpf(txt)
    exact = mp.hyp2f1(mp.mpf(1)/6, mp.mpf(5)/6, 1, t)
    rows.append({'t': txt, 'helper_F': mp.nstr(G.Fh(t), 60),
                 'fresh_parameter_F': mp.nstr(exact, 60),
                 'difference': mp.nstr(G.Fh(t)-exact, 30)})
out = {'status': 'NUMERICAL_PARAMETER_PRECISION_DIAGNOSTIC',
       'source_sha256': hashlib.sha256(Path(G.__file__).read_bytes()).hexdigest(),
       'import_dps': 15, 'evaluation_dps': 60,
       'one6_parameter_error': mp.nstr(G.one6-mp.mpf(1)/6, 40),
       'five6_parameter_error': mp.nstr(G.five6-mp.mpf(5)/6, 40),
       'rows': rows,
       'scope': 'Import-order precision floor in helper parameters; not a disproof of the analytic theorems.'}
(root / 'inherited_precision_diagnostic.json').write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps(out, indent=2))
