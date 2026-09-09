"""Validate the counter against the published Chen-Wang four-cycle field."""
import json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from counter import foci, scan, refine, equilibria, jacobian
import numpy as np

d1, d2 = 1e-2, 2e-5
COEFFS = (-d2, -3.0, 1.0 - d1, 1.0, 2.0 / 9.0, -3.0, 0.0)

t0 = time.perf_counter()
out = {'source': 'Chen-Wang visualization parameters (Yu-Zeng eq. 12/15)',
       'coeffs': dict(zip(('lam', 'l', 'm', 'n', 'a', 'b', 'c'), COEFFS)),
       'equilibria': [{'x': p[0], 'y': p[1],
                       'eigs': [str(v) for v in np.linalg.eigvals(jacobian(COEFFS, p))]}
                      for p in equilibria(COEFFS)],
       'nests': []}
for f in foci(COEFFS):
    s = scan(COEFFS, f, 1e-3, 5.0, n=90)
    roots = [r for r in (refine(COEFFS, f, b) for b in s['brackets']) if r is not None]
    out['nests'].append({'centre': list(f['centre']), 'sense': f['sense'],
                         'trace': f['trace'], 'sign_changes': s['sign_changes'],
                         'resolved_all': s['resolved_all'],
                         'domain_end': s['domain_end'], 'cycle_radii': roots})
    print(json.dumps(out['nests'][-1]), flush=True)
out['total_cycles'] = sum(len(n['cycle_radii']) for n in out['nests'])
out['distribution'] = sorted((len(n['cycle_radii']) for n in out['nests']), reverse=True)
out['wall_seconds'] = time.perf_counter() - t0
Path(__file__).resolve().parent.joinpath('data', 'validate.json').write_text(
    json.dumps(out, indent=2) + '\n')
print('total', out['total_cycles'], 'distribution', out['distribution'],
      '%.1fs' % out['wall_seconds'])
