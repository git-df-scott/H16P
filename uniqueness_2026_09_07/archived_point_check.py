"""Falsification check: does any archived (c,K) record lie in the certified region?

Reads only saved JSON from earlier strikes.  Zero ODE calls.
A single archived cycle inside the certified region would refute the theorem.
"""
import json, glob
from pathlib import Path
from fractions import Fraction
from decimal import Decimal, getcontext
getcontext().prec = 60

def min_theta_over_positive_x(c, K):
    """Coarse logarithmic scan of Theta/(21W); NUM-class screening only."""
    from decimal import Decimal as D
    c = D(str(c)); K = D(str(K))
    e = 11*c - 5; d = 16 - 10*c; m = 5*(K + 42)/e
    best = None
    for k in range(-600, 8001):
        x = D(10) ** (D(k) / D(20))
        u = 1 + x
        W = m + (2*m + 10)*x + (m + D('22.2'))*x*x + (D('12.2') - c)*x**3
        Th = m*(21 + d*x)*(u ** (c + 1)) - 21*W
        v = Th / (21*W)
        if best is None or v < best:
            best = v
    return best

def walk(o, acc):
    if isinstance(o, dict):
        if 'c' in o and 'K' in o:
            acc.append(o)
        for v in o.values():
            walk(v, acc)
    elif isinstance(o, list):
        for v in o:
            walk(v, acc)

files = sorted(glob.glob('fold_surface_2026_09_05/*.json') +
               glob.glob('fold_closure_2026_09_05/*.json') +
               glob.glob('staged_2026_09_05/*.json') +
               glob.glob('reversible_reseed/*.json'))
points = set()
for f in files:
    try:
        data = json.load(open(f))
    except Exception:
        continue
    acc = []
    walk(data, acc)
    for o in acc:
        try:
            c = float(o['c']); K = float(o['K'])
        except Exception:
            continue
        if 1.0 <= c < 1.6 and K > 0:
            points.add((round(c, 12), round(K, 12)))

inside = [(c, K) for c, K in sorted(points) if min_theta_over_positive_x(c, K) > 0]
out = {'files_scanned': len(files),
       'distinct_archived_points_in_strip': len(points),
       'max_archived_K_in_strip': max((K for _, K in points), default=None),
       'points_inside_certified_region': inside,
       'verdict': 'consistent' if not inside else 'CONTRADICTION'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
