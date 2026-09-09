"""Search quadratic fields for four limit cycles in a single nest.

Why this target: Zhang Pingguang's distribution theorem forces the cycles of
a quadratic system into a (1, i) pattern -- one nest holds exactly one cycle.
Five cycles therefore require FOUR around a single focus.  Bautin caps
small-amplitude cycles at three, and Li Chengzhi (1986) forbids a cycle
surrounding an exact third-order weak focus, so the fourth must be a large
cycle appearing at finite distance from that stratum -- plausibly shed by the
outer graphic that bounds the nest.

So the objective is: three inner cycles preserved, and D(r) returning to zero
once more between the outermost cycle and the end of the return domain.

Evidence class: NUM.  A fourth sign change here would be a candidate to
certify, not a result.
"""
import json, math, os, random, sys, time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from counter import foci, scan

ROOT = Path(__file__).resolve().parent
RTOL = 1e-9
NGRID = 45
TIME_BUDGET = float(os.environ.get('SEARCH_SECONDS', '1500'))

# (lam, l, m, n, a, b, c)
SEED_CW = (-2e-5, -3.0, 0.99, 1.0, 2.0 / 9.0, -3.0, 0.0)
SEED_SHI = (-1e-14, -10.0, 4.99, 1.0, 1.0, -3113751.0 / 125000.0, 0.0)


def evaluate(coeffs):
    """Cycle counts per nest, plus the outer fold metric of the richest nest."""
    try:
        fs = foci(coeffs)
    except Exception:
        return None
    if not fs:
        return None
    nests = []
    for f in fs:
        try:
            s = scan(coeffs, f, 1e-3, 5.0, n=NGRID, rtol=RTOL)
        except Exception:
            continue
        if not s['values']:
            continue
        vals = s['values']
        outer_min = None
        if s['brackets']:
            r_out = s['brackets'][-1][1]
            tail = [d for r, d in vals if r > r_out]
            if len(tail) >= 2:
                sgn = 1.0 if tail[0] > 0 else -1.0
                outer_min = float(min(sgn * d for d in tail))
        nests.append({'centre': list(f['centre']), 'cycles': s['sign_changes'],
                      'domain_end': s['domain_end'], 'resolved_all': s['resolved_all'],
                      'outer_min': outer_min, 'brackets': s['brackets']})
    if not nests:
        return None
    nests.sort(key=lambda d: -d['cycles'])
    return {'coeffs': list(coeffs), 'nests': nests,
            'max_nest': nests[0]['cycles'],
            'total': sum(d['cycles'] for d in nests),
            'outer_min': nests[0]['outer_min']}


def sample(rng, stage):
    """Structured perturbations of the two validated three-cycle seeds."""
    base = SEED_CW if rng.random() < 0.7 else SEED_SHI
    lam, l, m, n, a, b, c = base
    if stage == 'local':
        f = rng.uniform(0.02, 0.25)
    else:
        f = rng.uniform(0.25, 0.9)
    jitter = lambda v, s: v * (1.0 + rng.uniform(-f, f)) + rng.uniform(-f, f) * s
    return (lam * math.exp(rng.uniform(-4.0, 4.0)) * rng.choice([1.0, 1.0, 0.3]),
            jitter(l, 1.0), jitter(m, 0.5), 1.0, jitter(a, 0.2), jitter(b, 1.0), 0.0)


def worker(args):
    seed, stage = args
    rng = random.Random(seed)
    return evaluate(sample(rng, stage))


def main():
    t0 = time.perf_counter()
    results, hits, best = [], [], []
    n_done = 0
    with Pool(4) as pool:
        seed = 0
        while time.perf_counter() - t0 < TIME_BUDGET:
            batch = []
            for _ in range(64):
                seed += 1
                batch.append((seed, 'local' if seed % 2 else 'wide'))
            for res in pool.imap_unordered(worker, batch):
                n_done += 1
                if res is None:
                    continue
                if res['max_nest'] >= 4 or res['total'] >= 5:
                    hits.append(res)
                    print('HIT', json.dumps(res), flush=True)
                if res['max_nest'] >= 3 and res['outer_min'] is not None:
                    best.append(res)
            if time.perf_counter() - t0 > TIME_BUDGET:
                break
    best.sort(key=lambda d: d['outer_min'])
    out = {'evidence': 'NUM; floating-point return-map sign counting',
           'target': 'four limit cycles in one nest (Zhang: distribution is (1,i))',
           'samples_evaluated': n_done, 'seconds': time.perf_counter() - t0,
           'rtol': RTOL, 'grid_points': NGRID,
           'hits_four_in_one_nest': hits,
           'three_cycle_nests_found': len(best),
           'closest_outer_folds': best[:15]}
    (ROOT / 'data' / 'search.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({k: v for k, v in out.items()
                      if k not in ('closest_outer_folds', 'hits_four_in_one_nest')},
                     indent=2), flush=True)
    print('hits:', len(hits), '| three-cycle nests:', len(best), flush=True)
    if best:
        print('best outer fold metric:', best[0]['outer_min'],
              'at', best[0]['coeffs'], flush=True)


if __name__ == '__main__':
    main()
