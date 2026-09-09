"""Validated displacement enclosures for the seed field's cycle brackets.

Runs the validated Taylor integrator of taylor.py through the return map of
poincare.py at the endpoints of each candidate bracket, and reports whether
the SIGN of D is rigorous there.

What a rigorous sign change does and does not give
--------------------------------------------------
Strictly signed enclosures at s0 - delta and s0 + delta establish those two
signs rigorously. Concluding that a root lies between them additionally
requires D to be defined and continuous on the whole interval, which needs
the return map validated for EVERY initial condition in [s0-delta, s0+delta],
not just the two endpoints. That step is not completed here. So this file
reports rigorous signs, not a rigorous existence proof.
"""
import json, time
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mpmath import iv, mp
from taylor import set_precision, width
import poincare as P

ROOT = Path(__file__).resolve().parent
set_precision(40)

TAU = 1e-4
A, B = mp.mpf(-7) / 4, mp.mpf(1) / 3
E0 = -4 * TAU / 11
E1 = -31379 * TAU / 25000
E2 = -7517 * TAU / 5000

# The lower bracket is NOT attempted. At |y| ~ 9355 the field has
# x_dot ~ b*y^2 ~ 3e7, and the Cartesian Taylor method drives the validated
# step below any usable floor there; a validated logarithmic reformulation is
# needed and is not implemented. Its omission is a limitation of this
# integrator, not a statement about that orbit.
ROOTS = [('upper_inner', '0.9802426972516128', 1),
         ('upper_middle', '1.3899487822846022', 1),
         ('upper_outer', '2.0452069404076005', 1)]
DELTA = mp.mpf('1e-3')


def main():
    t0 = time.perf_counter()
    out = {'evidence': 'VALIDATED interval enclosures of D at bracket endpoints; '
                       'NOT an existence proof (continuity across each bracket '
                       'is not established)',
           'precision_dps': iv.dps,
           'field': {'a': '-7/4', 'b': '1/3', 'tau': TAU, 'e0': '-4*tau/11',
                     'e1': '-31379*tau/25000', 'e2': '-7517*tau/5000'},
           'delta': str(DELTA), 'brackets': []}
    for label, s_str, side in ROOTS:
        s0 = mp.mpf(s_str)
        entry = {'label': label, 's_root_float': s_str, 'side': side}
        ends = {}
        for tag, sgn in (('minus', -1), ('plus', +1)):
            try:
                r = P.displacement(A, B, E0, E1, E2, s0 + sgn * DELTA, side)
                ends[tag] = {'D_low': iv.nstr(r['D'].a, 25),
                             'D_high': iv.nstr(r['D'].b, 25),
                             'width': float(r['width']), 'sign': r['sign'],
                             'partner_s_end':
                                 iv.nstr(r['forward']['s_end'], 20)}
            except P.Unresolved as exc:
                ends[tag] = {'unresolved': str(exc)}
        entry['endpoints'] = ends
        signs = [ends[t].get('sign') for t in ('minus', 'plus')]
        entry['rigorous_sign_change'] = (
            None not in signs and 'indeterminate' not in signs
            and signs[0] != signs[1])
        out['brackets'].append(entry)
        print(json.dumps({'label': label, 'signs': signs,
                          'rigorous_sign_change': entry['rigorous_sign_change']}),
              flush=True)
    out['rigorous_sign_changes'] = sum(1 for b in out['brackets']
                                       if b['rigorous_sign_change'])
    out['unresolved'] = [b['label'] for b in out['brackets']
                         if not b['rigorous_sign_change']]
    out['wall_seconds'] = time.perf_counter() - t0
    (ROOT / 'data' / 'certify_seed.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'rigorous_sign_changes': out['rigorous_sign_changes'],
                      'unresolved': out['unresolved']}, indent=2))


if __name__ == '__main__':
    main()
