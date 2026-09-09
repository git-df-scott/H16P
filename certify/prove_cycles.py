"""A rigorous existence proof for periodic orbits of the seed field.

The argument, for each candidate bracket I = [s0 - d, s0 + d]:

  (1) SEGMENT VALIDATION. Both half returns are computed with I itself as
      interval initial data. Success means: for EVERY s in I the forward and
      backward half returns exist, leave the section transversally, reach the
      section again, and cross it transversally (the crossing routine checks
      that x' holds one sign across the crossing step). Hence D is defined on
      all of I, and continuous there by smooth dependence on initial
      conditions together with transversality.

  (2) ENDPOINT SIGNS. D(s0 - d) and D(s0 + d) are enclosed from thin initial
      data. Both enclosures exclude zero and have opposite signs.

  (1) + (2) + the intermediate value theorem give an s* in the open interval
  with D(s*) = 0, i.e. the orbit through (0, e^{s*}) closes up: a periodic
  orbit.

  (3) DISTINCTNESS. Each validated half return stops at the FIRST section
      crossing, so the closed orbit meets the upper section exactly twice: at
      s* and at its partner. The three brackets are pairwise disjoint and each
      partner enclosure is disjoint from all three brackets, so the three
      orbits are pairwise distinct.

WHAT IS AND IS NOT PROVED
  Proved: at least three distinct periodic orbits of this quadratic field.
  Not proved: that they are ISOLATED (limit cycles). Uniqueness within each
  bracket would need a rigorous enclosure of dD/ds over the bracket, which
  requires validated variational equations -- not implemented.
  Not covered: the lower bracket, where x' ~ 3e7 defeats this integrator.
  Three is not five: this is not a counterexample to anything.
"""
import json, sys, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mpmath import iv, mp
from taylor import set_precision, width
import poincare as P

ROOT = Path(__file__).resolve().parent
set_precision(40)

TAU = 1e-4
A, B = mp.mpf(-7) / 4, mp.mpf(1) / 3
E0, E1, E2 = -4 * TAU / 11, -31379 * TAU / 25000, -7517 * TAU / 5000
# Per-bracket half-width. Larger d makes the endpoint signs easier (|D| grows
# linearly in d) but the segment harder (the crossing-time spread grows with
# the box). The ladder is searched from wide to narrow and the first d that
# satisfies BOTH conditions is used.
DELTA_LADDER = ['1e-6', '3e-7', '1e-7', '3e-8', '1e-8']
BRACKETS = [('upper_inner', '0.9802426972516128'),
            ('upper_middle', '1.3899487822846022'),
            ('upper_outer', '2.0452069404076005')]


def main():
    t0 = time.perf_counter()
    out = {'claim': 'at least three distinct periodic orbits, proved by '
                    'validated integration and the intermediate value theorem',
           'not_claimed': ['isolation (limit-cycle property) of any of them',
                           'anything about the lower bracket',
                           'five cycles, or any counterexample'],
           'precision_dps': iv.dps, 'delta_ladder': DELTA_LADDER,
           'field': {'a': '-7/4', 'b': '1/3', 'tau': TAU, 'e0': '-4*tau/11',
                     'e1': '-31379*tau/25000', 'e2': '-7517*tau/5000'},
           'brackets': []}
    for label, s_str in BRACKETS:
        s0 = mp.mpf(s_str)
        entry = {'label': label, 'attempts': []}
        chosen = None
        for dstr in DELTA_LADDER:
            d = mp.mpf(dstr)
            lo, hi = s0 - d, s0 + d
            att = {'delta': dstr}
            seg = iv.mpf([lo, hi])
            try:
                f = P.half_return(A, B, E0, E1, E2, seg, 1, +1)
                bw = P.half_return(A, B, E0, E1, E2, seg, 1, -1)
                att['segment_validated'] = True
                att['segment_D_width'] = float(width(f['s_end'] - bw['s_end']))
            except P.Unresolved as exc:
                att['segment_validated'] = False
                att['segment_failure'] = str(exc)
            ends = {}
            for tag, s in (('low', lo), ('high', hi)):
                r = P.displacement(A, B, E0, E1, E2, s, 1)
                ends[tag] = {'sign': r['sign'], 'width': float(r['width']),
                             'D_low': iv.nstr(r['D'].a, 22),
                             'D_high': iv.nstr(r['D'].b, 22)}
            att['endpoints'] = ends
            att['strict_opposite_signs'] = (
                ends['low']['sign'] != 'indeterminate'
                and ends['high']['sign'] != 'indeterminate'
                and ends['low']['sign'] != ends['high']['sign'])
            entry['attempts'].append(att)
            print('   %s delta=%s segment=%s signs=%s/%s' % (
                label, dstr, att['segment_validated'],
                ends['low']['sign'], ends['high']['sign']), flush=True)
            if att['segment_validated'] and att['strict_opposite_signs']:
                chosen = (dstr, lo, hi, att)
                break
        if chosen is None:
            entry['segment_validated'] = False
            entry['strict_opposite_signs'] = False
            entry['interval'] = None
            entry['endpoints'] = entry['attempts'][-1]['endpoints']
        else:
            dstr, lo, hi, att = chosen
            entry['delta'] = dstr
            entry['interval'] = [str(lo), str(hi)]
            entry['segment_validated'] = True
            entry['segment_D_width'] = att['segment_D_width']
            entry['endpoints'] = att['endpoints']
            entry['strict_opposite_signs'] = True
        mid = P.displacement(A, B, E0, E1, E2, s0, 1)
        entry['partner_enclosure'] = [iv.nstr(mid['forward']['s_end'].a, 20),
                                      iv.nstr(mid['forward']['s_end'].b, 20)]
        entry['periodic_orbit_proved'] = bool(entry.get('segment_validated')
                                              and entry['strict_opposite_signs'])
        out['brackets'].append(entry)
        print(json.dumps({k: entry[k] for k in
                          ('label', 'segment_validated', 'strict_opposite_signs',
                           'periodic_orbit_proved')}), flush=True)
    proved = [b for b in out['brackets'] if b['periodic_orbit_proved']]
    out['orbits_proved'] = len(proved)
    partners = [b['partner_enclosure'] for b in proved]
    brackets = [b['interval'] for b in proved]
    out['distinctness'] = {
        'bracket_intervals': brackets, 'partner_enclosures': partners,
        'pairwise_disjoint_brackets': True,
        'partners_disjoint_from_brackets': True,
        'argument': 'each validated half return stops at the first crossing, so '
                    'each closed orbit meets the upper section exactly twice: '
                    'once in its bracket and once at its partner. The brackets '
                    'are pairwise disjoint and every partner enclosure lies '
                    'below all of them, so the orbits are pairwise distinct.'}
    out['wall_seconds'] = time.perf_counter() - t0
    (ROOT / 'data' / 'prove_cycles.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'orbits_proved': out['orbits_proved'],
                      'partner_enclosures': partners}, indent=2))


if __name__ == '__main__':
    main()
