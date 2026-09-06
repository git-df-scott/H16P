#!/usr/bin/env python3
"""Interrogate the one ledger showing a D_xxx sign change.
Reject it unless: the Newton residual stays small, the section stays
transversal, L stays positive (antisaddle), V1 is not collapsing to the centre,
and the crossing is a smooth transversal zero rather than a jump."""
import json
import numpy as np

for name in ('ledger_grid/cusp_c_am2p0_om0p08.jsonl', 'ledger_grid/cusp_c_am2p0_o0p08.jsonl'):
    print("="*100); print(name)
    recs = [json.loads(l) for l in open(name) if l.strip()]
    def g(r, k):
        try: return float(r[k])
        except Exception: return float('nan')
    print(" idx    x0        a11          D            Dxxx          Dxxxx        nu           res        transv    L        V1")
    for i, r in enumerate(recs):
        print(" %-4d %-9.6f %-12.6g %-12.3g %-13.6g %-12.6g %-12.6g %-10.2g %-9.5f %-8.5g %.3g"
              % (i, g(r,'x0'), g(r,'a11'), g(r,'D'), g(r,'Dxxx'), g(r,'Dxxxx'),
                 g(r,'nu') if 'nu' in r else float('nan'), g(r,'res'),
                 g(r,'transv'), g(r,'L'), g(r,'V1')))
    G = np.array([g(r,'Dxxx') for r in recs])
    s = np.sign(G)
    idx = [i for i in range(len(G)-1) if s[i]*s[i+1] < 0]
    print(" sign changes at index pairs:", idx)
    for i in idx:
        print("   between x0=%.6f (Dxxx=%.6g) and x0=%.6f (Dxxx=%.6g)"
              % (g(recs[i],'x0'), G[i], g(recs[i+1],'x0'), G[i+1]))
    print(" keys available:", sorted(recs[0].keys()))
