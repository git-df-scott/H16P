#!/usr/bin/env python3
"""Fate of the transverse separatrices of the antipodal infinite saddles.
Terminates when the orbit gets close to the antipodal saddle (heteroclinic
connection candidate) or when it settles near a finite antisaddle."""
import numpy as np
from scipy.integrate import solve_ivp
from sphere2 import make, inf_points, tangent_eigs

def fate(l, a, T=600.0):
    F = make(l, a); pts = inf_points(l, a)
    A = [p for u, p in pts if p[0] > 0][0]
    B = [p for u, p in pts if p[0] < 0][0]
    res = {}
    for name, s0, tgt, sgn in (("A_out", A, B, +1), ("B_in", B, A, -1)):
        w, v = tangent_eigs(l, a, s0)
        k = int(np.argmax(np.abs(v[2, :]))); d = v[:, k]/np.linalg.norm(v[:, k])
        if d[2] < 0: d = -d
        s = s0 + 1e-8*d; s /= np.linalg.norm(s)
        def hit(t, y, tgt=tgt): return np.linalg.norm(y - tgt) - 1e-4
        hit.terminal = True; hit.direction = -1
        sol = solve_ivp(F, (0, sgn*T), s, rtol=1e-11, atol=1e-13, events=hit)
        d_all = np.linalg.norm(sol.y - tgt[:, None], axis=0)
        e = sol.y[:, -1]
        res[name] = dict(tmin=sol.t[-1], s3=e[2], dmin=float(d_all.min()),
                         hit=len(sol.t_events[0]) > 0,
                         xy=(e[:2]/e[2]) if abs(e[2]) > 1e-9 else None)
    return res

for (l, a) in [(-10.0, 1.0), (-12.0, 0.8), (-8.0, 1.2), (-8.0, 0.8), (-12.0, 1.2), (-9.0, 1.0)]:
    r = fate(l, a)
    print(f"l={l:6} a={a:4} | " + " | ".join(
        f"{k}: t={v['tmin']:+8.2f} s3={v['s3']:+.5f} dmin={v['dmin']:.2e} hit={v['hit']} "
        f"xy={None if v['xy'] is None else np.round(v['xy'],4)}" for k, v in r.items()))
