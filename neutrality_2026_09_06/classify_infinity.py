#!/usr/bin/env python3
"""Step 2. For each (l,a) on the stratum, list the real infinite singular
directions and classify them.  lam1 = -3u^2-10au+2l+5, lam2 = -(u^2+5au+l).
A hyperbolic infinite singularity is a saddle iff lam1*lam2<0."""
import numpy as np

def infsing(l, a):
    # -u^3 -5a u^2 + 2l u + 5u + a = 0  ->  u^3 +5a u^2 -(2l+5)u - a = 0
    r = np.roots([1.0, 5*a, -(2*l+5), -a])
    out = []
    for v in r:
        if abs(v.imag) < 1e-9:
            uu = v.real
            L1 = -3*uu**2 - 10*a*uu + 2*l + 5
            L2 = -(uu**2 + 5*a*uu + l)
            out.append((uu, L1, L2))
    return sorted(out)

def kind(L1, L2):
    if L1*L2 < 0: return "saddle"
    if L1*L2 > 0: return "node"
    return "degenerate"

print("stratum box of Attack 2 and a wider sweep")
for l in (-12, -10, -8, -4, -2, -1.19, -1.0, 0.0, 1.0, 4.0):
    for a in (0.8, 1.0, 1.2, 2.0):
        s = infsing(l, a)
        desc = "  ".join(f"u={uu:+.4f} {kind(L1,L2)[0]} r={abs(min(L1,L2)/max(L1,L2)):.4f}"
                         if L1*L2 < 0 else f"u={uu:+.4f} {kind(L1,L2)[0]}"
                         for uu, L1, L2 in s)
        print(f"l={l:6} a={a:4}: {len(s)} real dir | {desc}")
