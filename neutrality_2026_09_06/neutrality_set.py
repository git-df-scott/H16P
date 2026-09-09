#!/usr/bin/env python3
"""Region of the order-3 stratum where infinity carries THREE real directions.
Then two of them are saddles that are NOT antipodes of one another, and a
graphic
    equator arc -> S1 -> plane orbit -> S2 -> equator arc -> ...
has first stability coefficient
    r = |lam_eq(u1)/lam_tr(u1)| * |lam_tr(u2)/lam_eq(u2)|
which is a genuine codimension-one condition.  We locate r=1 and compare its
sign requirement with eta_3 (Li-Cherkas)."""
import numpy as np

def data(l, a):
    r = np.roots([1.0, 5*a, -(2*l+5), -a])
    us = sorted(v.real for v in r if abs(v.imag) < 1e-9)
    out = []
    for u in us:
        le = -3*u*u - 10*a*u + 2*l + 5      # along the equator
        lt = -(u*u + 5*a*u + l)             # transverse, into the plane
        out.append((u, le, lt))
    return out

def eta3(l, a):
    return -25*a*(2*a**2+l+2)*(5*a**2*l+6*a**2-3*l**3-12*l**2-15*l-6)/64

def ratio_pair(l, a):
    d = data(l, a)
    if len(d) != 3: return None
    sad = [(u, le, lt) for u, le, lt in d if le*lt < 0]
    if len(sad) != 2: return None
    (u1, e1, t1), (u2, e2, t2) = sad
    r = abs(e1/t1) * abs(t2/e2)
    return r, u1, u2

print("scanning the three-direction region of the stratum for r=1")
found = []
LS = np.arange(-40.0, 20.01, 0.5)
AS = np.arange(0.1, 6.01, 0.05)
grid = np.full((len(LS), len(AS)), np.nan)
for i, l in enumerate(LS):
    for j, a in enumerate(AS):
        rp = ratio_pair(float(l), float(a))
        if rp is not None:
            grid[i, j] = np.log(rp[0])
ok = ~np.isnan(grid)
print("three-direction, two-saddle points: %d of %d" % (ok.sum(), grid.size))
print("log r range: [%.4f, %.4f]" % (np.nanmin(grid), np.nanmax(grid)))
# sign changes along a
cnt = 0
for i in range(len(LS)):
    row = grid[i]
    for j in range(len(AS)-1):
        if not np.isnan(row[j]) and not np.isnan(row[j+1]) and row[j]*row[j+1] < 0:
            # bisect in a
            lo, hi = AS[j], AS[j+1]
            for _ in range(60):
                mid = 0.5*(lo+hi); v = ratio_pair(float(LS[i]), float(mid))
                if v is None: break
                if np.log(v[0])*np.log(ratio_pair(float(LS[i]), float(lo))[0]) < 0: hi = mid
                else: lo = mid
            am = 0.5*(lo+hi); rp = ratio_pair(float(LS[i]), float(am))
            cnt += 1
            if cnt <= 20:
                print("  r=1 at l=%7.2f a=%.9f   r=%.12f  u1=%.6f u2=%.6f  eta3=%+.6g"
                      % (LS[i], am, rp[0], rp[1], rp[2], eta3(LS[i], am)))
print("total r=1 crossings found on the scanned lines: %d" % cnt)
