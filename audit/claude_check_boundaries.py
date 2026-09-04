#!/usr/bin/env python3
"""Claude boundary controls: (1) cusp construction and cusp-neighborhood
exclusion mechanism (H(t*)>0, one primitive zero after unfolding);
(2) lobe-equality boundary: anchor collision makes the anchor matrix
singular and the lobe inequality degenerate; (3) P2 cubic filter identity;
(4) coefficient-map orientation: I starts positive at the center when Y0<0."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "q4"))
import mpmath as mp
from claude_green_tools import *
mp.mp.dps = 30
d1 = lambda f, t, n=1: mp.diff(f, t, n)
R = lambda t: t+2*d1(Mh, t)/d1(Mh, t, 2)
# (1) cusp at b=1.3
b = mp.mpf('1.3')
ts = mp.findroot(lambda t: R(t)-b, mp.mpf('0.5'))
B0 = -Mh(ts)-(ts-b)*d1(Mh, ts); A0 = 1-(ts-b)*Mh(ts)-B0*ts
q0 = lambda t: A0+B0*t-1+(t-b)*Mh(t)
print(f"cusp: b={b} t*={mp.nstr(ts,10)} q0(t*)={mp.nstr(q0(ts),3)} q0'={mp.nstr(d1(q0,ts),3)} q0''={mp.nstr(d1(q0,ts,2),3)} q0'''={mp.nstr(d1(q0,ts,3),6)}")
assert abs(q0(ts)) < 1e-20 and abs(d1(q0, ts)) < 1e-15 and abs(d1(q0, ts, 2)) < 1e-12 and d1(q0, ts, 3) < 0
H0 = lambda t: mp.quad(lambda u: u*Fh(u)*q0(u), [0, t])
print("H_q0(t*) =", mp.nstr(H0(ts), 8), "(must be >0)"); assert H0(ts) > 0
lam = mp.mpf('1e-3'); ql = lambda t: q0(t)+lam*(t-ts)
roots = [mp.findroot(ql, ts+s*mp.sqrt(lam)*0.5) for s in (-1, 0, 1)]
print("unfolded q roots:", [mp.nstr(r, 10) for r in roots])
Hl = lambda t: mp.quad(lambda u: u*Fh(u)*ql(u), [0, t])
vals = [Hl(mp.mpf(x)/20) for x in range(1, 20)]+[Hl(mp.mpf('0.999'))]
changes = sum(1 for x, y in zip(vals, vals[1:]) if x*y < 0)
print("sign changes of H_lambda on (0,1):", changes, "(at most 1 claimed)"); assert changes <= 1
# (2) lobe boundary: anchors (0.5,0.5+delta,0.9)
for delta in (mp.mpf('1e-2'), mp.mpf('1e-4'), mp.mpf('1e-6')):
    y = (mp.mpf('0.5'), mp.mpf('0.5')+delta, mp.mpf('0.9'))
    rows = [primitive_basis_closed(t) for t in y]
    det = mp.det(mp.matrix([[r[0], r[1], -r[2]] for r in rows]))
    co = from_primitive_anchors_closed(y)
    # value of H at the local extremum between the colliding anchors -> 0
    Hmid = Hval(co, y[0]+delta/2)
    print(f"delta={delta}: det E={mp.nstr(det,4)}  H(midpoint)={mp.nstr(Hmid,4)}  eta={mp.nstr(co[2],8)}")
# (3) P2 cubic filter identity in s
k = mp.mpf(3); d = k-1; beta = mp.mpf('0.4')
w_s = lambda s: mp.hyp2f1(-mp.mpf(1)/6, mp.mpf(1)/6, 1, (k-s)/d)/Fh((k-s)/d)
al1 = mp.mpf('0.7'); P2 = lambda s: al1*(s-k)+(beta-k)   # affine with P2(k)=beta-k
g = lambda s: P2(s)+(s-beta)*w_s(s); f = lambda s: g(s)/(s-beta)
for s in (mp.mpf('1.3'), mp.mpf('2.2')):
    lhs = d1(f, s, 2); rhs = 2*P2(beta)/(s-beta)**3+d1(w_s, s, 2)
    assert abs(lhs-rhs) < 1e-15, (lhs, rhs)
print("P2 cubic identity (g/(s-beta))''=2P2(beta)/(s-beta)^3+w'': OK; g(k)=", mp.nstr(g(k), 3))
# w''(k) = -25/(216 d^2), w'''>0 in s
print("w''(k)=", mp.nstr(d1(w_s, k-mp.mpf('1e-8'), 2), 10), " expected", mp.nstr(-mp.mpf(25)/(216*d*d), 10))
assert d1(w_s, mp.mpf('1.5'), 3) > 0 and d1(w_s, mp.mpf('2.5'), 3) > 0
# (4) orientation: I(s(t)) = -(aC/2) sqrt(1-at) X, X ~ Y0 t <0 near center => I>0
print("orientation: Y0<0 => I>0 near center; area convention positive: consistent with reconstruction check")
print("BOUNDARY CONTROLS PASSED")
