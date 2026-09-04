#!/usr/bin/env python3
"""Claude hostile checks of the Strike-2 exclusion ingredients.
C_a<9/4 (and its actual supremum), r(a)<-1/2, y bounds, the 601/136136
ratio bound at extreme lobe-box corners, the E0 integral bound, the
P0<=0 sign chain, and a brute-force attempt to violate Z(p1)<0 with
tau1<=5/11 by direct evaluation on lobe-region points."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "q4"))
import mpmath as mp
mp.mp.dps = 30
from fractions import Fraction as Q
def r_of_a(a):
    k = 1/(1-a); x = mp.asinh(mp.sqrt(k-1))
    O = mp.mpf(3)/10*mp.sinh(5*x/3)+mp.mpf(3)/2*mp.sinh(x/3)
    Ox = (mp.cosh(5*x/3)+mp.cosh(x/3))/2
    return -mp.sqrt(a)*Ox/(2*O)
Ca = lambda a: mp.mpf(3)/2*(1+a)+r_of_a(a)
grid = [mp.mpf(i)/1000 for i in range(1, 1000)]+[mp.mpf('0.999999'), mp.mpf('0.9999999999')]
sup = max(Ca(a) for a in grid); inf_r = max(r_of_a(a) for a in grid)
print("sup C_a on grid =", mp.nstr(sup, 12), " (claimed < 9/4 = 2.25; limit a->1 is 13/6 =", mp.nstr(mp.mpf(13)/6, 12), ")")
print("max r(a) on grid =", mp.nstr(inf_r, 12), " (claimed < -1/2)")
assert sup < mp.mpf(9)/4 and inf_r < -mp.mpf(1)/2
# y bounds (1-t)^{5/6} < y < (1-t)^{1/2}
def y_of(a, t):
    k = 1/(1-a); s = k-(k-1)*t
    O = lambda x: mp.mpf(3)/10*mp.sinh(5*x/3)+mp.mpf(3)/2*mp.sinh(x/3)
    return O(mp.asinh(mp.sqrt(s-1)))/O(mp.asinh(mp.sqrt(k-1)))
for a in (mp.mpf('0.05'), mp.mpf('0.5'), mp.mpf('0.95')):
    for t in (mp.mpf('0.1'), mp.mpf('0.5'), mp.mpf('0.9'), mp.mpf('0.999')):
        y = y_of(a, t); assert (1-t)**(mp.mpf(5)/6) < y < mp.sqrt(1-t), (a, t, y)
print("y bounds (1-t)^(5/6)<y<(1-t)^(1/2): OK")
# ratio bound m/eta < 601/136136 on the lobe bounding box corners (A>1+eta/6, B>-1, eta in (1,54/31))
worst = Q(0)
for eta in (Q(1)+Q(1,10**9), Q(54,31)-Q(1,10**9), Q(5,4)):
    for A in (1+eta/6+Q(1,10**9), Q(85,31)):
        for B in (Q(-1)+Q(1,10**9), Q(-49,744)):
            Y0 = 3*(1326*A+864*B-2431*eta-102)/1361360
            if Y0 < 0: worst = max(worst, -Y0/eta)
print("max |Y0|/eta over box corners =", float(worst), " bound 601/136136 =", float(Q(601,136136)))
assert worst < Q(601, 136136)
# E0 integral bound int_0^{5/11}(1-t)^{-13/6} = 6/7[(11/6)^{7/6}-1] < 8/9
val = mp.quad(lambda t: (1-t)**(-mp.mpf(13)/6), [0, mp.mpf(5)/11])
print("int_0^{5/11}(1-t)^{-13/6} =", mp.nstr(val, 12), "< 8/9:", val < mp.mpf(8)/9); assert val < mp.mpf(8)/9
# Direct hostile evaluation: Z(p1) on lobe-region points with tau1<=5/11 for several a
from q4_threshold_path import from_primitive_anchors_closed as from_primitive_anchors
from q4_reconstruction import reconstruct, green_coordinates, initial_weighted_derivative
from scipy.optimize import brentq
import numpy as np
mp.mp.dps = 30
best = -1e9
for anchors in ((0.1, 0.5, 0.9), (0.3, 0.6, 0.9), (0.45, 0.7, 0.95), (0.4545, 0.5, 0.6), (0.44, 0.9, 0.99)):
    co = tuple(map(float, from_primitive_anchors(tuple(map(mp.mpf, anchors)))))
    for a in (0.13, 0.3, 0.6, 0.9, 0.99):
        P0 = initial_weighted_derivative(a, *co)
        if P0 <= 0: continue
        sol = reconstruct(a, *co, t_end=0.995)
        P = lambda t: float(green_coordinates(a, sol, t)[1])
        ts = np.linspace(1e-4, anchors[0], 200); vals = [P(t) for t in ts]
        idx = [i for i in range(len(ts)-1) if vals[i]*vals[i+1] < 0]
        if not idx: continue
        p1 = brentq(P, ts[idx[0]], ts[idx[0]+1]); Z = float(green_coordinates(a, sol, p1)[0])
        best = max(best, Z)
        print(f"anchors={anchors} a={a} P0={P0:.3e} p1={p1:.4f} Z(p1)={Z:.6f}")
print("max Z(p1) found with tau1<=5/11 (must be <0):", best)
assert best < 0
print("STRIKE-2 INGREDIENT CHECKS PASSED")
