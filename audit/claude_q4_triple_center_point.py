#!/usr/bin/env python3
"""Step 6: the unique triple-center coefficient point (Y0=Y1=Y2=0) is
(A,B,eta)=(1,-17/12,0). Scan X(1;a); if it vanishes at some a*, unfold with
hierarchical small (Y0,eta,q0) of alternating signs and tune a near a* to
look for 3 small interior zeros plus one near the loop (four interior)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "q4"))
import numpy as np
from q4_reconstruction import reconstruct, original_values
from scipy.optimize import brentq
TE = 1-1e-7
ts = np.concatenate([np.logspace(-9, -1, 500), np.linspace(0.1, 0.9, 200), 1-np.logspace(-1, -6.9, 300)])
A0, B0, E0 = 1.0, -17/12, 0.0
def X1(a, co=(A0, B0, E0)): return reconstruct(a, *co, t_end=TE).sol(TE)[3]
print("X(1;a) at the triple-center point:")
avals = np.linspace(0.05, 0.995, 20); xv = [X1(a) for a in avals]
for a, v in zip(avals, xv): print(f"   a={a:.3f}: X(1)={v:+.4e}")
roots = [brentq(X1, avals[i], avals[i+1]) for i in range(len(avals)-1) if xv[i]*xv[i+1] < 0]
print("roots a* of X(1;a):", roots)
def count(a, co):
    sol = reconstruct(a, *co, t_end=TE); I = original_values(a, sol, ts); s = np.sign(I)
    idx = np.nonzero(s[:-1]*s[1:] < 0)[0]
    return len(idx), [float(ts[i]) for i in idx]
for astar in roots:
    print(f"=== unfolding near a*={astar:.8f} ===")
    # hierarchy: coefficients of X/t: Y0, Y1/2 = -eta/384 (at Y0~0), Y2/6 ~ -q0/13824 ; need alternation with the next (fixed) term
    for sgn in (+1, -1):
        for (e0, e1, e2) in ((1e-12, 1e-8, 1e-4), (1e-11, 1e-7, 1e-3), (1e-9, 1e-6, 1e-3), (1e-10, 1e-6, 1e-2)):
            # Y0 = sgn*e0, eta with Y1/2 = -eta/384 having sign -sgn => eta = sgn*384*e1*... choose eta = sgn*e1 (then Y1/2 = -sgn*e1/384)
            # q0 = A-1-eta/6 = -sgn*e2 gives Y2/6 ~ +sgn*e2/13824 (alternating again)
            eta = sgn*e1; A = 1+eta/6-sgn*e2
            B = (sgn*e0*1361360/3+2431*eta+102-1326*A)/864.0
            for da in (0.0, 1e-4, -1e-4, 1e-3, -1e-3, 1e-2, -1e-2):
                a = astar+da
                if not 0 < a < 1: continue
                n, locs = count(a, (A, B, eta))
                if n >= 3: print(f"   sgn={sgn:+d} (Y0,eta,q0)~({sgn*e0:+.0e},{eta:+.0e},{-sgn*e2:+.0e}) a={a:.6f}: interior zeros={n} at {locs}  X(1)={X1(a,(A,B,eta)):+.2e}")
    print("   (rows with fewer than 3 interior zeros suppressed)")
