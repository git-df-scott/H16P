#!/usr/bin/env python3
"""Lane B: the two leading Dulac coefficients of the Q4 integral at the loop.
Near t=1: Y = Phi*y + P*y2 with y->0 like sqrt(1-t), y2(1) finite, and
P ~ -Omega0*H(1)*log(1-t). Hence
  c1 (coefficient of w log w in I)  is proportional to H(1)  [linear in (A,B,eta)],
  c0 = I at the loop = -(aC/2) sqrt(1-a) X(1), X(1)=int_0^1 Y/(1-au)^{3/2}  [affine].
Check: (i) c1=0 is the plane 9061A+6289B-2431eta=7242, which is a lobe-region
boundary (H(1)<0 strictly inside), so inside the lobe region the loop
cyclicity at first order is at most one; (ii) the map (A,B,eta)->(c0,c1) has
rank 2 at several kappa (loop functionals are independent); (iii) on the
c0=c1=0 line, how many interior zeros H and I have (numerical)."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "q4"))
import mpmath as mp, numpy as np
from claude_green_tools import *
from q4_reconstruction import reconstruct, original_values
mp.mp.dps = 25
for k in (2.0, 4.0, 9.0):
    a = 1-1/k
    def X1(co):
        sol = reconstruct(a, *map(float, co), t_end=0.999)
        # tail beyond 0.999: Y ~ log divergence, integrable; estimate crudely
        return sol.sol(0.999)[3]
    base = (1.2, -0.15, 1.2)
    grads = []
    for j in range(3):
        d = [0, 0, 0]; d[j] = 1e-4
        plus = X1([b+dd for b, dd in zip(base, d)]); minus = X1([b-dd for b, dd in zip(base, d)])
        grads.append((plus-minus)/2e-4)
    c1 = np.array([9061., 6289., -2431.])
    Mx = np.array([grads, c1])
    print(f"kappa={k}: grad X(1)={np.round(grads,6)}  rank of (c0,c1) map = {np.linalg.matrix_rank(Mx)}  angle cos={np.dot(grads,c1)/np.linalg.norm(grads)/np.linalg.norm(c1):.4f}")
    # solve c0=c1=0 with eta free: two linear equations in (A,B)
    X0 = X1(base) - np.dot(grads, base)
    for eta in (1.0, 1.2, 1.5):
        Mat = np.array([[grads[0], grads[1]], [9061., 6289.]])
        rhs = np.array([-X0-grads[2]*eta, 7242+2431*eta])
        A, B = np.linalg.solve(Mat, rhs)
        sol = reconstruct(a, A, B, eta, t_end=0.999)
        ts = np.linspace(0.02, 0.995, 400); I = original_values(a, sol, ts); H = sol.sol(ts)[0]
        zI = int(np.sum(I[:-1]*I[1:] < 0)); zH = int(np.sum(H[:-1]*H[1:] < 0))
        q0 = A-1-eta/6
        print(f"   c0=c1=0 line at eta={eta}: A={A:.5f} B={B:.5f}  q(0)={q0:.4f}  interior sign changes: H={zH} I={zI}")
