#!/usr/bin/env python3
"""Is the trace independently controllable alongside the cusp conditions?

cusp5's Newton stalled with the residual sitting almost exactly on the T
equation.  A stalled Newton is not evidence, so measure the matrix instead:

    B = d(D, D_s, D_ss, T) / d(a, a20, a11, a01, a10)     (4 x 5, scaled)

If rank B = 4 the trace IS controllable at fixed cusp conditions and the stall
is a path problem; if the fourth singular value is at round-off, it is not.
"""
import json
import mpmath as mp
from engine import Engine
from cusp import Cusp
import foldcont as FC
mp.mp.dps = 40
C = json.load(open("ledger_opus/cusp_end.json"))
MU = {p: mp.mpf(C[p]) for p in FC.PARAMS}; SC = mp.mpf(C["s_c"])
eng = Engine(); print("engine:", eng.banner)
c = Cusp(eng, MU["a"], MU["a20"], side=1)
sc = FC.scales(MU)
g = FC.grad_mu(c, MU, SC, which=("D", "Dx", "Dxx"))
gT = [mp.mpf(v)*sc[j] for j, v in enumerate([-2, 0, 1, 1, 0])]
B = mp.matrix(4, 5)
for j in range(5):
    B[0,j] = g["D"][j]; B[1,j] = g["Dx"][j]; B[2,j] = g["Dxx"][j]; B[3,j] = gT[j]
names = ("d_mu D", "d_mu D_s", "d_mu D_ss", "grad_mu T")
print("\nB (scaled coordinates a, a20, a11, a01, a10):")
for i in range(4):
    print("   %-10s %s" % (names[i], [mp.nstr(B[i,j], 8) for j in range(5)]))
def svals(M, rn):
    A = mp.matrix(M.rows, M.cols)
    for i in range(M.rows):
        n = mp.sqrt(sum(M[i,j]**2 for j in range(M.cols))) if rn else mp.mpf(1)
        for j in range(M.cols): A[i,j] = M[i,j]/n
    G = A*A.T
    Gm = mp.matrix([[G[i,j] for j in range(G.cols)] for i in range(G.rows)])
    ev = mp.eigsy(Gm, eigvals_only=True)
    return sorted([mp.sqrt(max(mp.mpf(0), e)) for e in ev], reverse=True)
for rn, lab in ((True, "row-normalised"), (False, "raw")):
    s = svals(B, rn)
    print("\n   %s singular values of B: %s" % (lab, [mp.nstr(v,8) for v in s]))
    print("      sigma_4/sigma_1 = %s" % mp.nstr(s[-1]/s[0], 8))
# the 4x4 sub-block actually used by cusp5 (a, a11, a01, a10)  -- drop column a20
B4 = mp.matrix(4,4); cols = [0,2,3,4]
for i in range(4):
    for k, j in enumerate(cols): B4[i,k] = B[i,j]
print("\n   det of the 4x4 block (a, a11, a01, a10) = %s" % mp.nstr(mp.det(B4), 10))
s4 = svals(B4, True)
print("   row-normalised singular values: %s   ratio %s"
      % ([mp.nstr(v,8) for v in s4], mp.nstr(s4[-1]/s4[0], 8)))
# how much does a20 help?
print("\n   (a20 column, scaled) = %s" % [mp.nstr(B[i,1],8) for i in range(4)])
print("\ncalls:", eng.ncalls); eng.close()
