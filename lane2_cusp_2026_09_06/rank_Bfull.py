#!/usr/bin/env python3
"""The independent-unfolding condition: rank(B_full) = 2, where

    B_full = [ d_mu D(s1) ;  d_mu D(s2) ]   over mu = (a, a20, a11, a01, a10).

If the two rows are nearly parallel, the displacement at the two stationary
points cannot be moved independently: one control raises/lowers BOTH, and no
parameter pair can drive both folds to zero at finite amplitude.  That is the
structural reason a solver escapes to the focus instead."""
import mpmath as mp
from engine import Engine
from cusp import Cusp
mp.mp.dps = 40
PARAMS = ("a", "a20", "a11", "a01", "a10")

ROWS = {
 8: dict(a="1.04", a20=-120, a11="1.51997", a01="1.56", a10="-79.6", s1="1.20", s2="3.48"),
 7: dict(a=mp.mpf(8)/11, a20=-12, a11="2.1502", a01=mp.mpf(67)/220, a10="-26.5", s1="1.28", s2="3.72"),
 2: dict(a="1.5", a20=-15, a11="0.79993", a01="3.2", a10="9.17", s1="1.12", s2="2.72"),
}
eng = Engine(); print("engine:", eng.banner)
for rid, P in ROWS.items():
    c = Cusp(eng, P["a"], P["a20"], side=1)
    mu = dict(a=mp.mpf(P["a"]), a20=mp.mpf(P["a20"]), a11=mp.mpf(P["a11"]),
              a01=mp.mpf(P["a01"]), a10=mp.mpf(P["a10"]))
    s1, s2 = mp.mpf(P["s1"]), mp.mpf(P["s2"])
    def val(m, s):
        c.a, c.a20 = m["a"], m["a20"]
        return c.val([m["a11"], m["a01"], m["a10"]], s)
    B = [[mp.mpf(0)]*5 for _ in range(2)]
    for j, nm in enumerate(PARAMS):
        h = mp.mpf("1e-9")*max(mp.mpf(1), abs(mu[nm]))
        mp_, mm_ = dict(mu), dict(mu); mp_[nm] += h; mm_[nm] -= h
        for i, s in enumerate((s1, s2)):
            rp, rm = val(mp_, s), val(mm_, s)
            if rp["status"] != "OK" or rm["status"] != "OK": B[i][j] = mp.nan; continue
            B[i][j] = (rp["D"] - rm["D"])/(2*h)
    print("\n=== row %d   s1=%s s2=%s ===" % (rid, mp.nstr(s1,6), mp.nstr(s2,6)))
    print("   d D(s1)/d mu = %s" % [mp.nstr(v, 5) for v in B[0]])
    print("   d D(s2)/d mu = %s" % [mp.nstr(v, 5) for v in B[1]])
    # normalise rows and measure the angle between them
    n0 = mp.sqrt(sum(v**2 for v in B[0])); n1 = mp.sqrt(sum(v**2 for v in B[1]))
    cosang = sum(B[0][j]*B[1][j] for j in range(5))/(n0*n1)
    M = mp.matrix(2, 5)
    for i in range(2):
        for j in range(5): M[i, j] = B[i][j]/(n0 if i == 0 else n1)
    G = M*M.T
    ev = mp.eigsy(mp.matrix([[G[0,0], G[0,1]],[G[1,0], G[1,1]]]), eigvals_only=True)
    sv = [mp.sqrt(abs(e)) for e in ev]
    print("   |row1| = %s   |row2| = %s" % (mp.nstr(n0, 5), mp.nstr(n1, 5)))
    print("   cos(angle between rows) = %s" % mp.nstr(cosang, 10))
    print("   singular values of the row-normalised B = %s" % [mp.nstr(v, 6) for v in sv])
    print("   effective rank-2 conditioning sigma_min/sigma_max = %s" %
          mp.nstr(min(sv)/max(sv), 6))
print("\ncalls:", eng.ncalls); eng.close()
